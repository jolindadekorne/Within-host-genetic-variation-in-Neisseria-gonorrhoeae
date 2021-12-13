IDS, = glob_wildcards("raw_data/{id}_fastpout_1.fastq.gz")

rule all:
        input:
                expand("fastp_out/{sample}_fastp.json", sample = IDS),
                "multiqc_fastp_out",
                expand("skesa_out/{sample}.fasta", sample = IDS),
                expand("quast_out/{sample}", sample=IDS),
                "multiqc_quast_out",
		expand("coverage_FA1090/{sample}_cov.txt", sample = IDS),
		expand("snippy_refFA1090_out/{sample}", sample = IDS),
		"snippy_refFA1090_out/clean.full.aln",
                "maskrc_gubbins_FA1090_out/masked_gubbins_snippyrefFA1090.aln"

rule fastp:
        input:
                fw = "raw_data/{sample}_R1.fastq.gz",
                rv = "raw_data/{sample}_R2.fastq.gz"
        output:
                fw = "trimmed_illumina/{sample}_fastpout_1.fastq.gz",
                rv = "trimmed_illumina/{sample}_fastpout_2.fastq.gz",
                json = "fastp_out/{sample}_fastp.json",
                html = "fastp_out/{sample}_fastp.html"
	conda: 
		"envs/fastp.yml"
	params:
		general = "--disable_length_filtering",
		compression_level = 9
	log:
		"logs/fastp/fastp_{sample}.log"
	threads: 4
	shell:
		"""
		fastp -w {threads} -z {params.compression_level} -i {input.fw} -o {output.fw} -I {input.rv} -O {output.rv} {params.general} --html {output.html} --json {output.json} 2>&1>{log}
		"""

rule multiqc_fastp:
	input:
		expand("fastp_out/{sample}_fastp.json", sample=IDS)
	output:
		directory("multiqc_fastp_out")
	conda: 
		"envs/multiqc.yml"
	log:
		"logs/multiqc/multiqc_fastp.log"
	shell:
		"""
		mkdir -p {output}
		multiqc {input} --outdir {output} 2>&1>{log}
		"""

rule skesa:
	input:
		fw = "trimmed_illumina/{sample}_fastpout_1.fastq.gz",
		rv = "trimmed_illumina/{sample}_fastpout_2.fastq.gz"
	output:
		assembly = "skesa_out/{sample}.fasta"
	params:
		min_length = "500"
	conda: 
		"envs/skesa.yml"
	log:
		"logs/skesa/skesa_{sample}.log"
	threads: 6
	shell:
		"""
		skesa --fastq {input.fw},{input.rv} --use_paired_ends --min_contig {params.min_length} 1> {output.assembly} 2>{log}
		"""

rule quast:
	input:
		assembly = "skesa_out/{sample}.fasta"
	output:
		directory("quast_out/{sample}")
	conda:
		"envs/quast.yml"
	log:
		"logs/quast/quast_{sample}.log"
	shell:
		"""
		quast.py -o {output} {input.assembly}
		"""

rule multiqc_quast:
	input:
		expand("quast_out/{sample}", sample=IDS)
	output:
		directory("multiqc_quast_out")
	params:
		input = expand("quast_out/{sample}/report.tsv", sample=IDS)
	conda:
		"envs/multiqc.yml"
	log:
		"logs/multiqc_quast/multiqc_quast.log"
	shell:
		"""
		multiqc {params.input} -o {output} 2>&1>{log}
		"""

rule coverage:
	input:
		fw = "trimmed_illumina/{sample}_fastpout_1.fastq.gz",
		rv = "trimmed_illumina/{sample}_fastpout_2.fastq.gz"
	output:
		cov = "coverage_FA1090/{sample}_cov.txt"
	log:
		bam = "logs/bwamem_samtools/{sample}.log",
		cov = "logs/coverage_FA1090/coverage_FA1090_{sample}.log"
	conda:
		"envs/samtools.yml"
	params:
		ref = "NgRefFA1090",
		path = "/home/jdkorne/samtools_1.11/bin",
		bam_out = "temp_bam/{sample}_temp_sorted.bam"
	threads: 6
	shell:
		"""
		mkdir -p temp_bam
		mkdir -p coverage_FA1090
		/home/jdkorne/bwa-mem2/bwa-mem2 index -p {params.ref} {params.ref}.fna
		/home/jdkorne/bwa-mem2/bwa-mem2 mem -t {threads} {params.ref} {input.fw} {input.rv} | {params.path}/samtools sort --threads {threads} -o {params.bam_out} 2>&1>{log.bam}
		{params.path}/samtools coverage -o {output.cov} {params.bam_out} 2>&1>{log.cov}
		rm {params.bam_out}
		"""

rule snippy_FA1090:
	input:
		fw = "trimmed_illumina/{sample}_fastpout_1.fastq.gz",
		rv = "trimmed_illumina/{sample}_fastpout_2.fastq.gz"
	output:
		directory("snippy_refFA1090_out/{sample}")
	conda:	
		"envs/snippy.yml"
	params:
		outdir = directory("snippy_refFA1090_out"),
		ref = "NgRefFA1090.gbk"
	log:    "logs/snippy_refFA1090/snippy_refFA1090_{sample}.log"
	threads: 16
	shell:
		"""
		mkdir -p {params.outdir}
		snippy --cpus {threads} --ref {params.ref} --R1 {input.fw} --R2 {input.rv} --outdir {output} 2>&1>{log}
		"""

rule snippy_aln_FA1090:
	input:
		corein = expand("snippy_refFA1090_out/{sample}", sample = IDS)
	output:
		cleanout = "snippy_refFA1090_out/clean.full.aln"
	conda:	
		"envs/snippy.yml"
	params:
		outdir = directory("snippy_refFA1090_out"),
		ref = "NgRefFA1090.gbk",
		prefix = "core",
		cleanin = "core.full.aln"
	shell:
		"""
		snippy-core {input.corein} --ref {params.ref} --prefix {params.prefix}
		snippy-clean_full_aln {params.cleanin} > {output.cleanout}
		"""

rule gubbins:
	input:
		snippy_FA1090 = "snippy_refFA1090_out/clean.full.aln"
	output:
		FA1090 = "maskrc_gubbins_FA1090_out/masked_gubbins_snippyrefFA1090.aln"
	conda:
		"envs/gubbins.yml"
	params:	
		gubbins_out_FA1090 = "gubbins_snippyrefFA1090_out",
		maskrc_out_FA1090 = "maskrc_gubbins_FA1090_out",
		model = "GTRGAMMA",
		prefix_FA1090 = "gubbins_snippyrefFA1090"
	log:	gubbins_FA1090 = "logs/gubbins/gubbins_snippyrefFA1090.log",
		maskrc_FA1090 = "logs/maskrc_gubbins_snippyrefFA1090.log"
	shell:
		"""
		mkdir -p {params.gubbins_out_FA1090} {params.maskrc_out_FA1090} 
		run_gubbins.py {input.snippy_FA1090} --prefix {params.prefix_FA1090} --raxml_model {params.model} 2>&1>{log.gubbins_FA1090}
		python3 scripts/maskrc-svg.py --gubbins --aln {input.snippy_FA1090} --out {output.FA1090} {params.prefix_FA1090} 2>{log.maskrc_FA1090}
	        mv {params.prefix_FA1090}.* {params.gubbins_out_FA1090}
		"""

rule snp_dists_snippy:
	input:
		snippy_FA1090 = "snippy_refFA1090_out/clean.full.aln",
		masked_gubbins_FA1090 = "maskrc_gubbins_FA1090_out/masked_gubbins_snippyrefFA1090.aln"
	output:
		dists_snippy_FA1090 = "snp_dists/snp_dists_snippy_FA1090.tsv",
		dists_masked_gubbins_FA1090 = "snp_dists/snp_dists_masked_gubbins_FA1090.tsv"
	params:
		dir = "snp_dists"
	conda:
		"envs/snp_dists.yml"
	shell:
		"""
		mkdir -p {params.dir}
		snp-dists -m {input.snippy_FA1090}>{output.dists_snippy_FA1090} 
		snp-dists -m {input.masked_gubbins_FA1090}>{output.dists_masked_gubbins_FA1090}
		"""

