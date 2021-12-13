#!/bin/bash

## Script to create input for Artemis to visualize SNPs found between isolates in pairs
## Written by Jolinda de Korne (2021)
##Script needs:
## - directory with snippy output, called 'snippy_refFA1090_out'
## - directory 'pairwise_aln' with file 'comps.txt' including all comparisons
## - directory 'ann_vcf' with files 'annotations.hdr' & 'annotations.tab.gz' & 'annotations.tab.gz.tbi' & 'gene_list.txt'

CORE_ALN=snippy_refFA1090_out/clean.full.aln
LIN_CORE_ALN=snippy_refFA1090_out/clean.full.aln.lin
PAIRWISE_ALN=pairwise_aln
ANN_VCF=ann_vcf

##Make a linear alignment file to be able to extract isolates from the file
awk '/^>/ {printf("\n%s\n",$0);next; } { printf("%s",$0);}  END {printf("\n");}' < "$CORE_ALN" > "$LIN_CORE_ALN"

##Go to directory for pairwise alignments
cd "$PAIRWISE_ALN"

##Input file with all comparisons to make, called comps.txt. Extract for each comparison the core genome alignments from the involved isolates.
while read col1 col2 col3;
do
        grep -A 1 -e "$col1" -e "$col2" ../"$LIN_CORE_ALN" > "$col3"_"$col1"_"$col2".full.aln.fas ;
done < comps.txt

echo "Pairwise aligments were created"

##Run snp-sites to pairwisely compare isolates and get SNPs between isolates in a comparison in a  vcf file.
for file in *.full.aln.fas;
do
	echo "$file"
        snp-sites -v -o "$file".vcf "$file";
done

##Zip vcf files
for file in *.vcf;
do
        bgzip -f "$file";
done

##Create index file
for file in *.vcf.gz;
do
        tabix -p vcf "$file";
done

echo "VCF files were created, gzipped and indexed"

##Go to previous directory
cd ../"$ANN_VCF"
echo "This is the current directory ${PWD}"

##Annotate vcf files to create annotated vcfs as input for Artemis. Files annotations.tab.gz and annotations.hdr should be provided and originate from the reference genome FA1090 gff and bed files
for file in ../"$PAIRWISE_ALN"/*.vcf.gz;
do
	echo "$file will be annotated"
	bcftools annotate -a annotations.tab.gz -h annotations.hdr -c CHROM,FROM,TO,ID "$file" >> "$file".ann
	mv ../"$PAIRWISE_ALN"/*.ann ${PWD}
done

echo "VCF files were annotated"

##Search for each gene in the gene_list.txt in all files and get the number of samples that gene was annotated in (meaning these samples had SNPs in this gene). Gene_list.txt should be provided.
while read line in file;
do
	echo "$line"; grep -F "$line" *.ann | sort -k1,1 -u | wc -l ;
done < gene_list.txt > counts.txt

echo "SNPs in each gene were counted"
echo "Thanks for using this script!"
