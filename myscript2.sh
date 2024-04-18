#!/usr/bin/env bash

echo "We'll call the previous pipeline here."

sh ~/myscript2.shecho "We'll call the previous pipeline here."
sh ~/myscript2.sh


# we are going to run the previous myscript.sh here, hence we do not need to run this in the terminal

bash pipeline.sh ~/bioinformatics_assignment/data/untrimmed_fastq/NGS0001.fastq.gz \
        ~/bioinformatics_assignment/data/untrimmed_fastq/NGS0002.fastq.gz

cd ~/bioinformatics_assignment/results/fastqc_trimmed_reads/

ls
trimmed_data_1P_fastqc.html  trimmed_data_1P_fastqc.zip  trimmed_data_2P_fastqc.html  trimmed_data_2P_fastqc.zip

# We are going to convert to text format the reference (as required bysamtools faidx), index it with samtools faidx, call variants with Freebayes, compress the resulting variant file (VCF) and index the VCF with tabix:

zcat ~/bioinformatics_assignment/data/reference/hg19.fa.gz > ~/bioinformatics_assignment/data/reference/hg19.fa

samtools faidx ~/bioinformatics_assignment/data/reference/hg19.fa

freebayes --bam ~/bioinformatics_assignment/data/aligned_data/NGS0001_sorted_filtered.bam --fasta-reference ~/bioinformatics_assignment/data/reference/hg19.fa --vcf ~/bioinformatics_assignment/results/NGS0001.>

bgzip ~/bioinformatics_assignment/results/NGS0001.vcf

tabix -p vcf ~/bioinformatics_assignment/results/NGS0001.vcf.gz

                # Filtering the VCF using Freebayes

conda install vcflib #if you did not install vcflib before. vcffilter is part of the vcflib suite

vcffilter -f "QUAL > 1 & QUAL / AO > 10 & SAF > 0 & SAR > 0 & RPR > 1 & RPL > 1" \ ~/bioinformatics_assignment/results/NGS0001.vcf.gz > ~/bioinformatics_assignment/results/NGS0001_filtered.vcf

# As we run Freebayes with default parameters, the resulting VCF contains a large number of bad calls. The Freebayes Information Fields we will use for filtering are:

# QUAL=probability that there is a polymorphism at the loci described by the record: 1 - P(locus is homozygous given the data). 
#AO=Alternate allele observations, with partial observations recorded fractionally AF=Number of alternate observations on the forward strand SAR=Number of alternate observations on the reverse strand RPL=Reads 
#Placed Left: number of reads supporting the alternate balanced to the left (5’) of -the alternate allele RPR=Reads 
#Placed Right: number of reads supporting the alternate balanced to the right (3’) of the alternate allele 
#What calls do we “know” are poor (w.r.t. freebayes VCF annotations)?

# low-quality calls (QUAL < N)
# also, maybe QUAL is high but QUAL / AO is low
# loci with low read depth (DP < N)
# alleles that are only seen on one strand
# remove by “SAF > 0 & SAR > 0”
# alleles that are only observed by reads placed to the left or right
# remove by “RPL > 0 & RPR > 0”

# We will apply the following suggested freebayes hard filter for human diploid sequencing and use vcffilter to apply it:

# QUAL > 1 & QUAL / AO > 10 & SAF > 0 & SAR > 0 & RPR > 1 & RPL > 1

# QUAL > 1: removes horrible sites QUAL / AO > 10 : additional contribution of each obs should be 10 log units (~ Q10 per read) SAF > 0 & SAR > 0 : reads on both strands RPR > 1 & RPL > 1 : at least two reads "balanced" to each side of the site



                # Filtering the VCF!!!

bedtools intersect -header -wa -a ~/bioinformatics_assignment/results/NGS0001_filtered.vcf ../annotation.bed \
        > ~/bioinformatics_assignment/results/NGS0001_filtered.vcf_annotation.vcf

 bgzip ~/bioinformatics_assignment/results/NGS0001_filtered.vcf_annotation.vcf

 tabix -p vcf ~/bioinformatics_assignment/results/NGS0001_filtered.vcf_annotation.vcf.gz


                # ANNOTATION
# Download and setup Annovar database for variant annotation
mkdir ~/annovar
cd ~/annovar
tar -zxvf annovar.latest.tar.gz # N.B. this file must be downloaded from the Annovar website as it is not opensource
# Application form link = http://download.openbioinformatics.org/annovar_download_form.php

# setup databases
./annotate_variation.pl -buildver hg19 -downdb -webfrom annovar knownGene humandb/
./annotate_variation.pl -buildver hg19 -downdb -webfrom annovar refGene humandb/
./annotate_variation.pl -buildver hg19 -downdb -webfrom annovar ensGene humandb/
./annotate_variation.pl -buildver hg19 -downdb -webfrom annovar clinvar_20180603 humandb/
./annotate_variation.pl -buildver hg19 -downdb -webfrom annovar exac03 humandb/
./annotate_variation.pl -buildver hg19 -downdb -webfrom annovar dbnsfp31a_interpro humandb/
# ./annotate_variation.pl -buildver hg19 -downdb -webfrom annovar snp138 humandb/
# due to the size of the dbSNP 138 database, I had to download and run this step in my local environment

cd ~/

# So now we will annotate and prioritize the variant

# Convert VCF to annovar input
~/tools/annovar/convert2annovar.pl -format vcf4 $workdir/results/NGS0001_filtered_in_bedfile.vcf.gz  \
  > $workdir/results/NGS0001_filtered_in_bedfile.avinput

# Run Annovar to annotate variants with database frequencies/functional consequences.
#       Output is a .csv that can be opened in MS Excel.
~/tools/annovar/table_annovar.pl $workdir/results/NGS0001_filtered_in_bedfile.avinput ~/tools/annovar/humandb/ -buildver hg19 \
  -out $workdir/results/NGS0001_filtered_in_bedfile -remove \
  -protocol refGene,ensGene,clinvar_20180603,exac03,dbnsfp31a_interpro -operation g,g,f,f,f -otherinfo -nastring . -csvout

# # Line for including dbSNP 138 database (due to size of database we've commented it out)
# ~/tools/annovar/table_annovar.pl $workdir/results/NGS0001_filtered_in_bedfile.avinput ~/tools/annovar/humandb/ -buildver hg19 \
#   -out $workdir/results/NGS0001_filtered_in_bedfile -remove \
#   -protocol refGene,ensGene,clinvar_20180603,exac03,dbnsfp31a_interpro,snp138 -operation g,g,f,f,f,f -otherinfo -nastring . -csvout
