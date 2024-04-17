#!/bin/bash

# Advanced Bioinformatics assignment
# Candidate ID: m2306450/k23165732
# Date: 16/04/2024

#Intial set-up of BAM files,  Post Alignment Q.C. and Filtering and Alignment Statistics


# Preliminary steps to make directories:

## Installing required tools for the setup for the analysis

          #Setting up the project structure
          # We make the working directory bioinformatics_assignment and make sub directories data meta results logs
mkdir bioinformatics_assignment
cd ~/bioinformatics_assignment
mkdir data meta results logs
cd ~/bioinformatics_assignment/data #moving to the data folder we make more directories for untrimmed and trimmed fastq data 
mkdir trimmed_fastq
mkdir untimmed_fastq

            ## download FASTQ raw read data
wget https://s3-eu-west-1.amazonaws.com/workshopdata2017/NGS0001.R1.fastq.qz
wget https://s3-eu-west-1.amazonaws.com/workshopdata2017/NGS0001.R2.fastq.qz
wget https://s3-eu-west-1.amazonaws.com/workshopdata2017/annotation.bed
mv *fastq.qz untrimmed_fastq
mv annotation.bed ~/bioinformatics_assignment/data

#we make a reference sub directory to save the reference file to map the data when performing the alignment
mkdir ../reference
cd ../reference
wget http://hgdownload.cse.ucsc.edu/goldenPath/hg19/bigZips/hg19.fa.gz

mv hg19.fa.gz ~/bioinformatics/data/

                # Installing all the pre-required libraries

cd ~/

                # Installing anaconda and its required packages

wget https://repo.anaconda.com/archive/Anaconda3-2022.10-Linux-x86_64.sh

chmod +x ./Anaconda3-2022.10-Linux-x86_64.sh

bash ./Anaconda3-2022.10-Linux-x86_64.sh

source ~/.bashrc

conda config --add channels defaults
conda config --add channels bioconda
conda config --add channels conda-forge

conda install samtools
conda install samtools
conda install bwa
conda install freebayes
conda install picard
conda install bedtools
conda install trimmomatic
conda install fastqc
conda install vcflib
conda install any_other_tool_you_need


# Create a simple README description of the project
touch README.md
echo "
Advanced Bioinformatics Course Assignment 2024\
Candidate number: m2306450/k23165732\
This project runs a standard Bioinformatics NGS pipeline to perform read \
alignment, variant discovery and annotation on some raw sequencing data." \
> README.md

                # To check the directory structure visually
# sudo apt-get install tree
tree bioinformatics_assignments #used this to separate data,logs,reference genomes, trimmed and untrimmed fastq data

                ##Navigate to the directory with the fastq files
cd ~/bioinformatics_assignment/data/untrimmed_fastq

#
fastqc -t 4 *.fastq.gz


mkdir ~/bioinformatics_assignment/results/fastqc_untrimmed_reads

mv *fastqc * ~/bioinformatics_assignment/results/fastqc_untrimmed_reads/

cd ~/bioinformatics_assignment/results/fastqc_untrimmed_reads/

                # unzipping processed untrimmed reads
        for zip in *.zip
        do
        unzip $zip
        done


        # **Exit v.m. to extract from home to extract html of the untrimmed fastq files

        #scp -i <path to .key file on your PC root@<ip address>:~/home/ubuntu/bioinformatics_assignment/data/untrimmed_fastq/*.html ~/Desktop



                #Trimming the raw sequence using Trimmomatic

cd ~/bioinformatics_assignment/data/untrimmed_fastq

# Performming the read trimming with Trimmomatic on raw sequencing data - to trim away adapters and filter out poor quality score reads
# Drops reads below 50 bp in length, and trims bases from the end of reads, if below a threshold quality (25).
# Also removes Illumina adapter sequences.

trimmomatic PE  \
-threads 4 \
-phred33 \
~/bioinformatics_assignment/data/untrimmed_fastq/NGS0001.R1.fastq.gz \ ~/bioinformatics_assignment/data/untrimmed_fastq/NGS0001.R2.fastq.gz \
-baseout ~/bioinformatics_assignment/data/trimmed_fastq/NGS0001.R1_trimmed_R \
  ILLUMINACLIP:/home/ubuntu/anaconda3/pkgs/trimmomatic-0.39-hdfd78af_2/share/trimmomatic-0.39-2/adapters/NexteraPE-PE.fa:2:30:10  \ TRAILING:25 MINLEN:50


                #Fastqc analysis *
# run FastQC to generate quality metrics on trimmed reads

fastqc -t 4 /home/ubuntu/bioinformatics_assignment/data/trimmed_fastq/trimmed_data_1P \
        /home/ubuntu/bioinformatics_assignment/data/trimmed_fastq/trimmed_data_2P

mkdir ~/bioinformatics_assignment/results/fastqc_trimmed_reads

mv ~/bioinformatics_assignment/data/trimmed_fastq/*fastqc * ~/bioinformatics_assignment/results/fastqc_trimmed_reads/


-- Alignment Duplication and Marking--

# Marking the index
mkdir -p ~/bioinformatics_assignment/data/reference

mv ~/bioinformatics_assignment/data/hg19.fa.gz ~/bioinformatics_assignment/data/reference/

 bwa index ~/bioinformatics_assignment/data/reference/hg19.fa.gz


# Running BWA mem with the RG information

mkdir ~/bioinformatics_assignment/data/aligned_data

 bwa mem -t 4 -v 1 -R '@RG\tID:HWI-D0011.50.H7AP8ADXX.1.WES01\tSM:WES01\tPL:ILLUMINA\tLB:nextera-wes01-blood\tDT:2017-02-23\tPU:HWI-D00119' -I 250,50  ~/bioinformatics_assignment/data/reference/hg19.fa.gz ~/bi>



       # Changing directories to aligned_folder

cd ~/bioinformatics_assignment/data/aligned_data


#Converting sam file to bam format, sort it and generate index using samtools


samtools view -h -b NGS0001.sam > NGS0001.bam

samtools sort NGS0001.bam > NGS0001_sorted.bam

#To generate .bai index file!
samtools index NGS0001_sorted.bam


                ## Post Alignment Q.C. and Filtering

workdir= ~/bioinformatics_assignment
alignment_dir= ~/bioinformatics_assignment/data/aligned_data #directory path has been shortened for easy calling for subsequenct analysis

        # this step locates and marks duplicated reads in the BAM file
picard MarkDuplicates I=$alignment_dir/NGS0001_sorted.bam O=$alignment_dir/NGS0001_sorted_marked.bam \
  M=$alignment_dir/marked_dup_metrics.txt
samtools index $alignment_dir/NGS0001_sorted_marked.bam # generate index
rm $alignment_dir/NGS0001_sorted.bam # remove uneeded data to save disk space

        # Filter reads below a minimum MAPQ score (-q 20). The samtools flag -F 1796 filters reads that are
        # unmapped, not in the primary alignment, that fail platform/vendor quality checks, or that are
        # PCR/optical duplicates. For more info see: https://broadinstitute.github.io/picard/explain-flags.html
samtools view -F 1796  -q 20 -o $alignment_dir/NGS0001_sorted_marked_filtered.bam \
  $alignment_dir/NGS0001_sorted_marked.bam

samtools index $alignment_dir/NGS0001_sorted_marked_filtered.bam # generate index



                ##Generate Alignment Statistics


mkdir $alignment_dir/alignment_stats # make a new directory for these stats
stats_dir=$alignment_dir/alignment_stats

                # Generate flagstats
samtools flagstats $alignment_dir/NGS0001_sorted_marked_filtered.bam \
  > $stats_dir/flagstats_output.txt

                # view the BAM file in the command line
# samtools view -h $alignment_dir/NGS0001_sorted_marked_filtered.bam | less

                # Generate alignment statistics per chromosome with idxstats
samtools idxstats  $alignment_dir/NGS0001_sorted_marked_filtered.bam \
  > $stats_dir/idxstats_output.txt

                # Determine the distribution of insert sizes between read pairs with Picard tools
picard CollectInsertSizeMetrics -H $stats_dir/CollectInsertSizeMetrics_histogram.pdf \
        -I $alignment_dir/NGS0001_sorted_marked_filtered.bam -O $stats_dir/CollectInsertSizeMetrics_output.txt

                # Calculate Depth of coverage.
                
# First get the all BAM regions that overlap with the input, then use this output to calculate the coverage.
#       (This saves memory usage)
bedtools intersect -bed -a NGS0001_sorted_marked_filtered.bam -b $workdir/data/annotation.bed \
  | bedtools coverage -d -a $workdir/data/annotation.bed -b - \
  > $stats_dir/coverageBed_output.txt

rm $alignment_dir/NGS0001_sorted_marked.bam # remove unneeded data to save disk space
