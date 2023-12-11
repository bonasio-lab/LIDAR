#!/bin/bash

species=$1
threads=$2

# usage: ./index_setup.sh <species (hsap, mmus)> < # threads>

if command -v bowtie2 >/dev/null 2>&1 ; then
    : ; else
    echo "error -  bowtie2 not found" ;
    echo 'suggested conda environment to fix this issue: conda create -n LIDAR -c bioconda cutadapt trim-galore umi_tools bowtie2 samtools=1.16 gmap=2019.02.26 STAR=2.7.10a bbmap bedtools'
    echo 'then use conda activate LIDAR to activate the environment'
    exit 1 ; fi
if command -v STAR >/dev/null 2>&1 ; then
    : ; else
    echo "error - STAR not found"
    echo 'suggested conda environment to fix this issue: conda create -n LIDAR -c bioconda cutadapt trim-galore umi_tools bowtie2 samtools=1.16 gmap=2019.02.26 STAR=2.7.10a bbmap bedtools'
    echo 'then use conda activate LIDAR to activate the environment'
    exit 1 ; fi



gunzip ../annotation/*gz
gunzip ../index/$species/*gz
gunzip ../index/$species/tRNA_GSNAP/*gz
gunzip ../index/$species/tRNA_SNP/*gz

if [ ! -d ../index/$species/genome_STAR ]
then	
    mkdir ../index/$species/genome_STAR
fi

if [ ! -d ../index/$species/genome_bowtie ]
then	
    mkdir ../index/$species/genome_bowtie
fi

if [ ! -d ../index/$species/rRNA_STAR ]
then
    mkdir ../index/$species/rRNA_STAR
fi


if [ $species == "mmus" ]
then
    cat ../index/mmus/mm39.fa ../index/mmus/rRNA_tRNA.fa  > ../index/mmus/genome_rRNA_tRNA.fa
fi

if [ $species == "hsap" ]
then
    cat ../index/hsap/hg38.fa ../index/hsap/rRNA_tRNA.fa > ../index/hsap/genome_rRNA_tRNA.fa
fi



bedtools maskfasta -fi ../index/$species/genome_rRNA_tRNA.fa -bed ../index/$species/rRNArep_tRNA_mask.bed -fo ../index/$species/genome_rRNA_tRNA.mask.fa

STAR --runMode genomeGenerate --runThreadN $threads --genomeDir ../index/$species/genome_STAR/ --genomeFastaFiles ../index/$species/genome_rRNA_tRNA.mask.fa
STAR --runMode genomeGenerate --runThreadN $threads --genomeDir ../index/$species/rRNA_STAR/ --genomeFastaFiles ../index/$species/rRNA.fasta
bowtie2-build --threads $threads ../index/$species/genome_rRNA_tRNA.mask.fa ../index/$species/genome_bowtie/genome_rRNA_tRNA


rm ../index/$species/genome_rRNA_tRNA.mask.fa







    
