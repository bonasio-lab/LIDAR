#!/bin/bash

# process fastq files with/out UMI and map to the genome
# requirements: STAR, bowtie2, COPE, TrimGalore, UMI_tools, cutadapt

end_type=$1
ncRNA_genome_dir=$2
ncRNA=$3
sample_name=$4
threads=$5

if [ $end_type == "paired" ]
then
	reads1=$6
	reads2=$7
	
	bowtie2 --very-sensitive -x ${ncRNA_genome_dir}/$ncRNA -1 $reads1 -2 $reads2 > ${sample_name}.map.sam
	
	# keep any alignments that have insert size >16, mapped in proper pair, only primary alignments
	awk 'substr($0,1,1)=="@" || ($9>15) || ($9<-15)' ${sample_name}.map.sam | samtools view -@ $threads -F 260 -f 2 -b | samtools sort -o ${sample_name}.paired.${ncRNA}.bam
	
	# extract any alignments with insert size <= 15
	awk 'substr($0,1,1)=="@" || ($9<0 && $9>=-15) || ($9>0 && $9<=15)' ${sample_name}.map.sam | samtools view -@ $threads -b | samtools sort -o ${sample_name}.short.bam -
	
	# any reads unmapped or mapped with discordant pairs, etc
	samtools view -@ $threads -F 2 -b ${sample_name}.map.sam | samtools sort -o ${sample_name}.unmapped.bam - 
	
	# merge the short and unmapped reads
	samtools merge -@ $threads ${sample_name}.remap.unsorted.bam ${sample_name}.unmapped.bam ${sample_name}.short.bam
	samtools sort -@ $threads -n -o ${sample_name}.remap.bam ${sample_name}.remap.unsorted.bam
	
	samtools bam2fq -1 ${sample_name}.paired.${ncRNA}.remap.R1.gz -2 ${sample_name}.paired.${ncRNA}.remap.R2.gz -s ${sample_name}.singleton.gz -n ${sample_name}.remap.bam
	
	
	rm ${sample_name}.map.sam ${sample_name}.unmapped.bam ${sample_name}.short.bam ${sample_name}.remap.unsorted.bam ${sample_name}.remap.bam ${sample_name}.singleton.gz
fi


if [ $end_type == "single" ]
then
	reads=$6
	
	bowtie2 --very-sensitive -x ${ncRNA_genome_dir}/$ncRNA -U $reads -S ${sample_name}.map.sam

	# keep any alignments that are > 15, only primary alignments
	samtools view -F 4 -e 'length(seq)>15' -b ${sample_name}.map.sam | samtools sort -o ${sample_name}.single.${ncRNA}.bam
	
	# extract any alignments that are < 15
	samtools view -F 4 -e 'length(seq)<=15' -b ${sample_name}.map.sam | samtools sort -o ${sample_name}.short.bam
	
	# extract any reads that don't map
	samtools view -f 4 -b ${sample_name}.map.sam | samtools sort -o ${sample_name}.unmapped.bam
	
	# merge the short and unmapped reads
	samtools merge -@ $threads ${sample_name}.remap.unsorted.bam ${sample_name}.unmapped.bam ${sample_name}.short.bam
	samtools sort -@ $threads -n -o ${sample_name}.remap.bam ${sample_name}.remap.unsorted.bam
	
	samtools bam2fq -n ${sample_name}.remap.bam > ${sample_name}.single.${ncRNA}.remap.single
	
	rm ${sample_name}.map.sam ${sample_name}.unmapped.bam ${sample_name}.short.bam ${sample_name}.remap.unsorted.bam ${sample_name}.remap.bam 
	
fi
