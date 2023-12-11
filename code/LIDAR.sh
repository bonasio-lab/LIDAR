#!/bin/bash

threads=1
flip=0
UMI="none"
UMI_preset="none"
del=1
cope_path="-1"
length_cutoff0=100
intron_max=100000
sample_name="sample"
mode="basic"
species="mmus"
R2_3prime_adapter="none"
tRNA_analysis=1
count_analysis=1
single=0

set -e
trap "trap - SIGTERM && kill -- -$$" SIGINT EXIT SIGTERM

code_dir=$( cd -- "$( dirname -- "${BASH_SOURCE[0]}" )" &> /dev/null && pwd )
index_dir=${code_dir}/../index/
annotation_dir=${code_dir}/../annotation/

current_dir=$(pwd)

start_time=$(date +%F.%H.%M)


while test $# -gt 0; do
    case "$1" in
	-h)
	    shift
	    echo
	    echo 'Mapping of fastq with/out UMIs'
	    echo
	    echo
	    echo 'Usage: ./LIDAR.sh [-h help] [-keep keep files (no)] [-flip flip orientation (no)] [-p threads (1)] [-c /path/to/cope if COPE not in path][--mode <basic, process, analyze_only> (basic)] [-species <mmus, hsap> (mmus)] [-no_tRNA skip tRNA analysis (no)] [-no_counts skip making count tables (no)] [-single single-end mode (no)] -samples sample_table'
	    echo
	    echo
	    echo 'OPTIONAL PARAMETERS:'
	    echo
	    echo '[-keep keep files (no)]                                        set -keep to keep temp files'
	    echo '[-flip flip orientation (no)]                                  set -flip if bam files need to be flipped so that reads are on the correct strand (Illumina pA RNA-seq)'
	    echo '[-p threads (1)]                                               number of threads'
	    echo '[-c /path/to/cope ()] - path to cope executable]               path to COPE, if not already in PATH and not placed in code folder'
	    echo '[-mode <basic, tRNA, process> (basic)]                         choose mode: basic = main analysis pipeline, process = only process fastqs, analyze_only = fastqs provided are already fully processed'
	    echo '[-species <hsap, mmus> (mmus)]                                 species (hsap = human hg38, mmus= mouse mm39), default is mouse'
	    echo '[-no_tRNA skip tRNA analysis (no)]                             set -no_tRNA if downstream tRNA analysis should be skipped'
	    echo '[-no_counts skip making count tables (no)]                     set -no_counts if downstream count table computation should be skipped'
	    echo '[-single single-end mode (no)]                                 set -single if running on single-end reads'
	    echo
	    echo 'REQUIRED PARAMETERS'
	    echo
	    echo '-samples                                                       tab separated sample table with columns for sample name, R1, R2, and library method (NEB, LIDAR, none)'
	    echo ''

	    echo 'SOFTWARE REQUIREMENTS:'
	    echo 'suggested conda environment: conda create -n LIDAR -c bioconda cutadapt trim-galore umi_tools bowtie2 samtools=1.16 gmap=2019.02.26 STAR=2.7.10a bbmap bedtools'
	    echo 'cope must be downloaded and location placed in PATH, path provided with -c <path/to/cope/, or the cope executable placed in the code folder'
	    echo ''
	    echo '-cope'
	    echo '-trim_galore'
	    echo '-bedtools'
	    echo '-samtools ver > 1.16'
	    echo '-umi_tools'
	    echo '-cutadapt'
	    echo '-STAR ver 2.7.10a'
	    echo '-GMAP ver 2019.02.26'
	    echo '-bowtie2'
	    echo '-bbmap'
	    echo ''
	    echo 'Required R packages:'
	    echo 'GenomicRanges'
	    echo 'GenomicAlignments'
	    echo 'rtracklayer'
	    echo 'Rsamtools'
	    echo 'parallel'
	    echo 'ggplot2'
	    echo 'Biostrings'
	    echo
	    shift
	    exit 1
	    ;;
	-mode)
	    shift
	    mode=$1
	    shift
	    ;;
	-p) shift
	    threads=$1
	    shift
	    ;;
	-flip)
	    shift
	    flip=1
	    ;;
	-c)
	    shift
	    cope_path=$1
	    PATH=$PATH:/$cope_path
	    shift
	    ;;
	-keep)
	    shift
	    del=0
	    ;;
	-samples)
	    shift
	    sample_table=$1
	    shift
	    ;;
	-species)
	    shift
	    species=$1
	    shift
	    ;;
	-no_tRNA)
	    shift
	    tRNA_analysis=0
	    ;;
	-no_counts)
	    shift
	    count_analysis=0
	    ;;
	-single)
	    shift
	    single=1
	    ;;
	*) echo "$1 is not a recognized flag! run with -h to see options"
	   exit 1;
	   ;;
    esac
done


if [ ! -f $sample_table ]
then
    echo "error - provide a sample table with -samples and  make sure this file exists"
fi


all_names=()
R1s=()
R2s=()
libs=()
single_reads=()

# add a newline to sample table if there isn't one
[[ $(tail -c1 $sample_table) && -f $sample_table ]]&&echo ''>>$sample_table



if [ $single == 0 ]
then
    while IFS=$'\t' read -r name R1 R2 lib single_read;do
	names+=($name)
	R1s+=($R1)
	R2s+=($R2)
	libs+=($lib)
	single_reads+=($single_read)
	if [ $lib != "LIDAR" ] && [ $lib != "NEB" ] && [ $lib != "none" ]
	then
	    echo "error - library types must be either LIDAR, NEB, or none, or you can provide processed fastqs with analyze_only mode"
	    exit 1
	fi
    done < $sample_table
    
    num_samples="${#names[@]}"
    if [ "${#names[@]}" != $num_samples ] || [  "${#R1s[@]}" != $num_samples ] || [ "${#R2s[@]}" != $num_samples ] || [ "${#libs[@]}" != $num_samples ]
    then
	echo "error - check your sample table, one of your columns may have fewer samples than the others, or you may have the wrong experiment type (single/paired)"
	exit 1
    fi
    
fi

																
if [ $single == 1 ]
then
    while IFS=$'\t' read -r name R1 lib;do
	names+=($name)
	R1s+=($R1)
	libs+=($lib)
	if [ $lib != "LIDAR" ] && [ $lib != "NEB" ] && [ $lib != "none" ] 
	then
	    echo "error - library types must be either LIDAR, NEB, or none, or you can provide processed fastqs with analyze_only mode. check your sample table"
	    exit 1
	fi
    done < $sample_table
    
    num_samples="${#names[@]}"

    if [ "${#names[@]}" != $num_samples ] || [  "${#R1s[@]}" != $num_samples ] || [ "${#libs[@]}" != $num_samples ]
    then
	echo "error - check your sample table, one of your columns may have fewer samples than the others"
	exit 1
    fi
    
fi

# check that all sample names are unique
len_unique=$(cut -f1 $sample_table | sort | uniq | wc -l)
len_all=$(cut -f1 $sample_table | wc -l)

if [ $len_unique != $len_all ]
then
    echo "error - in your sample table, two of the sames have the same names. all samples must have unique name"
    exit 1
fi

if [ $species != "mmus" ] && [ $species != "hsap" ]
then
    echo "error - species must be either mmus or hsap"
    exit 1
fi

if [ -f $code_dir/cope ]
then
    PATH=$PATH:$code_dir/
fi


if [ $mode != "analyze_only" ]
then    
    if command -v cope >/dev/null 2>&1 ; then
	: ; else
	echo "error - cope not found. install COPE from here https://sourceforge.net/projects/coperead/ and either put its directory your PATH, specify with -c /path/to/cope/, or place the cope executable in the code folder" ; exit 1 ; fi
    if command -v trim_galore >/dev/null 2>&1 ; then
	: ; else
	echo "error - trim_galore not found"
	echo 'suggested conda environment to fix this issue: conda create -n LIDAR -c bioconda cutadapt trim-galore umi_tools bowtie2 samtools=1.16 gmap=2019.02.26 STAR=2.7.10a bbmap'
	echo 'then use conda activate LIDAR to activate the environment'
	exit 1 ; fi
    if command -v cutadapt >/dev/null 2>&1 ; then
	: ; else
	echo "error! cutadapt not found"
	echo 'suggested conda environment to fix this issue: conda create -n LIDAR -c bioconda cutadapt trim-galore umi_tools bowtie2 samtools=1.16 gmap=2019.02.26 STAR=2.7.10a bbmap'
	echo 'then use conda activate LIDAR to activate the environment'
	    exit 1 ; fi
fi

if command -v bowtie2 >/dev/null 2>&1 ; then
    : ; else
    echo "error -  bowtie2 not found" ;
    echo 'suggested conda environment to fix this issue: conda create -n LIDAR -c bioconda cutadapt trim-galore umi_tools bowtie2 samtools=1.16 gmap=2019.02.26 STAR=2.7.10a bbmap'
    echo 'then use conda activate LIDAR to activate the environment'
    exit 1 ; fi
if command -v STAR >/dev/null 2>&1 ; then
    : ; else
    echo "error - STAR not found"
    echo 'suggested conda environment to fix this issue: conda create -n LIDAR -c bioconda cutadapt trim-galore umi_tools bowtie2 samtools=1.16 gmap=2019.02.26 STAR=2.7.10a bbmap'
    echo 'then use conda activate LIDAR to activate the environment'
    exit 1 ; fi

if command -v umi_tools >/dev/null 2>&1 ; then
    : ; else
    echo "error - umi_tools not found"
    echo 'suggested conda environment to fix this issue: conda create -n LIDAR -c bioconda cutadapt trim-galore umi_tools bowtie2 samtools=1.16 gmap=2019.02.26 STAR=2.7.10a bbmap'
    echo 'then use conda activate LIDAR to activate the environment'
    exit 1 ; fi


for ((c=0;c<=$num_samples;c++))
do
    if test -f ${R1s[$c]} ; then
	: ; else
	echo "error - one of your reads1 does not exist"
	echo "the file name was given as ${R1s[$c]} - check that this file is in the same folder as the sample table"
	exit 1 ; fi
    if [ $single == 0 ]
    then	
	if test -f ${R2s[$c]} ; then
	    : ; else
	    echo "error - one of your reads2 does not exist"
	    echo "the file name was given as ${R2s[$c]} - check that this file is in the same folder as the sample table"
	    exit 1 ; fi
    fi
done

if [ $count_analysis == 1 ] || [ $tRNA_analysis == 1 ]
then
    Rscript $code_dir/check_packages.R $current_dir
    if [ -f $current_dir/LIDAR.R.need_packages ]
    then
	echo "Error - for downstream analysis, one or more R packages are missing:"
	cat $current_dir/LIDAR.R.need_packages
	rm $current_dir/LIDAR.R.need_packages
	exit 1
    fi
fi




exec > LIDAR.log
exec 2>&1

date
echo "starting LIDAR analysis pipeline"


echo "parameters for all analysis:"

echo species: $species
echo delete files after running: $del
echo total threads provided: $threads
echo mode: $mode
echo cope path: $cope_path
echo invert bam files: $flip
echo single-end mode: $single

threads_per=$(($threads/$num_samples))

echo threads per sample: $threads_per

echo ""
echo ""

sequential=0
if [ "$threads" -lt "$num_samples" ]
then
    sequential=1
    echo "fewer threads provided than samples - running mapping jobs sequentially"
fi

num_samples=("${#names[@]}"-1)

date
if [ $sequential == 0 ]
then
    echo "starting mapping jobs"
    if [ $single == 0 ]
    then
	if [ $mode != "analyze_only" ]
	then
	    for ((c=0;c<=$num_samples;c++))
	    do
		$code_dir/LIDAR_mapping.sh -mode $mode -p $threads_per -flip2 $flip -c $cope_path -keep2 $del -1 ${R1s[$c]} -2 ${R2s[$c]} -name ${names[$c]} -species $species -lib ${libs[$c]} & done
	    wait
	fi
       
	if [ $mode == "analyze_only" ]
	then
	    for ((c=0;c<=$num_samples;c++))
	    do
		$code_dir/LIDAR_mapping.sh -mode $mode -p $threads_per -flip2 $flip -c $cope_path -keep2 $del -1 ${R1s[$c]} -2 ${R2s[$c]} -name ${names[$c]} -species $species -lib ${libs[$c]} -singles ${single_reads[$c]} & done
	    wait
	fi
    fi
    
	

    if [ $single == 1 ]
    then
	for ((c=0;c<=$num_samples;c++))
	do
	    $code_dir/LIDAR_mapping_SE.sh -mode $mode -p $threads_per -flip2 $flip -c $cope_path -keep2 $del -1 ${R1s[$c]} -name ${names[$c]} -species $species -lib ${libs[$c]} &
	    wait
	done
    fi
fi

	    
    
	
	    

if [ $sequential == 1 ]
then
    echo "starting mapping jobs"
    if [ $single == 0 ]
    then	
	for ((c=0;c<=$num_samples;c++))
	do
	    if [ $mode != "analyze_only" ]
	    then
		$code_dir/LIDAR_mapping.sh -mode $mode -p $threads -flip2 $flip -c $cope_path -keep2 $del -1 ${R1s[$c]} -2 ${R2s[$c]} -name ${names[$c]} -species $species -lib ${libs[$c]}
	    fi

	    if [ $mode == "analyze_only" ]
	    then
		echo "here!"
		$code_dir/LIDAR_mapping.sh -mode $mode -p $threads -flip2 $flip -c $cope_path -keep2 $del -1 ${R1s[$c]} -2 ${R2s[$c]} -singles ${single_reads[$c]} -name ${names[$c]} -species $species -lib ${libs[$c]}
	    fi    
	done
    fi

    if [ $single == 1 ]
    then
	for ((c=0;c<=$num_samples;c++))
	do
	    $code_dir/LIDAR_mapping_SE.sh -mode $mode -p $threads_per -flip2 $flip -c $cope_path -keep2 $del -1 ${R1s[$c]} -name ${names[$c]} -species $species -lib ${libs[$c]}
	done
    fi
fi


date


current_dir=$(pwd)
bam_dir=$current_dir/bam/
log_dir=$current_dir/logs/

for i in $bam_dir/*bam;do samtools index $i & done
wait

for i in $bam_dir/*bam;do samtools idxstats $i > ${i%.bam}.idxstats & done
wait

if [ ! -d $bam_dir/idxstats ]
then
    mkdir $bam_dir/idxstats
fi
mv $bam_dir/*.idxstats $bam_dir/idxstats/



if [ ! -d counts ] && [ $count_analysis == 1 ]
then
    mkdir counts
fi

if [ ! -d figures ] && ( [ $count_analysis == 1 ] || [ $tRNA_analysis == 1 ] )
then
    mkdir figures
fi


if [ $count_analysis == 1 ]
then
    date
    echo "finished with mapping all samples - starting count table analysis"
    echo "bam_dir" $bam_dir
    echo "annotation_dir" $annotation_dir
    echo "index_dir" $index_dir
    
    Rscript $code_dir/make_count_tables.R $bam_dir $annotation_dir $index_dir $threads $species $single
    
    echo "done with count table analysis"
fi

if [ $tRNA_analysis == 1 ]
then
    date
   echo "running tRNA analysis:"
    if [ ! -d $bam_dir/tRNA/ ]
    then
	mkdir $bam_dir/tRNA
    fi
	
    for i in $bam_dir/*bam;do samtools view -b -L $annotation_dir/tRNA/$species.tRNA.bed $i > ${i%.bam}.tRNA.bam & done
    wait
    mv $bam_dir/*tRNA.bam $bam_dir/tRNA/

    Rscript $code_dir/tRNA_analysis.R $bam_dir $annotation_dir $species $single
fi

echo "done with all analysis"
date

mkdir LIDAR_analysis_$start_time
mv LIDAR.log bam logs LIDAR_analysis_$start_time

if [ $del == 0 ]
then
    mv *MAP LIDAR_analysis_$start_time
fi


if [ -d figures ];then mv figures LIDAR_analysis_$start_time;fi
if [ -d counts ];then mv counts LIDAR_analysis_$start_time;fi


trap  - SIGINT EXIT SIGTERM
