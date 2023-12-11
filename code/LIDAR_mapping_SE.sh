#!/bin/bash

# process single-end fastq files with/out UMI and map to the genome

threads=1
flip=0
UMI="none"
del=1
cope_path="-1"
length_cutoff0=100
intron_max=100000
sample_name="sample"
mode="basic"
species="mmus"


code_dir=$( cd -- "$( dirname -- "${BASH_SOURCE[0]}" )" &> /dev/null && pwd )
index_dir=${code_dir}/../index/

while test $# -gt 0; do
    case "$1" in
	-h)
	   shift
	   echo
	   echo 'Mapping of fastq with/out UMIs'
	   echo
	   echo
	   echo 'Usage: ./LIDAR_ver3.sh [-h help] [-keep keep files (no)] [-flip flip orientation (no)] [-p threads (1)] [-c /path/to/cope if COPE not in path] [-lib lib_type <LIDAR, 8bp, NEB> (LIDAR)] [-name sample_name (sample)] [--mode <basic, process> (basic)] [-species <mmus, hsap> (mmus)] -1 R1_reads'
	   echo
	   echo
	   echo 'OPTIONAL PARAMETERS:'
	   echo
	   echo '[-keep keep files (no)]                                        set -k to keep temp files'
	   echo '[-flip flip orientation (no)]                                  set -f if bam files need to be flipped so that reads are on the correct strand (Illumina pA RNA-seq)'
	   echo '[-p threads (1)]                                               number of threads'
	   echo '[-c /path/to/cope ()] - path to cope executable]               path to COPE, if not already in PATH'              
	   echo '[-lib (LIDAR)]                  lib type (LIDAR - LIDAR UMI, 8bp - straight 8bp UMI, NEB - no UMI, other - no UMI)'
	   echo '[-name sample name (sample)]                                   base name for sample'
	   echo '[-mode <basic, tRNA, process> (basic)]                         choose mode: basic = main analysis pipeline, process = only process fastqs'
	   echo '[-species <hsap, mmus> (mmus)]                                 species (hsap = human hg38, mmus= mouse mm39), default is mouse'
	   echo
	   echo 'REQUIRED PARAMETERS'
	   echo
	   echo '-1                                                             R1 reads'
	   echo ''
 
	   echo 'SOFTWARE REQUIREMENTS:'
	   echo 'suggested conda environment: conda create -n LIDAR -c bioconda cutadapt trim-galore umi_tools bowtie2 samtools=1.16 gmap=2019.02.26 STAR=2.7.10a bbmap'
	   echo 'cope must be downloaded and placed in PATH, or path provided with -c <path/to/cope/'
	   echo ''
	   echo '-cope'
	   echo '-trim_galore'
	   echo '-samtools ver > 1.16'
	   echo '-umi_tools'
	   echo '-cutadapt'
	   echo '-STAR ver 2.7.10a'
	   echo '-GSNAP ver 2019.02.26'
	   echo '-bowtie2'
	   echo '-cutadapt'
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
	-flip2)
	    shift
	    flip=$1
	    shift
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
	-keep2)
	    shift
	    del=$1
	    shift
	    ;;
	-1)
	    shift
	    reads1=$1
	    shift
	   ;;
	-name)
	    shift
	    sample_name=$1
	    shift
	   ;;
	-species)
	    shift
	    species=$1
	    shift
	    ;;
	-lib)
	    shift
	    lib=$1
	    shift
	   ;;
	*) echo "$1 is not a recognized flag! run with -h to see options"
	   exit 1;
	   ;;
    esac
done


if [ $lib == "NEB" ]
then
    UMI="none"
fi

if [ $lib == "other" ]
then
    UMI="none"
fi


if [ $lib == "LIDAR" ]
then
    UMI="LIDAR"
fi

if [ $lib == "8bp" ]
then
    UMI="8bp"
fi


if [[ $name == NEB* ]] && [ $UMI != "none" ]
then
    echo "warning! you are running an NEB sample with UMI set!"
fi

if test -f "$reads1" ; then
    : ; else
    echo "error - reads1 does not exist!"
    echo "the file name was given as $reads1 - check that this file exists"
    exit 1 ; fi



shift $(($OPTIND-1))


genome_bowtie=$index_dir/$species/genome_bowtie/$species.LIDAR
genome_STAR=$index_dir/$species/genome_STAR/
rRNA_STAR=$index_dir/$species/rRNA_STAR/
ncRNA_bowtie=$index_dir/$species/ncRNA_bowtie/
tRNA_GSNAP_D=$index_dir/$species/tRNA_GSNAP/
tRNA_SNP_V=$index_dir/$species/tRNA_SNP/


if [ $species == "mmus" ]
then
    tRNA_SNP_v="mm39_modificationSNPs"
    tRNA_GSNAP_d="Mmus_tRNAgenome"
fi

if [ $species == "hsap" ]
then
    tRNA_SNP_v="hg38_modificationSNPs"
    tRNA_GSNAP_d="Hsap_tRNAgenome"
fi


length_cutoff=$((length_cutoff0))

if [ ! -d $sample_name.MAP ]
then
    mkdir $sample_name.MAP
fi


exec > ${sample_name}.bash.log
exec 2>&1
date

echo "parameters:"

echo reads1 file: $reads1
echo same name: $sample_name
echo library type: $lib
echo UMI type from library: $UMI
echo species: $species
echo delete files after running: $del
echo threads for this sample: $threads
echo mode: $mode
echo cope path: $cope_path
echo invert bam file: $flip




#######################
#### PROCESS READS ####
#######################

br='


##########################################'
br2='#########################################

'


echo "$br"
echo 'trimming with trim_galore'
echo "$br2"

# trim_galore: trimming of auto-detected adapters. new version: keep reads where the r2 gets trimmed because of length, but R1 is long enough to contain the UMI and an insert (set to 20 - any short reads will be removed later anyway). don't want to keep any R2 only reads, since they won't have UMI, so set the r2 length threshold very high. 
trim_galore -j $threads --length 5 -o $sample_name.MAP --basename ${sample_name} --dont_gzip $reads1 

cd $sample_name.MAP


if [ $UMI != "none" ]
then
    echo "$br"
    echo 'extracting UMIs and trimming remainder of UMI construct, if needed'
    echo "$br2"

    $code_dir/process_files_single.sh ${sample_name}_trimmed.fq ${sample_name} $UMI
  
    # stats on UMI files
    
    grep "@" ${sample_name}.UMI | wc -l > ${sample_name}.03a.UMI.single.count
    sed -i "s/^/$sample_name\tUMI\tsingle\tall\t/" ${sample_name}.03a.UMI.single.count
    
    
    # clean up of big files
    rm ${sample_name}.UMI
    rm ${sample_name}.to_map0.single
       
fi

# do some stats and then clean up - large files that we don't need after doing stats

zcat ../$reads1 | grep "@" | wc -l > ${sample_name}.00.original.count
sed -i "s/^/$sample_name\toriginal\tall\tall\t/" ${sample_name}.00.original.count

trimmed_count="$(zcat ${sample_name}_trimmed.fq | grep "@" | wc -l | awk '{print $1}')"


echo $((trimmed_count)) > ${sample_name}.01.trimmed.count
sed -i "s/^/${sample_name}\ttrimmed\tall\tall\t/" ${sample_name}.01.trimmed.count


if [ $UMI == "none" ]
then
    mv ${sample_name}_trimmed.fq ${sample_name}.to_map.single
fi

# stats on to map files
grep "@" ${sample_name}.to_map.single | wc -l > ${sample_name}.04a.to_map.single.count
sed -i "s/^/$sample_name\ttoMap\tsingle\tall\t/" ${sample_name}.04a.to_map.single.count


if [ $mode == "process" ]
then
    mv ${sample_name}.to_map.single ../

    if [ ! -d ../logs ];then mkdir ../logs;fi
    mv ${sample_name}*log ../logs
    
    cd ../
    mv ${sample_name}.bash.log logs
    if [ $del == 1 ]
    then
	rm -r ${sample_name}.MAP
    fi
    exit 1;
fi


# do some cleanup - large files
rm ${sample_name}_trimmed.fq


#####################################################
############### ALIGN TO GENOME #####################
#####################################################

if [ $mode == "basic" ]
   then
       
       #########################
       #####SINGLE END READS####
       #########################

       echo "$br"
       echo "mapping single reads to rRNA with STAR"
       echo "br2"
	   # map to rRNA with STAR
       STAR --runThreadN $threads --readFilesIn ${sample_name}.to_map.single --genomeDir $rRNA_STAR --outSAMtype BAM Unsorted --outFileNamePrefix ${sample_name}.single.rRNA. --outReadsUnmapped Fastx --outSAMattributes All  --outFilterMismatchNoverLmax 0.05 --outFilterMatchNmin 16 --outFilterScoreMinOverLread 0.9  --outFilterMatchNminOverLread 0.9 --alignIntronMax 1


       mv ${sample_name}.single.rRNA.Log.final.out ${sample_name}.single.rRNA.STAR.log
       
       samtools view -@ $threads -F 260 -b ${sample_name}.single.rRNA.Aligned.out.bam | samtools sort -o ${sample_name}.single.rRNA.bam -

       # compress the fastqs
       gzip ${sample_name}.to_map.single

       echo "$br"
       echo "mapping single reads to snoRNA with bowtie2"
       echo "$br2"
       
       /home/emily/data/EJ25/code/map_ncRNA.sh single ${ncRNA_bowtie} snoRNA $sample_name $threads ${sample_name}.single.rRNA.Unmapped.out.mate1

       echo "$br"
       echo "mapping single reads to snRNA with bowtie2"
       echo "$br2"
       
       /home/emily/data/EJ25/code/map_ncRNA.sh single ${ncRNA_bowtie} snRNA $sample_name $threads ${sample_name}.single.snoRNA.remap.single

       echo "$br"
       echo "mapping single reads to piRNA with bowtie2"
       echo "$br2"
       
       /home/emily/data/EJ25/code/map_ncRNA.sh single ${ncRNA_bowtie} piRNA $sample_name $threads ${sample_name}.single.snRNA.remap.single

       gzip ${sample_name}.single.rRNA.Unmapped.out.mate1
       gzip ${sample_name}.single.snoRNA.remap.single
       gzip ${sample_name}.single.snRNA.remap.single

       echo "$br"
       echo "mapping single reads to tRNA with GSNAP"
       echo "$br2"
       
       
       gsnap -D $tRNA_GSNAP_D -d $tRNA_GSNAP_d -V $tRNA_SNP_V -v $tRNA_SNP_v -t $threads --split-output ${sample_name}.gsnap.single --ignore-trim-in-filtering 1 --format sam --genome-unk-mismatch 0 --failed-input ${sample_name}.gsnap.single.unmapped.single --md-lowercase-snp --max-mismatches 0.075 ${sample_name}.single.piRNA.remap.single

       samtools view -b -F 260 ${sample_name}.gsnap.single.unpaired_uniq | samtools sort -o ${sample_name}.gsnap.single.unpaired_uniq.bam
       samtools view -b -F 260 ${sample_name}.gsnap.single.unpaired_mult | samtools sort -o ${sample_name}.gsnap.single.unpaired_mult.bam
       samtools merge ${sample_name}.single.tRNA.bam ${sample_name}.gsnap.single.unpaired_uniq.bam ${sample_name}.gsnap.single.unpaired_mult.bam
	
       # map any remaining reads to the genome
       echo "$br"
       echo "mapping single reads to genome with bowtie2"
       echo "$br2"
       

       bowtie2 -p $threads -x $genome_bowtie -U ${sample_name}.gsnap.single.unmapped.single --very-sensitive 1> ${sample_name}.s.bowtie.Aligned.out.sam 2> ${sample_name}.single.genome.bowtie.log
       samtools view -@ $threads -b ${sample_name}.s.bowtie.Aligned.out.sam -F 256 > ${sample_name}.s.bowtie.bam
       
       samtools view -@ $threads -F 4 -b ${sample_name}.s.bowtie.bam > ${sample_name}.single.genome.bam

       # add read group tags so that origin of alignment is known
       samtools addreplacerg -r ID:rRNA_STAR -o ${sample_name}.single.rRNA.tag.bam ${sample_name}.single.rRNA.bam
       samtools addreplacerg -r ID:snoRNA_bowtie -o ${sample_name}.single.snoRNA.tag.bam ${sample_name}.single.snoRNA.bam
       samtools addreplacerg -r ID:snRNA_bowtie -o ${sample_name}.single.snRNA.tag.bam ${sample_name}.single.snRNA.bam
       samtools addreplacerg -r ID:piRNA_bowtie -o ${sample_name}.single.piRNA.tag.bam ${sample_name}.single.piRNA.bam
       samtools addreplacerg -r ID:tRNA_GSNAP -o ${sample_name}.single.tRNA.tag.bam ${sample_name}.single.tRNA.bam
       samtools addreplacerg -r ID:genome_bowtie -o ${sample_name}.single.bowtie.tag.bam ${sample_name}.single.genome.bam

       rm ${sample_name}.single.rRNA.bam ${sample_name}.single.snoRNA.bam ${sample_name}.single.snRNA.bam ${sample_name}.single.piRNA.bam ${sample_name}.single.tRNA.bam ${sample_name}.single.genome.bam
       
       # merge all mapped bams (rRNA, snoRNA, snRNA, tRNA genome bowtie) together
       samtools merge -@ $threads ${sample_name}.single.unsorted.bam ${sample_name}.single.rRNA.tag.bam ${sample_name}.single.snoRNA.tag.bam ${sample_name}.single.snRNA.tag.bam ${sample_name}.single.piRNA.tag.bam ${sample_name}.single.tRNA.tag.bam ${sample_name}.single.bowtie.tag.bam

       samtools sort -@ $threads -o ${sample_name}.single.bam ${sample_name}.single.unsorted.bam
       
       ##########################
       #### FINAL PROCESSING ####
       ##########################

       # if the sample has UMI, deduplicate them
       if [ $UMI != "none" ]
       then	   
	   if [ $del == 0 ]
	   then
	       cp ${sample_name}.single.bam ${sample_name}.single.before_dedup.bam
	       samtools index ${sample_name}.single.before_dedup.bam
	   fi
	   	   
	   samtools index ${sample_name}.single.bam
	   samtools idxstats ${sample_name}.single.bam > ${sample_name}.single.before_dedup.idxstats

	   echo "$br"
	   echo "deduplicating single alignments"
	   echo "$br2"
	   
	   umi_tools dedup --method unique -L ${sample_name}.single.dedup.log -I ${sample_name}.single.bam -S ${sample_name}.single.dedup.bam    
	   rm ${sample_name}.single.bam*
	   mv ${sample_name}.single.dedup.bam ${sample_name}.single.bam
	   
       fi
      
       # if bam needs to be flipped
       if [ $flip == 1 ]
       then
	   echo "$br"
	   echo "inverting bam strands"
	   echo "$br2"
	  	   
	   samtools view -h ${sample_name}.single.bam | awk 'BEGIN{FS=OFS="\t"} $1~/^@/ {print;next} and($2,16) {$2=$2-16; print;next} !and($2,16) {$2=$2+16;print}' | samtools view -@ $threads -bS - 1> ${sample_name}.s.inverted.bam
	   rm ${sample_name}.single.bam*
	   mv ${sample_name}.s.inverted.bam ${sample_name}.single.bam
	   
       fi
       
       mv ${sample_name}.single.bam ${sample_name}.single.unsorted.bam
       samtools sort -@ $threads -o ${sample_name}.unfiltered.single.bam ${sample_name}.single.unsorted.bam

       samtools view -h ${sample_name}.unfiltered.single.bam | awk '/^@/ || length($10) >= 15' | samtools view -b > ${sample_name}.single.bam
       samtools index ${sample_name}.single.bam
	
       ###############################
       ######## do some stats ########
       ###############################

       if [ $UMI != "none" ]
       then
     
	   grep -v "BK000964.3" ${sample_name}.single.before_dedup.idxstats | awk '{total +=$3} END {print total}' > ${sample_name}.05a.single.main.before_dedup.count
	   sed -i "s/^/$sample_name\tbeforeDedup\tsingle\tmain\t/" ${sample_name}.05a.single.main.before_dedup.count

	   grep "BK000964.3" ${sample_name}.single.before_dedup.idxstats | cut -f 3 > ${sample_name}.05b.single.rRNA.before_dedup.count
	   sed -i "s/^/$sample_name\tbeforeDedup\tsingle\trRNA\t/" ${sample_name}.05b.single.rRNA.before_dedup.count
       fi
       
       samtools idxstats ${sample_name}.single.bam > ${sample_name}.single.idxstats
 
       grep -v "BK000964.3" ${sample_name}.single.idxstats | awk '{total +=$3} END {print total}' > ${sample_name}.06a.single.main.count
       sed -i "s/^/$sample_name\tmap\tsingle\tmain\t/" ${sample_name}.06a.single.main.count

       grep "BK000964.3" ${sample_name}.single.idxstats | cut -f 3 > ${sample_name}.06b.single.rRNA.count
       sed -i "s/^/$sample_name\tmap\tsingle\trRNA\t/" ${sample_name}.06b.single.rRNA.count

       cat ${sample_name}*count > ${sample_name}.stats
       
       if [ ! -d ../bam ];then mkdir ../bam;fi
       if [ ! -d ../logs ];then mkdir ../logs;fi
       
       mv ${sample_name}.single.bam ../bam
       mv ${sample_name}*log ../logs
       mv ${sample_name}.stats ../logs
       
fi

       
date
cd ../
mv ${sample_name}.bash.log logs

if [ $del == 1 ]
then
    rm -r ${sample_name}.MAP
fi

