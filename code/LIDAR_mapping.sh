#!/bin/bash

# process fastq files with/out UMI and map to the genome

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
R2_3prime_adapter="none"
readS=""

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
	   echo 'Usage: ./LIDAR.sh [-h help] [-keep keep files (no)] [-flip flip orientation (no)] [-p threads (1)] [-c /path/to/cope if COPE not in path] [-lib lib_type <LIDAR, NEB, none> (LIDAR)] [-name sample_name (sample)] [--mode <basic, process, analyze_only> (basic)] [ -singles single ends reads (analyze_only)] [-species <mmus, hsap> (mmus)] -1 R1_reads -2 R2_reads'
	   echo
	   echo
	   echo 'OPTIONAL PARAMETERS:'
	   echo
	   echo '[-keep keep files (no)]                                        set -k to keep temp files'
	   echo '[-flip flip orientation (no)]                                  set -f if bam files need to be flipped so that reads are on the correct strand (Illumina pA RNA-seq)'
	   echo '[-p threads (1)]                                               number of threads'
	   echo '[-c /path/to/cope ()] - path to cope executable]               path to COPE, if not already in PATH'              
	   echo '[-lib (LIDAR)]                                                 lib type (LIDAR - LIDAR UMI, NEB - no UMI and with R2 adapter trimming,  none - no UMI)'
	   echo '[-name sample name (sample)]                                   base name for sample'
	   echo '[-mode <basic, tRNA, process> (basic)]                         choose mode: basic = main analysis pipeline, process = only process fastqs, analyze_only = you are providing processed fastqs with all adapter/UMI steps completed'
	   echo '[-species <hsap, mmus> (mmus)]                                 species (hsap = human hg38, mmus= mouse mm39), default is mouse'
	   echo '[-singles single reads (none)]                                 processed single-end reads - can only be used in analyze_only mode, otherwise provide as single in table'
	   echo
	   echo 'REQUIRED PARAMETERS'
	   echo
	   echo '-1                                                             R1 reads'
	   echo '-2                                                             R2 reads'
	   echo ''
 
	   echo 'SOFTWARE REQUIREMENTS:'
	   echo 'suggested conda environment: conda create -n LIDAR -c bioconda cutadapt trim-galore umi_tools bowtie2 samtools=1.16 gmap=2019.02.26 STAR=2.7.10a bbmap bedtools'
	   echo 'cope must be downloaded and placed in PATH, or path provided with -c <path/to/cope/'
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
	-2)
	    shift
	    reads2=$1
	    shift
	    ;;
	-singles)
	    shift
	    readS=$1
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
    R2_3prime_adapter="GATCGTCGG"
    UMI="none"
fi

if [ $lib == "none" ]
then
    UMI="none"
fi


if [ $lib == "LIDAR" ]
then
    UMI="LIDAR"
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
if test -f "$reads2" ; then
    : ; else
    echo "error! reads2 does not exist!"
    echo "the file name was given as $reads2 - check that this file exists"
    exit 1 ; fi


shift $(($OPTIND-1))


genome_bowtie=$index_dir/$species/genome_bowtie/genome_rRNA_tRNA
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
echo reads2 file: $reads2
echo same name: $sample_name
echo library type: $lib
echo UMI type from library: $UMI
echo species: $species
echo delete files after running: $del
echo threads for this sample: $threads
echo mode: $mode
echo cope path: $cope_path
echo invert bam file: $flip


if [ $mode != "analyze_only" ]
then
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
    trim_galore --paired -j $threads --length 5 --retain_unpaired --basename ${sample_name} -r1 20 -r2 100 -o $sample_name.MAP $reads1 $reads2
    
    cd $sample_name.MAP
    
    
    if [ $R2_3prime_adapter != "none" ]
    then
	echo "$br"
	echo "trimming R2 3' adapter"
	echo "$br2"
	
	mv ${sample_name}_val_1.fq.gz ${sample_name}_pretrim.R1.gz
	mv ${sample_name}_val_2.fq.gz ${sample_name}_pretrim.R2.gz
	
	mv ${sample_name}_R1_unpaired_1.fq.gz ${sample_name}_pretrim.single.gz
	# extra trimming step for remaining 3prime R2 adapter
	cutadapt -j $threads -A $R2_3prime_adapter --minimum-length 5 --pair-filter=any -o ${sample_name}_val_1.fq.gz -p ${sample_name}_val_2.fq.gz ${sample_name}_pretrim.R1.gz ${sample_name}_pretrim.R2.gz > ${sample_name}.cutadapt.3prime.log
	cutadapt -j $threads -a $R2_3prime_adapter --minimum-length 5 -o ${sample_name}_R1_unpaired_1.fq.gz ${sample_name}_pretrim.single.gz > ${sample_name}.cutadapt.unpaired.3prime.log
    fi
    
    echo "$br"
    echo 'collapsing reads with cope'
    echo "$br2"
    
    
    # connect overlapping reads with COPE
    cope -a ${sample_name}_val_1.fq.gz -b ${sample_name}_val_2.fq.gz -o ${sample_name}.COPE -s 33 -m 0 -l 8 > ${sample_name}.COPE.log ${sample_name}.COPE.error
    
    awk '{split($0,a,"_");print a[1]}' ${sample_name}.COPE.connect.fq > ${sample_name}.COPE.connect.fixed0.fq
    
    
    # combine any unpaired reads from the trimming with the COPE connected reads
    unpigz ${sample_name}_R1_unpaired_1.fq.gz
    
    
    cat ${sample_name}.COPE.connect.fixed0.fq ${sample_name}_R1_unpaired_1.fq > ${sample_name}.COPE.connect.fixed.fq
    
    if [ $UMI != "none" ]
    then
	echo "$br"
	echo 'extracting UMIs and trimming remainder of UMI construct, if needed'
	echo "$br2"
	$code_dir/process_files_single.sh ${sample_name}.COPE.connect.fixed.fq ${sample_name} $UMI
	$code_dir/process_files_paired.sh ${sample_name}.COPE.unConnect_1.fq ${sample_name}.COPE.unConnect_2.fq ${sample_name} $UMI
	
	# stats on UMI files
	
	grep "@" ${sample_name}.UMI | wc -l > ${sample_name}.03a.UMI.single.count
	sed -i "s/^/$sample_name\tUMI\tsingle\tall\t/" ${sample_name}.03a.UMI.single.count
	
	grep "@" ${sample_name}.UMI.R1 | wc -l > ${sample_name}.03b.UMI.paired.count
	sed -i "s/^/$sample_name\tUMI\tpaired\tall\t/" ${sample_name}.03b.UMI.paired.count
	
	# clean up of big files
	rm ${sample_name}.UMI.R1 ${sample_name}.UMI.R2 ${sample_name}.UMI
	rm ${sample_name}.to_map0.R1 ${sample_name}.to_map0.R2 ${sample_name}.to_map0.single
	rm ${sample_name}.to_map1.R1 ${sample_name}.to_map1.R2
	
    fi
    
    # do some stats and then clean up - large files that we don't need after doing stats
    
    zcat ../$reads1 | grep "@" | wc -l > ${sample_name}.00.original.count
    sed -i "s/^/$sample_name\toriginal\tall\tall\t/" ${sample_name}.00.original.count
    
    trimmed_count="$(zcat ${sample_name}_val_1.fq.gz | grep "@" | wc -l | awk '{print $1}')"
    unpaired_count="$(grep "@"  ${sample_name}_R1_unpaired_1.fq | wc -l | awk '{print $1}')"
    
    echo $((trimmed_count + unpaired_count)) > ${sample_name}.01.trimmed.count
    sed -i "s/^/${sample_name}\ttrimmed\tall\tall\t/" ${sample_name}.01.trimmed.count
    
    grep "@" ${sample_name}.COPE.connect.fixed.fq | wc -l > ${sample_name}.02a.COPE.single.count
    sed -i "s/^/$sample_name\tCOPE\tsingle\tall\t/" ${sample_name}.02a.COPE.single.count
    
    grep "@" ${sample_name}.COPE.unConnect_1.fq | wc -l > ${sample_name}.02b.COPE.paired.count
    sed -i "s/^/$sample_name\tCOPE\tpaired\tall\t/" ${sample_name}.02b.COPE.paired.count
    
    if [ $UMI == "none" ]
    then
	
	mv ${sample_name}.COPE.connect.fixed.fq ${sample_name}.to_map.single
	mv ${sample_name}.COPE.unConnect_1.fq ${sample_name}.to_map.R1
	mv ${sample_name}.COPE.unConnect_2.fq ${sample_name}.to_map.R2
    fi
    
    # stats on to map files
    grep "@" ${sample_name}.to_map.single | wc -l > ${sample_name}.04a.to_map.single.count
    sed -i "s/^/$sample_name\ttoMap\tsingle\tall\t/" ${sample_name}.04a.to_map.single.count
    
    grep "@" ${sample_name}.to_map.R1 | wc -l > ${sample_name}.04b.to_map.paired.count
    sed -i "s/^/$sample_name\ttoMap\tpaired\tall\t/" ${sample_name}.04b.to_map.paired.count
    
    if [ $mode == "process" ]
    then
	mv ${sample_name}.to_map.single ../
	mv ${sample_name}.to_map.R1 ../
	mv ${sample_name}.to_map.R2 ../
	
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
fi

if [ $mode == "analyze_only" ]
then
    cp $reads1 $sample_name.MAP/${sample_name}.to_map.R1
    cp $reads2 $sample_name.MAP/${sample_name}.to_map.R2
    if [ $readS != "" ]
    then
	cp $readS $sample_name.MAP/${sample_name}.to_map.single
    fi
	    
    cd $sample_name.MAP
  
    if [[ $reads1 =~ \.gz$ ]];
    then
	mv ${sample_name}.to_map.R1 ${sample_name}.to_map.R1.gz
	unpigz ${sample_name}.to_map.R1.gz
    fi

    if [[ $reads2 =~ \.gz$ ]];
    then
	mv ${sample_name}.to_map.R2 ${sample_name}.to_map.R2.gz
	unpigz ${sample_name}.to_map.R2.gz
    fi

    if [ $readS != "" ]
    then
	if [[ $readS =~ \.gz$ ]];
	then
	    mv ${sample_name}.to_map.single ${sample_name}.to_map.single.gz
	    unpigz ${sample_name}.to_map.single.gz
	fi
    fi   
fi


# do some cleanup - large files
#rm ${sample_name}*val*
#rm ${sample_name}.COPE.*fq


#exit 1

#####################################################
############### ALIGN TO GENOME #####################
#####################################################

if [ $mode == "basic" ] || [ $mode == "analyze_only" ]
   then

       ######################
       #### PAIRED READS ####
       ######################

       echo "$br"
       echo "mapping paired reads to rRNA with STAR"
       echo "$br2"
       
       # map to rRNA with STAR
       STAR --runThreadN $threads --readFilesIn ${sample_name}.to_map.R1 ${sample_name}.to_map.R2 --genomeDir $rRNA_STAR --outFileNamePrefix ${sample_name}.paired.rRNA. --outSAMtype BAM Unsorted --outReadsUnmapped Fastx --outSAMattributes All --outFilterMatchNmin 16 --alignIntronMax 1 --outFilterScoreMinOverLread 0.9 --outFilterMatchNminOverLread 0.9 --outFilterMismatchNoverLmax 0.05

       mv ${sample_name}.paired.rRNA.Log.final.out ${sample_name}.paired.rRNA.STAR.log
       
       samtools view -@ $threads -F 260 -f 2 -b ${sample_name}.paired.rRNA.Aligned.out.bam | samtools sort -o ${sample_name}.paired.rRNA.bam - 

       # compress the fastqs
       gzip ${sample_name}.to_map.R1
       gzip ${sample_name}.to_map.R2
       gzip ${sample_name}.paired.rRNA.Unmapped.out.mate1
       gzip ${sample_name}.paired.rRNA.Unmapped.out.mate2

       echo "$br"
       echo "mapping paired reads to snoRNA with bowtie2"
       echo "$br2"
       
       # map to ncRNA with bowtie
       $code_dir/map_ncRNA.sh paired ${ncRNA_bowtie} snoRNA $sample_name $threads ${sample_name}.paired.rRNA.Unmapped.out.mate1.gz ${sample_name}.paired.rRNA.Unmapped.out.mate2.gz

       echo "$br"
       echo "mapping paired reads to snRNA with bowtie2"
       echo "$br2"
       $code_dir/map_ncRNA.sh paired ${ncRNA_bowtie} snRNA $sample_name $threads ${sample_name}.paired.snoRNA.remap.R1.gz ${sample_name}.paired.snoRNA.remap.R2.gz

       echo "$br"
       echo "mapping paired reads to piRNA with bowtie2"
       echo "$br2"
       $code_dir/map_ncRNA.sh paired ${ncRNA_bowtie} piRNA $sample_name $threads ${sample_name}.paired.snRNA.remap.R1.gz ${sample_name}.paired.snRNA.remap.R2.gz

       echo "$br"
       echo "mapping paired reads to tRNA with GSNAP"
       echo "$br2"
       
       # map remaining reads to tRNA using gsnap
       gsnap --gunzip -D $tRNA_GSNAP_D -d $tRNA_GSNAP_d -V $tRNA_SNP_V -v $tRNA_SNP_v -t $threads --split-output ${sample_name}.gsnap.paired --ignore-trim-in-filtering 1 --format sam --genome-unk-mismatch 0 --failed-input ${sample_name}.gsnap.paired.unmapped --md-lowercase-snp --max-mismatches 0.075 ${sample_name}.paired.piRNA.remap.R1.gz ${sample_name}.paired.piRNA.remap.R2.gz

       samtools view -b -F 260 ${sample_name}.gsnap.paired.concordant_uniq | samtools sort -o ${sample_name}.gsnap.paired.concordant_uniq.bam -
       samtools view -b -F 260 ${sample_name}.gsnap.paired.concordant_mult | samtools sort -o ${sample_name}.gsnap.paired.concordant_mult.bam -
       samtools merge ${sample_name}.paired.tRNA.bam ${sample_name}.gsnap.paired.concordant_uniq.bam ${sample_name}.gsnap.paired.concordant_mult.bam

        
       # get any reads that are still unmapped and map to the genome  - get reads from halfmapped uniq and and halfmapped mult, cat with the unmapped 
       samtools view -b ${sample_name}.gsnap.paired.halfmapping_uniq | samtools sort -@ $threads -o ${sample_name}.gsnap.paired.halfmapping_uniq.bam -
       samtools view -b ${sample_name}.gsnap.paired.halfmapping_mult | samtools sort -@ $threads -o ${sample_name}.gsnap.paired.halfmapping_mult.bam -
       
       samtools merge -@ $threads ${sample_name}.gsnap.paired.halfmapping.unsorted.bam ${sample_name}.gsnap.paired.halfmapping_uniq.bam ${sample_name}.gsnap.paired.halfmapping_mult.bam 
       samtools sort -n -o ${sample_name}.gsnap.paired.halfmapping.bam ${sample_name}.gsnap.paired.halfmapping.unsorted.bam
       samtools bam2fq -1 ${sample_name}.paired.gsnap.halfmapping.R1 -2 ${sample_name}.paired.gsnap.halfmapping.R2 -s ${sample_name}.paired.gsnap.singleton.gz -n ${sample_name}.gsnap.paired.halfmapping.bam
       
       # finally cat these with the unmapped, and get rid of the /1 and /2 that gsnap added to the unmapped
       cat ${sample_name}.gsnap.paired.unmapped.1 ${sample_name}.paired.gsnap.halfmapping.R1 | sed 's/\/1//' > ${sample_name}.tRNA.unmapped.R1
       cat ${sample_name}.gsnap.paired.unmapped.2 ${sample_name}.paired.gsnap.halfmapping.R2 | sed 's/\/2//' > ${sample_name}.tRNA.unmapped.R2

       echo "$br"
       echo "mapping paired reads to genome with STAR"
       echo "$br2"
       
       STAR --peOverlapNbasesMin 5 --runThreadN $threads --alignIntronMin 20 --alignIntronMax $intron_max --readFilesIn ${sample_name}.tRNA.unmapped.R1 ${sample_name}.tRNA.unmapped.R2 --genomeDir $genome_STAR --outFileNamePrefix ${sample_name}.STAR. --outReadsUnmapped Fastx --outSAMtype BAM Unsorted --outFilterMismatchNoverLmax 0.05 --outFilterScoreMinOverLread 0.9 --outFilterMatchNminOverLread 0.9 --outFilterMatchNmin 16 --outSAMunmapped Within --outSAMattributes All

       mv ${sample_name}.STAR.Log.final.out ${sample_name}.paired.genome.STAR.log
	   
       # get the reads that map with STAR
       samtools view -b -F 256 -f 2 ${sample_name}.STAR.Aligned.out.bam > ${sample_name}.STAR.all_sizes.bam
       
       # get the reads that don't map - will remap with bowtie
       samtools view -b -f 4 ${sample_name}.STAR.Aligned.out.bam > ${sample_name}.STAR.unmapped.bam
       
       # split the paired-end reads into large (>= cutoff) and small (<cutoff)
       neg_length_cutoff=$((0-length_cutoff))
       
       samtools view -h ${sample_name}.STAR.all_sizes.bam | awk -v cutoff=$length_cutoff -v neg_cutoff=$neg_length_cutoff 'substr($0,1,1)=="@" || ($9< cutoff && $9>0) || ($9> neg_cutoff && $9<0)' | samtools view -b > ${sample_name}.STAR.small.bam
       
       samtools view -h ${sample_name}.STAR.all_sizes.bam | awk -v cutoff=$length_cutoff -v neg_cutoff=$neg_length_cutoff 'substr($0,1,1)=="@" || ($9 >= cutoff) || ($9 <= neg_cutoff)' | samtools view -b > ${sample_name}.STAR.large.bam
       
       
       # add the short insert size to the unmapped ones
       samtools merge -@ $threads ${sample_name}.to_remap.bam ${sample_name}.STAR.small.bam ${sample_name}.STAR.unmapped.bam
       
       # convert the ones with short insert size to fastq and remap using bowtie2
       samtools sort -n ${sample_name}.to_remap.bam -o ${sample_name}.to_remap.namesort.bam
       samtools bam2fq -1 ${sample_name}.to_remap.R1.gz -2 ${sample_name}.to_remap.R2.gz -s ${sample_name}.to_remap.singleton.gz -n ${sample_name}.to_remap.namesort.bam

       echo "$br"
       echo "mapping unmapped/short insert paired reads to genome with bowtie2"
       echo "$br2"
       
       # mapping with bowtie for reads with small insert size
       bowtie2 -p $threads -x $genome_bowtie -1 ${sample_name}.to_remap.R1.gz -2 ${sample_name}.to_remap.R2.gz --very-sensitive 1> ${sample_name}.bowtie.Aligned.out.sam 2> ${sample_name}.paired.genome.bowtie.log --un-conc-gz ${sample_name}.bowtie.paired.unmapped.gz
       samtools view -@ $threads -b ${sample_name}.bowtie.Aligned.out.sam -F 256 -f 2 > ${sample_name}.bowtie.bam

       # add read group tags so that origin of alignment is known
       samtools addreplacerg -r ID:rRNA_STAR -o ${sample_name}.paired.rRNA.tag.bam ${sample_name}.paired.rRNA.bam
       samtools addreplacerg -r ID:snoRNA_bowtie -o ${sample_name}.paired.snoRNA.tag.bam ${sample_name}.paired.snoRNA.bam
       samtools addreplacerg -r ID:snRNA_bowtie -o ${sample_name}.paired.snRNA.tag.bam ${sample_name}.paired.snRNA.bam
       samtools addreplacerg -r ID:piRNA_bowtie -o ${sample_name}.paired.piRNA.tag.bam ${sample_name}.paired.piRNA.bam
       samtools addreplacerg -r ID:tRNA_GSNAP -o ${sample_name}.paired.tRNA.tag.bam ${sample_name}.paired.tRNA.bam
       samtools addreplacerg -r ID:genome_STAR -o ${sample_name}.paired.STAR.tag.bam ${sample_name}.STAR.large.bam
       samtools addreplacerg -r ID:genome_bowtie -o ${sample_name}.paired.bowtie.tag.bam ${sample_name}.bowtie.bam

       # remove the untagged bams 
       rm ${sample_name}.paired.rRNA.bam ${sample_name}.paired.snoRNA.bam ${sample_name}.paired.snRNA.bam ${sample_name}.paired.piRNA.bam ${sample_name}.paired.tRNA.bam ${sample_name}.STAR.large.bam ${sample_name}.bowtie.bam
                    
       # merge all mapped bams (rRNA, snoRNA, snRNA, tRNA, genome STAR, genome bowtie) together
       samtools merge -@ $threads ${sample_name}.paired.unsorted.bam ${sample_name}.paired.rRNA.tag.bam ${sample_name}.paired.snoRNA.tag.bam ${sample_name}.paired.snRNA.tag.bam ${sample_name}.paired.piRNA.tag.bam ${sample_name}.paired.tRNA.tag.bam ${sample_name}.paired.STAR.tag.bam ${sample_name}.paired.bowtie.tag.bam
       
       samtools sort -@ $threads -o ${sample_name}.paired.bam ${sample_name}.paired.unsorted.bam
       
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
       
       $code_dir/map_ncRNA.sh single ${ncRNA_bowtie} snoRNA $sample_name $threads ${sample_name}.single.rRNA.Unmapped.out.mate1

       echo "$br"
       echo "mapping single reads to snRNA with bowtie2"
       echo "$br2"
       
       $code_dir/map_ncRNA.sh single ${ncRNA_bowtie} snRNA $sample_name $threads ${sample_name}.single.snoRNA.remap.single

       echo "$br"
       echo "mapping single reads to piRNA with bowtie2"
       echo "$br2"
       
       $code_dir/map_ncRNA.sh single ${ncRNA_bowtie} piRNA $sample_name $threads ${sample_name}.single.snRNA.remap.single

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
	   	   
	   samtools index ${sample_name}.paired.bam
	   samtools idxstats ${sample_name}.paired.bam > ${sample_name}.paired.before_dedup.idxstats
	   
	   
	   if [ $del == 0 ]
	   then
	       cp ${sample_name}.paired.bam ${sample_name}.paired.before_dedup.bam
	       cp ${sample_name}.single.bam ${sample_name}.single.before_dedup.bam
		   
		   samtools index ${sample_name}.single.before_dedup.bam
		   samtools index ${sample_name}.paired.before_dedup.bam
	   fi

	   echo "$br"
	   echo "deduplicating paired alignments"
	   echo "$br2"
	   
	   
	   umi_tools dedup --paired --method unique -L ${sample_name}.paired.dedup.log -I ${sample_name}.paired.bam -S ${sample_name}.paired.dedup.bam
	   rm ${sample_name}.paired.bam*
	   mv ${sample_name}.paired.dedup.bam ${sample_name}.paired.bam
	   
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
	   
	   samtools view -h ${sample_name}.paired.bam | awk 'BEGIN{FS=OFS="\t"} $1~/^@/ {print;next} and($2,128) {$2=$2-64; print;next} and($2,64) {$2=$2+64;print}' | samtools view -@ $threads -bS - 1> ${sample_name}.inverted.bam
	   rm ${sample_name}.bam*
	   mv ${sample_name}.inverted.bam ${sample_name}.paired.bam
	   
	   samtools view -h ${sample_name}.single.bam | awk 'BEGIN{FS=OFS="\t"} $1~/^@/ {print;next} and($2,16) {$2=$2-16; print;next} !and($2,16) {$2=$2+16;print}' | samtools view -@ $threads -bS - 1> ${sample_name}.s.inverted.bam
	   rm ${sample_name}.single.bam*
	   mv ${sample_name}.s.inverted.bam ${sample_name}.single.bam
	   
       fi
       
       mv ${sample_name}.paired.bam ${sample_name}.paired.unsorted.bam
       mv ${sample_name}.single.bam ${sample_name}.single.unsorted.bam
       
       samtools sort -@ $threads -o ${sample_name}.unfiltered.paired.bam ${sample_name}.paired.unsorted.bam 
       samtools sort -@ $threads -o ${sample_name}.unfiltered.single.bam ${sample_name}.single.unsorted.bam

       samtools view -h ${sample_name}.unfiltered.paired.bam | awk 'substr($0,1,1)=="@" || ($9>=15) || ($9<=-15)' | samtools view -b > ${sample_name}.paired.bam
       samtools view -h ${sample_name}.unfiltered.single.bam | awk '/^@/ || length($10) >= 15' | samtools view -b > ${sample_name}.single.bam
                   
       samtools index ${sample_name}.paired.bam
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
	   
	   grep -v "BK000964.3" ${sample_name}.paired.before_dedup.idxstats | awk '{total +=$3} END {print total}' > ${sample_name}.paired.main.before_dedup.count0
	   awk -v c=2 '{print $1/c}' ${sample_name}.paired.main.before_dedup.count0 > ${sample_name}.05c.paired.main.before_dedup.count
	   sed -i "s/^/$sample_name\tbeforeDedup\tpaired\tmain\t/" ${sample_name}.05c.paired.main.before_dedup.count

	   grep "BK000964.3" ${sample_name}.paired.before_dedup.idxstats | cut -f 3 > ${sample_name}.paired.rRNA.before_dedup.count0
	   awk -v c=2 '{print $1/c}' ${sample_name}.paired.rRNA.before_dedup.count0 > ${sample_name}.05d.paired.rRNA.before_dedup.count
	   sed -i "s/^/$sample_name\tbeforeDedup\tpaired\trRNA\t/" ${sample_name}.05d.paired.rRNA.before_dedup.count
       fi
       
       samtools idxstats ${sample_name}.single.bam > ${sample_name}.single.idxstats
       samtools idxstats ${sample_name}.paired.bam > ${sample_name}.paired.idxstats

       grep -v "BK000964.3" ${sample_name}.single.idxstats | awk '{total +=$3} END {print total}' > ${sample_name}.06a.single.main.count
       sed -i "s/^/$sample_name\tmap\tsingle\tmain\t/" ${sample_name}.06a.single.main.count

       grep "BK000964.3" ${sample_name}.single.idxstats | cut -f 3 > ${sample_name}.06b.single.rRNA.count
       sed -i "s/^/$sample_name\tmap\tsingle\trRNA\t/" ${sample_name}.06b.single.rRNA.count

       grep -v "BK000964.3" ${sample_name}.paired.idxstats | awk '{total +=$3} END {print total}' > ${sample_name}.paired.main.count0
       awk -v c=2 '{print $1/c}' ${sample_name}.paired.main.count0 > ${sample_name}.06c.paired.main.count
       sed -i "s/^/$sample_name\tmap\tpaired\tmain\t/" ${sample_name}.06c.paired.main.count

       grep "BK000964.3" ${sample_name}.paired.idxstats | cut -f 3 > ${sample_name}.paired.rRNA.count0
       awk -v c=2 '{print $1/c}' ${sample_name}.paired.rRNA.count0> ${sample_name}.06d.paired.rRNA.count
       sed -i "s/^/$sample_name\tmap\tpaired\trRNA\t/" ${sample_name}.06d.paired.rRNA.count

       cat ${sample_name}*count > ${sample_name}.stats
       
       if [ ! -d ../bam ];then mkdir ../bam;fi
       if [ ! -d ../logs ];then mkdir ../logs;fi
       
       mv ${sample_name}.paired.bam ../bam
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

