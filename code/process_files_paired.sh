#!/bin/bash


reads1=$1
reads2=$2
sample_name=$3
UMI_type=$4


if [ $UMI_type == "LIDAR" ]
then
    sed -i '2~4s/^/XXX/' $reads1
    sed -i '4~4s/^/AAA/' $reads1

    umi_tools extract --extract-method=regex --bc-pattern=".*(?<umi_1>.{7})CG(?P<umi_2>.{2})AG(?P<umi_3>.{2})GGG" -L $sample_name.paired_UMI_extract.log -I $reads1 --stdout ${sample_name}.UMI.R1 --read2-in=$reads2 --read2-out=${sample_name}.UMI.R2
    
    cutadapt --pair-filter=any --minimum-length 5 -g CGAGGGG -o ${sample_name}.to_map0.R1 -p ${sample_name}.to_map0.R2 ${sample_name}.UMI.R1 ${sample_name}.UMI.R2 > ${sample_name}.paired_cutadapt_UMI_pattern_trim.log
   
fi


cutadapt -e 3 -b CGTCAGATGTGTATAAGAGACAG --pair-filter=any --minimum-length 10 --discard-trimmed -o ${sample_name}.to_map1.R1 -p ${sample_name}.to_map1.R2 ${sample_name}.to_map0.R1 ${sample_name}.to_map0.R2 > ${sample_name}.cutadapt.paired_adapter1.log
cutadapt -e 3 -b CTGTCTCTTATACACATCTGACG --pair-filter=any --minimum-length 10 --discard-trimmed -o ${sample_name}.to_map.R1 -p ${sample_name}.to_map.R2 ${sample_name}.to_map1.R1 ${sample_name}.to_map1.R2 > ${sample_name}.cutadapt.paired_adapter2.log
   
