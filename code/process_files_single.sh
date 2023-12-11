#!/bin/bash

reads=$1
sample_name=$2
UMI_type=$3


if [ $UMI_type == "LIDAR" ]
then    
    sed -i '2~4s/^/XXX/' $reads
    sed -i '4~4s/^/AAA/' $reads
    
    umi_tools extract --extract-method=regex --bc-pattern=".*(?<umi_1>.{7})CG(?P<umi_2>.{2})AG(?P<umi_3>.{2})GGG" -L $sample_name.single_UMI_extract.log -I $reads --stdout ${sample_name}.UMI
    cutadapt --minimum-length 5 -g CGAGGGG -o ${sample_name}.to_map0.single ${sample_name}.UMI > ${sample_name}.single_cutadapt_UMI_pattern_trim.log
fi
    

cutadapt --revcomp -e 3 -b CGTCAGATGTGTATAAGAGACAG --discard-trimmed -o ${sample_name}.to_map.single ${sample_name}.to_map0.single > ${sample_name}.cutadapt.single_adapter.log


