# LIDAR
Mapping and basic analysis to accompany LIDAR-seq

A ligation-independent sequencing method reveals tRNA-derived RNAs with blocked 3' termini

Scacchetti, A., Shields, E.J., Trigg, N.A., Wilusz, J.E., Conine, C.C., and Bonasio, R (2023). A ligation-independent sequencing method reveals tRNA-derived RNAs with blocked 3' termini. Preprint at bioRxiv, 10.1101/2023.06.06.543899.


-------------------------------
Installation instructions:

No compilation is required, but some setup is required for the indexes (see below). The following dependencies are required:<br />
-COPE (https://sourceforge.net/projects/coperead/) - place in PATH or provide its location with -c <path/to/cope/><br />
-trim_galore<br />
-bedtools
-samtools ver > 1.16<br />
-UMI_tools<br />
-cutadapt<br />
-bbmap<br />
-STAR ver 2.7.10a<br />
-GMAP ver 2019.02.26<br />
-bowtie2<br />
<br />
Other than COPE, these can be installed with conda:<br />
conda create -n LIDAR -c bioconda cutadapt trim-galore umi_tools bowtie2 samtools=1.16 gmap=2019.02.26 STAR=2.7.10a bbmap bedtools<br />
<br />
<br />
R packages for downstream analysis (can be run without downstream analysis)<br />
-GenomicRanges<br />
-GenomicAlignments<br />
-rtracklayer<br />
-Rsamtools<br />
-parallel<br />
-ggplot2<br />
-Biostrings<br />
<br />
-------------------------------<br />
<br />
Setup instructions:
Genome indexes for the genome must be created before you run for the first time, which will take some time. First, download the genome fasta for mm39 mouse or hg38 human from the appropriate bigZips directory and place in the index/mmus/ or index/hsap folder.<br />
mmus: https://hgdownload.soe.ucsc.edu/goldenPath/mm39/bigZips/mm39.fa.gz   <br />
hsap: https://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/hg38.fa.gz  <br />

Then, with STAR 2.7.10a, bowtie2, and bedtools available, run the index setup code for the species of interest using code/index_setup.sh  <mmus/hsap> <# threads>  < br />

<br />
-------------------------------<br />
<br />
Running instructions:<br />

A sample table is required, with one row for each column and 4 (or 3 in single-end mode) columns:<br />
sample_name    R1_file	    R2_file    library_type<br />
<br />
Example running command, with mouse data and 10 threads provided:<br />
./LIDAR.sh -samples sample_table -species mmus -p 10<br />
<br />
To see all options, run ./LIDAR.sh -h<br />
<br />
After run is complete, logs and # of reads mapped at each step (*.stats) can be found in the logs folder<br />
<br />
<br />
-------------------------------<br />
<br />
Sample data:<br />
<br />
A small subset of data from the LIDAR publication is provided as a test set, along with a sample table. This data comes from mES serum, 50nt fraction cells using the LIDAR or NEB library preparation.<br />
<br />






