# LIDAR
Mapping and basic analysis to accompany LIDAR-seq

A ligation-independent sequencing method reveals tRNA-derived RNAs with blocked 3' termini

Scacchetti, A., Shields, E.J., Trigg, N.A., Wilusz, J.E., Conine, C.C., and Bonasio, R (2023). A ligation-independent sequencing method reveals tRNA-derived RNAs with blocked 3' termini. Preprint at bioRxiv, 10.1101/2023.06.06.543899.


-------------------------------
Installation instructions:

No complilation is required, simply download and run with instructions below. The following dependencies are required:
-COPE (https://sourceforge.net/projects/coperead/) - place in PATH or provide its location with -c <path/to/cope/>
-trim_galore
-samtools ver > 1.16
-UMI_tools
-cutadapt
-bbmap
-STAR ver 2.7.10a
-GMAP ver 2019.02.26
-bowtie2

Other than COPE, these can be installed with conda:
conda create -n LIDAR -c bioconda cutadapt trim-galore umi_tools bowtie2 samtools=1.16 gmap=2019.02.26 STAR=2.7.10a bbmap


R packages for downstream analysis (can be run without downstream analysis)
-GenomicRanges
-GenomicAlignments
-rtracklayer
-Rsamtools
-parallel
-ggplot2
-Biostrings

-------------------------------

Running instructions:

A sample table is required, with one row for each column and 4 (or 3 in single-end mode) columns:
sample_name    R1_file	    R2_file    library_type

Example running command, with mouse data and 10 threads provided:
./LIDAR.sh -samples sample_table -species mmus -p 10

To see all options, run ./LIDAR.sh -h

After run is complete, logs and # of reads mapped at each step (*.stats) can be found in the logs folder


-------------------------------

Sample data:

A small subset of data from the LIDAR publication is provided as a test set, along with a sample table. This data comes from mES serum, 50nt fraction cells using the LIDAR or NEB library preparation.







