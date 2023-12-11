args <- commandArgs(trailingOnly = TRUE)

current_dir=args[1]

need_packages=c()

if(suppressMessages(suppressWarnings(!require(GenomicRanges)))==T){
  need_packages=c(need_packages,"GenomicRanges")
}

if(suppressMessages(suppressWarnings(!require(GenomicAlignments)))==T){
  need_packages=c(need_packages,"GenomicAlignments")
}

if(suppressMessages(suppressWarnings(!require(rtracklayer)))==T){
  need_packages=c(need_packages,"rtracklayer")
}

if(suppressMessages(suppressWarnings(!require(Rsamtools)))==T){
  need_packages=c(need_packages,"Rsamtools")
}
if(suppressMessages(suppressWarnings(!require(parallel)))==T){
  need_packages=c(need_packages,"parallel")
}
if(suppressMessages(suppressWarnings(!require(ggplot2)))==T){
  need_packages=c(need_packages,"ggplot2")
}
if(suppressMessages(suppressWarnings(!require(Biostrings)))==T){
  need_packages=c(need_packages,"Biostrings")
}



if(length(need_packages)>0){
  #print("The following R packages need to be installed:")
  #print(paste(need_packages,collapse=", "))
  
  write.table(paste(need_packages,collapse=", "),paste0(current_dir,"/LIDAR.R.need_packages"),col.names=F,row.names=F,quote=F,sep="\t")
}





