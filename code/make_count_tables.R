#libraries and directories

suppressPackageStartupMessages(library(GenomicRanges))
suppressPackageStartupMessages(library(GenomicAlignments))
suppressPackageStartupMessages(library(rtracklayer))
suppressPackageStartupMessages(library(Rsamtools))
suppressPackageStartupMessages(library(parallel))
suppressPackageStartupMessages(library(ggplot2))



options(stringsAsFactors = FALSE)
theme_set(theme_classic())

get_split<-function(x,split_char,index){return(unlist(strsplit(x,split_char,fixed=T))[index])}
get_sem<-function(x){return(sd(x)/sqrt(length(x)))}

make_data_frame<-function(data,x_levels=NULL,y_levels=NULL,x_label="x",y_label="y",value_label="value",y_flip=F,x_flip=F){
  if(is.null(x_levels)){x_levels=colnames(data)}
  if(is.null(y_levels)){y_levels=rownames(data)}
  
  if(y_flip){y_levels=rev(y_levels)}
  if(x_flip){x_levels=rev(x_levels)}  
  df<-data.frame(x=factor(c(sapply(colnames(data),rep,nrow(data))),
                          levels=x_levels),
                 y=factor(rep(rownames(data),ncol(data)),levels=y_levels),
                 value=c(data))
  if(x_label!="x"){colnames(df)[1]=x_label}
  if(y_label!="y"){colnames(df)[2]=y_label}
  if(value_label!="value"){colnames(df)[3]=value_label}
  return(df)
  
}

args <- commandArgs(trailingOnly = TRUE)

bam_dir=args[1]
annotation_dir=args[2]
index_dir=args[3]
threads=as.numeric(args[4])
species=args[5]
single=as.numeric(args[6])


if(species=="mmus"){
  annotation.genes=readRDS(paste0(annotation_dir,"/GRCm39.genes.rds"))
  annotation.regions=readRDS(paste0(annotation_dir,"/GRCm39.regions.rds"))
  rRNA.name="BK000964.3"
}
if(species=="hsap"){
  annotation.genes=readRDS(paste0(annotation_dir,"/GRCh38.genes.rds"))
  annotation.regions=readRDS(paste0(annotation_dir,"/GRCh38.regions.rds"))
  rRNA.name="U13369.1"
}




# split by gene
genes.split=split(annotation.genes,annotation.genes$gene)


# import bams and set bam parameters
single.bams<-list.files(bam_dir,pattern=".single.bam$",full.names=T)
names(single.bams)<-gsub(".single.bam","",sapply(single.bams,get_split,"//",2))

if(single==0){
  paired.bams<-list.files(bam_dir,pattern=".paired.bam$",full.names=T)
  names(paired.bams)<-gsub(".paired.bam","",sapply(paired.bams,get_split,"//",2))
}

#use only primary alignments
param=ScanBamParam(flag=scanBamFlag(isSecondaryAlignment = FALSE))
stranded=T

#function to count reads

count_reads<-function(reads,annotation,stranded=stranded,singleEnd,param=param){
  return(assay(summarizeOverlaps(annotation,reads,ignore.strand=!stranded,singleEnd=singleEnd,param=param,mode="IntersectionNotEmpty")))
}


# # need to split the bams into chunks to run the counting - if you don't the parlapply cl will crash with mysterious "error writing to connection." #From my tests it can only handle 20-30 files at once, sometimes less. I think it's a memory thing, I have been trying to fix for 6 months. #Changing the number of cl doesn't help. This was the simplest way I could get around it. 

get_chunk<- function(x,n) split(x, cut(seq_along(x), n, labels = FALSE))
if(length(single.bams)==1){
  chunks=1
}else{
  chunks=get_chunk(1:length(single.bams),max(length(single.bams),length(single.bams)/threads))
}

cl<-makeCluster(threads)
clusterExport(cl, list("summarizeOverlaps","assay","chunks"))

single.counts0<-lapply(chunks,function(x){
  return(do.call("cbind",parLapply(cl=cl,single.bams[x],count_reads,annotation=genes.split,stranded=stranded,singleEnd=T,param=param)))
})
single.counts<-do.call(cbind,single.counts0)

if(single==0){
  paired.counts0<-lapply(chunks,function(x){
    return(do.call("cbind",parLapply(cl=cl,paired.bams[x],count_reads,annotation=genes.split,stranded=stranded,singleEnd=F,param=param)))
  })
  paired.counts<-do.call(cbind,paired.counts0)
}

stopCluster(cl)

colnames(single.counts)<-names(single.bams)

if(single==0){
  colnames(paired.counts)<-names(paired.bams)
  all.counts1=paired.counts+single.counts[,colnames(paired.counts)]
}else{
  all.counts1=single.counts
}

### if there is only one sample, just put two columns of the counts to make it easier, then will only return one of them
if(length(single.bams)==1){
  all.counts1=cbind(all.counts1,all.counts1)
}


## counts on the snRNA, snoRNA, piRNA, tRNA (separate scaffolds)

single.idxstats0<-list.files(paste0(bam_dir,"idxstats/"),pattern="single.idxstats",full.names=T)
names(single.idxstats0)<-gsub(".single.idxstats","",sapply(single.idxstats0,get_split,"//",2))
single.idxstats<-lapply(single.idxstats0,read.delim,header=F,row.names=1)


# don't include E. coli chr
chrs<-setdiff(rownames(single.idxstats[[1]]),"NC_000913.3")

single.idxstats.counts<-sapply(single.idxstats,function(x){return(x[chrs,]$V3)})
rownames(single.idxstats.counts)<-chrs

# same as above: if there is only one sample, just copy it to make two to make everything downstream easier
if(length(single.idxstats)==1){
  single.idxstats.counts=cbind(single.idxstats.counts,single.idxstats.counts)
}





if(single==0){
  paired.idxstats0<-list.files(paste0(bam_dir,"idxstats/"),pattern="paired.idxstats",full.names=T)
  names(paired.idxstats0)<-gsub(".paired.idxstats","",sapply(paired.idxstats0,get_split,"//",2))
  paired.idxstats<-lapply(paired.idxstats0,read.delim,header=F,row.names=1)
  paired.idxstats.counts<-sapply(paired.idxstats,function(x){return(x[chrs,]$V3/2)})
  rownames(paired.idxstats.counts)<-chrs
  if(length(paired.idxstats)==1){
    paired.idxstats.counts=cbind(paired.idxstats.counts,paired.idxstats.counts)
  }
  
}



piRNA.names<-read.delim(paste0(index_dir,species,"/ncRNA_bowtie/piRNA.fai"),header=F)$V1
snoRNA.names<-read.delim(paste0(index_dir,species,"/ncRNA_bowtie/snoRNA.fai"),header=F)$V1
snRNA.names<-read.delim(paste0(index_dir,species,"/ncRNA_bowtie/snRNA.fai"),header=F)$V1

tRNA.contig=list.files(paste0(index_dir,species,"/tRNA_GSNAP/"))[grepl(".contig$",list.files(paste0(index_dir,species,"/tRNA_GSNAP/")))]
tRNA.names<-read.delim(paste0(index_dir,species,"/tRNA_GSNAP/",tRNA.contig),header=F)$V1


if(single==0){
  piRNA.counts<-single.idxstats.counts[piRNA.names,]+paired.idxstats.counts[piRNA.names,]
  snoRNA.counts0<-single.idxstats.counts[snoRNA.names,]+paired.idxstats.counts[snoRNA.names,]
  snRNA.counts0<-single.idxstats.counts[snRNA.names,]+paired.idxstats.counts[snRNA.names,]
  tRNA.counts<-single.idxstats.counts[tRNA.names,]+paired.idxstats.counts[tRNA.names,]
  rRNA.counts<-single.idxstats.counts[rRNA.name,]+paired.idxstats.counts[rRNA.name,]
}else{
  piRNA.counts<-single.idxstats.counts[piRNA.names,]
  snoRNA.counts0<-single.idxstats.counts[snoRNA.names,]
  snRNA.counts0<-single.idxstats.counts[snRNA.names,]
  tRNA.counts<-single.idxstats.counts[tRNA.names,]
  rRNA.counts<-single.idxstats.counts[rRNA.name,]
}

# combine genes that have a dash
snRNA.count.names=sapply(rownames(snRNA.counts0),function(x){
  x.split=get_split(x,"-")
  return(paste(x.split[1:(length(x.split)-1)],collapse="-"))
})
snRNA.counts=t(sapply(unique(snRNA.count.names),function(x){
  select=snRNA.counts0[snRNA.count.names==x,]
  if(is.null(nrow(select))){
    return(select)
  }
  return(colSums(select))
}))


# combine genes that have a dash
snoRNA.count.names=sapply(rownames(snoRNA.counts0),function(x){
  x.split=get_split(x,"-")
  return(paste(x.split[1:(length(x.split)-1)],collapse="-"))
})
snoRNA.counts=t(sapply(unique(snoRNA.count.names),function(x){
  select=snoRNA.counts0[snoRNA.count.names==x,]
  if(is.null(nrow(select))){
    return(select)
  }
  return(colSums(select))
}))


all.counts2<-rbind(piRNA.counts,snoRNA.counts,snRNA.counts,tRNA.counts)

## combine annotation counts with ncRNA

all.counts=rbind(all.counts1,all.counts2[,colnames(all.counts1)])

# there are some genes that are listed as snoRNA with GENCODE but pseudogenes in refseq, so they would be both in the annotation and listed as snoRNA. Combine their counts.

all.counts.dup_genes=unique(rownames(all.counts)[duplicated(rownames(all.counts))])

all.counts.no_dup=all.counts[setdiff(rownames(all.counts),all.counts.dup_genes),]
all.counts.dup=t(sapply(all.counts.dup_genes,function(x){
  select.rows=all.counts[rownames(all.counts)==x,]
  return(colSums(select.rows))
}))

all.counts.final=rbind(all.counts.no_dup,all.counts.dup)

# sort by name
all.counts.final=all.counts.final[order(rownames(all.counts.final)),]

# if only one sample, take just the first column since second is dup of first
if(length(single.idxstats)==1){
  all.counts.final=cbind(all.counts.final[,1])
  colnames(all.counts.final)=names(single.idxstats)
}

## save files
saveRDS(all.counts.final,paste0(bam_dir,"../counts/",format(Sys.Date(),"%y%m%d"),".gene.counts.rds"))
write.table(all.counts.final,paste0(bam_dir,"../counts/",format(Sys.Date(),"%y%m%d"),".gene.counts.txt"),col.names=NA,row.names=T,quote=F,sep="\t")


################# by biotypes

cl<-makeCluster(threads)
clusterExport(cl, list("summarizeOverlaps","assay"))

single.regions0<-lapply(chunks,function(x){
  return(do.call("cbind",parLapply(cl=cl,single.bams[x],count_reads,annotation=annotation.regions,stranded=stranded,singleEnd=T,param=param)))
})
single.regions<-do.call(cbind,single.regions0)

if(single==0){
  paired.regions0<-lapply(chunks,function(x){
    return(do.call("cbind",parLapply(cl=cl,paired.bams[x],count_reads,annotation=annotation.regions,stranded=stranded,singleEnd=F,param=param)))
  })
  paired.regions<-do.call(cbind,paired.regions0)
}
stopCluster(cl)


rownames(single.regions)=names(annotation.regions)
colnames(single.regions)=names(single.bams)

if(single==0){
  rownames(paired.regions)=names(annotation.regions)
  colnames(paired.regions)=names(paired.bams)
  all.regions1=single.regions+paired.regions[,colnames(single.regions)]
}else{
  all.regions1=single.regions
}


if(length(single.bams)==1){
  all.regions1=cbind(all.regions1,all.regions1)
}

all.regions2=rbind(piRNA=colSums(piRNA.counts),
                   snoRNA=colSums(snoRNA.counts),
                   snRNA=colSums(snRNA.counts),
                   tRNA=colSums(tRNA.counts[!grepl("mito",rownames(tRNA.counts)),]),
                   mt_tRNA=colSums(tRNA.counts[grepl("mito",rownames(tRNA.counts)),]),
                   rRNA=rRNA.counts)

all.regions3=all.regions1[!rownames(all.regions1) %in% c("snRNA-annotation","snoRNA-annotation",
                                                         "tRNA-annotation","Mt_tRNA-annotation",
                                                         "miRNA_primary_transcript","rRNA"),] 

all.regions2["snoRNA",]=all.regions2["snoRNA",]+all.regions1["snoRNA-annotation",colnames(all.regions2)]
all.regions2["snRNA",]=all.regions2["snRNA",]+all.regions1["snRNA-annotation",colnames(all.regions2)]
all.regions2["tRNA",]=all.regions2["tRNA",]+all.regions1["tRNA-annotation",colnames(all.regions2)]
all.regions2["mt_tRNA",]=all.regions2["mt_tRNA",]+all.regions1["Mt_tRNA-annotation",colnames(all.regions2)]
all.regions2["rRNA",]=all.regions2["rRNA",]+all.regions1["rRNA",colnames(all.regions2)]


all.regions<-rbind(all.regions3,all.regions2)


# if only one sample, take just the first column since second is dup of first
if(length(single.bams)==1){
  all.regions=cbind(all.regions[,1])
  colnames(all.regions)=names(single.bams)
}

## save files
saveRDS(all.regions,paste0(bam_dir,"../counts/",format(Sys.Date(),"%y%m%d"),".biotype.counts.rds"))
write.table(all.regions,paste0(bam_dir,"../counts/",format(Sys.Date(),"%y%m%d"),".biotype.counts.txt"),col.names=NA,row.names=T,quote=F,sep="\t")


#### make biotype plot, with protein-coding and rRNA

keep_biotypes=c("protein_coding","lncRNA","snoRNA","miRNA","mt_tRNA","tRNA","snRNA","piRNA","misc_RNA","rRNA")

if(length(single.bams)>1){
  all.regions2<-all.regions[setdiff(rownames(all.regions),c("non_annotation")),]
  all.regions3<-rbind(all.regions2[keep_biotypes,],other=colSums(all.regions2[setdiff(rownames(all.regions2),keep_biotypes),]))
  all.regions.norm<-100*(t(t(all.regions3)/colSums(all.regions3)))
  regions.frame<-make_data_frame(all.regions.norm,x_label="sample",y_label="biotype",y_levels=c("other",rev(keep_biotypes)))
}

if(length(single.bams)==1){
  all.regions2<-all.regions[setdiff(rownames(all.regions),c("non_annotation")),]
  all.regions3<-c(all.regions2[keep_biotypes],other=sum(all.regions2[setdiff(names(all.regions2),keep_biotypes)]))
  all.regions.norm<-100*(t(t(all.regions3)/sum(all.regions3)))
  regions.frame=data.frame(sample=names(single.bams),biotype=factor(rownames(all.regions.norm),
                                                                    levels=c("other",rev(keep_biotypes))),value=all.regions.norm[,1])
}


colors=c("white","gray","#7A0403FF","#CB2A04FF","#F66B19FF","#FABA39FF","#C7EF34FF","#72FE5EFF","#1AE4B6FF","#4662D7FF","#30123BFF")

num_samples=length(single.bams)
if(num_samples<10){
  fig.width=5
}
if(num_samples>=10 & num_samples<20){
  fig.width=10
}
if(num_samples>=20){
  fig.width=25
}


pdf(paste0(bam_dir,"../figures/",format(Sys.Date(),"%y%m%d"),".all_biotypes.pdf"),height=10,width=fig.width)
ggplot(regions.frame,aes(x=sample,y=value,fill=biotype))+geom_bar(stat="identity",color="black")+ylab("% of mapped reads")+xlab("")+ggtitle("Biotypes of mapped reads")+scale_fill_manual(values=colors)+theme(axis.text.x=element_text(angle=45,hjust=1))
invisible(dev.off())


######


keep_biotypes=c("lncRNA","snoRNA","miRNA","mt_tRNA","tRNA","snRNA","misc_RNA")


if(length(single.bams)>1){
  all.regions2<-all.regions[setdiff(rownames(all.regions),c("non_annotation","protein_coding","rRNA")),]
  all.regions3<-rbind(all.regions2[keep_biotypes,],other=colSums(all.regions2[setdiff(rownames(all.regions2),keep_biotypes),]))
  all.regions.norm<-100*(t(t(all.regions3)/colSums(all.regions3)))
  regions.frame<-make_data_frame(all.regions.norm,x_label="sample",y_label="biotype",y_levels=c("other",rev(keep_biotypes)))
}

if(length(single.bams)==1){
  all.regions2<-all.regions[setdiff(rownames(all.regions),c("non_annotation","protein_coding","rRNA")),]
  all.regions3<-c(all.regions2[keep_biotypes],other=sum(all.regions2[setdiff(names(all.regions2),keep_biotypes)]))
  all.regions.norm<-100*(t(t(all.regions3)/sum(all.regions3)))
  regions.frame=data.frame(sample=names(single.bams),biotype=factor(rownames(all.regions.norm),
                                                                    levels=c("other",rev(keep_biotypes))),value=all.regions.norm[,1])
}

colors<-c("white","#7A0403FF","#DB3A07FF","#FE9B2DFF","#D2E935FF","#62FC6BFF","#1BD0D5FF","#4777EFFF")

pdf(paste0(bam_dir,"../figures/",format(Sys.Date(),"%y%m%d"),".ncRNA_biotypes.pdf"),height=10,width=fig.width)
ggplot(regions.frame,aes(x=sample,y=value,fill=biotype))+geom_bar(stat="identity",color="black")+ylab("% of mapped, non-protein-coding/rRNA reads")+xlab("")+ggtitle("Biotypes of mapped reads, non-protein-coding/rRNA")+scale_fill_manual(values=colors)+theme(axis.text.x=element_text(angle=45,hjust=1))
invisible(dev.off())



