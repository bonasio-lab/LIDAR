library(GenomicAlignments)
library(GenomicRanges)
library(rtracklayer)
library(pheatmap)
library(Rsamtools)
library(parallel)
library(DESeq2)
library(ggplot2)


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


```


#import bams and set bam parameters
```{r,eval=F}

all.info<-readRDS("/home/emily/data/EJ25/code/sample_info/sample.info.rds")


single.bams<-all.info$singleBam
names(single.bams)<-rownames(all.info)

paired.bams<-all.info$pairedBam
names(paired.bams)<-rownames(all.info)


#use only primary alignments
param=ScanBamParam(flag=scanBamFlag(isSecondaryAlignment = FALSE))
stranded=T
#use only primary alignments and filter out duplicates
#param=ScanBamParam(flag=scanBamFlag(isSecondaryAlignment = FALSE,isDuplicate = FALSE))


```


## split the genome into biotypes
note: a GRangesList needs for each of its elements to have the same metadata or it will behave strangely, so here I just wiped the metadata for each of the components. Also going to count only the annotation/miRNA, then count the other scaffolds separately, as this will make it quicker. 
```{r}
miRNA.annotation.gr<-as(import.gff("/home/emily/data/genome_annotation/annotation/GRCm39/miRNA/mmu.GRCm39.gff3"),"GRanges")
miRNA.annotation.gr<-miRNA.annotation.gr[miRNA.annotation.gr$type=="miRNA"]
mcols(miRNA.annotation.gr)<-mcols(miRNA.annotation.gr)[,c("type","Name")]
colnames(mcols(miRNA.annotation.gr))<-c("gene_biotype","gene")

annotation.gr<-import.gff("/home/emily/data/genome_annotation/annotation/GRCm39/refseq/GRCm39.converted.gff")

annotation.gene<-annotation.gr[annotation.gr$type=="gene" | annotation.gr$type=="pseudogene"]
annotation.biotypes=annotation.gene$gene_biotype
names(annotation.biotypes)=annotation.gene$gene

annotation.exons<-annotation.gr[annotation.gr$type=="exon"]
annotation.exons$gene_biotype=annotation.biotypes[annotation.exons$gene]

annotation.exons$gene_biotype<-gsub("snRNA","snRNA-annotation",annotation.exons$gene_biotype)
annotation.exons$gene_biotype<-gsub("snoRNA","snoRNA-annotation",annotation.exons$gene_biotype)
annotation.exons$gene_biotype<-gsub("tRNA","tRNA-annotation",annotation.exons$gene_biotype)

annotation.exons<-annotation.exons[annotation.exons$gene_biotype!="miRNA"]

mcols(annotation.exons)<-mcols(annotation.exons)[,c("gene_biotype","gene")]

# now fix the biotypes so that all the Mt_tRNA are marked properly
annotation.exons[as.vector(seqnames(annotation.exons))=="chrM" & annotation.exons$gene_biotype=="tRNA-annotation"]$gene_biotype="Mt_tRNA-annotation"



# #also just get rid of any transcripts in the gencode that overlaps with an annotated miRNA from miRbase - this seems necessary, for example see Gm37327 (chr12 113432809-113432837), it's annotated as Ig_D gene but overlaps with miRNA. 
annotation.miRNA.overlaps=findOverlaps(annotation.exons,miRNA.annotation.gr)
annotation.exons<-annotation.exons[setdiff(1:length(annotation.exons),queryHits(annotation.miRNA.overlaps))]

all.gr<-c(annotation.exons,miRNA.annotation.gr)

# split by biotype and then overlap them - find ambiguous regions
annotation.exons.by_biotype=split(all.gr,all.gr$gene_biotype)

annotation.exons.by_biotype.reduced=sapply(annotation.exons.by_biotype,function(x){
  mcols(x)<-NULL
  return(reduce(x))
})

all.reduced=unlist(GRangesList(annotation.exons.by_biotype.reduced))
all.reduced.plus=all.reduced[as.vector(strand(all.reduced))=="+"]
all.reduced.minus=all.reduced[as.vector(strand(all.reduced))=="-"]

all.reduced.plus.coverage=coverage(all.reduced.plus)
all.reduced.minus.coverage=coverage(all.reduced.minus)

ambiguous.regions.plus<-sapply(names(all.reduced.plus.coverage),function(x0){
  x=all.reduced.plus.coverage[[x0]]
  x.select=which(x>1)
  if(length(x.select)==0){
    return()
  }
  x.select.gr=reduce(GRanges(seqnames=x0,
                             ranges=IRanges(start=x.select,
                                            end=x.select),
                             strand="+"))
  return(x.select.gr)
  
})

ambiguous.regions.plus=unlist(GRangesList(ambiguous.regions.plus[!sapply(ambiguous.regions.plus,is.null)]))


ambiguous.regions.minus<-sapply(names(all.reduced.minus.coverage),function(x0){
  x=all.reduced.minus.coverage[[x0]]
  x.select=which(x>1)
  if(length(x.select)==0){
    return()
  }
  x.select.gr=reduce(GRanges(seqnames=x0,
                             ranges=IRanges(start=x.select,
                                            end=x.select),
                             strand="+"))
  return(x.select.gr)
  
})

ambiguous.regions.minus=unlist(GRangesList(ambiguous.regions.minus[!sapply(ambiguous.regions.minus,is.null)]))


ambiguous.regions=c(ambiguous.regions.plus,ambiguous.regions.minus)


# find ambiguous regions that involve protein-coding - this is ~ half
ambiguous.pc_regions<-subsetByOverlaps(ambiguous.regions,annotation.exons.by_biotype.reduced$protein_coding,ignore.strand=F)

annotation.exons.by_biotype.reduced.unique=sapply(annotation.exons.by_biotype.reduced,function(x){
  return(setdiff(x,ambiguous.regions,ignore.strand=F))
})

regions.biotype=c(annotation.exons.by_biotype.reduced.unique,
                  ambiguous_non_coding=subsetByOverlaps(ambiguous.regions,annotation.exons.by_biotype.reduced$protein_coding,ignore.strand=F,invert=T),
                  ambiguous_coding=subsetByOverlaps(ambiguous.regions,annotation.exons.by_biotype.reduced$protein_coding,ignore.strand=F,invert=F))

regions.biotype=GRangesList(regions.biotype)

saveRDS(regions.biotype,"/home/emily/data/EJ25/code/regions.biotype_refseq3.rds")
```

```{r}
regions.biotype<-readRDS("/home/emily/data/EJ25/code/regions.biotype_refseq3.rds")

```

```{r}
count_reads<-function(reads,annotation,stranded=stranded,singleEnd,param=param){
  return(assay(summarizeOverlaps(annotation,reads,ignore.strand=!stranded,singleEnd=singleEnd,param=param,mode="IntersectionNotEmpty")))
}



# need to split the bams into chunks to run the counting - if you don't the parlapply cl will crash with mysterious "error writing to connection." From my tests it can only handle 20-30 files at once, sometimes less. I think it's a memory thing, I have been trying to fix for 6 months. Changing the number of cl doesn't help. This was the simplest way I could get around it. 
get_chunk<- function(x,n) split(x, cut(seq_along(x), n, labels = FALSE)) 
chunks=get_chunk(1:length(single.bams),length(single.bams)/20)

cl<-makeCluster(20)
clusterExport(cl, list("summarizeOverlaps","assay"))

# single.counts<-do.call("cbind",parLapply(cl=cl,single.bams,count_reads,annotation=regions.biotype,stranded=stranded,singleEnd=T,param=param))
# 
# paired.counts<-do.call("cbind",parLapply(cl=cl,paired.bams,count_reads,annotation=regions.biotype,stranded=stranded,singleEnd=F,param=param))

single.counts0<-lapply(chunks,function(x){
  print(x)
  return(do.call("cbind",parLapply(cl=cl,single.bams[x],count_reads,annotation=regions.biotype,stranded=stranded,singleEnd=T,param=param)))
})
single.counts<-do.call(cbind,single.counts0)

paired.counts0<-lapply(chunks,function(x){
  print(x)
  return(do.call("cbind",parLapply(cl=cl,paired.bams[x],count_reads,annotation=regions.biotype,stranded=stranded,singleEnd=F,param=param)))
})
paired.counts<-do.call(cbind,paired.counts0)

stopCluster(cl)


rownames(single.counts)=names(regions.biotype)
colnames(single.counts)=names(single.bams)

rownames(paired.counts)=names(regions.biotype)
colnames(paired.counts)=names(paired.bams)

all.counts1=single.counts+paired.counts[,colnames(single.counts)]


```


## counts on the snRNA, snoRNA, piRNA, tRNA, rRNA (separate scaffolds)
```{r}
paired.idxstats0<-list.files("/home/emily/data/EJ25/data/bam/current/filt15/idxstats/",pattern="paired.idxstats",full.names=T)
single.idxstats0<-list.files("/home/emily/data/EJ25/data/bam/current/filt15/idxstats/",pattern="single.idxstats",full.names=T)


names(paired.idxstats0)<-gsub(".paired.idxstats","",sapply(paired.idxstats0,get_split,"//",2))
names(single.idxstats0)<-gsub(".single.idxstats","",sapply(single.idxstats0,get_split,"//",2))

paired.idxstats<-lapply(paired.idxstats0,read.delim,header=F,row.names=1)
single.idxstats<-lapply(single.idxstats0,read.delim,header=F,row.names=1)

# don't include E. coli chr
chrs<-setdiff(rownames(paired.idxstats[[1]]),"NC_000913.3")

paired.idxstats.counts<-sapply(paired.idxstats,function(x){return(x[chrs,]$V3/2)})
rownames(paired.idxstats.counts)<-chrs

single.idxstats.counts<-sapply(single.idxstats,function(x){return(x[chrs,]$V3)})
rownames(single.idxstats.counts)<-chrs

piRNA.names<-read.delim("/home/emily/data/genome_annotation/genomes/bowtie_indexes/EJ25_4P/files_piRNA/piRNA_only.unique.reduced.GRCm39.gold.fa.fai",header=F)$V1
snoRNA.names<-read.delim("/home/emily/data/genome_annotation/genomes/bowtie_indexes/EJ25_4P/files_snoRNA/snoRNA.fa.fai",header=F)$V1
snRNA.names<-read.delim("/home/emily/data/genome_annotation/genomes/bowtie_indexes/EJ25_4P/files_snRNA/snRNA.fa.fai",header=F)$V1
tRNA.names<-read.delim("/home/emily/data/genome_annotation/annotation/GRCm39/tRNA_mim/mimseq_tRNA.fa.fai",header=F)$V1

ncRNA.genes=list(piRNA=piRNA.names,snoRNA=snoRNA.names,snRNA=snRNA.names,tRNA=tRNA.names[!grepl("mito",tRNA.names)],mt_tRNA=tRNA.names[grepl("mito",tRNA.names)],rRNA="BK000964.3")

all.counts2<-t(sapply(ncRNA.genes,function(x){
  select.single=single.idxstats.counts[x,]
  select.paired=paired.idxstats.counts[x,]
  if(length(x)==1){
    return(select.single+select.paired[names(select.single)])
  }
  return(colSums(select.single)+colSums(select.paired[,colnames(select.single)]))
}))


```


```{r}
# add the ncRNA counts from bowtie and genome. The # of tRNA and Mt_tRNA is extremely small, maybe just reads that map over this location to protein-coding or something? just remove them. There are also no reads for miRNA_primary_annotation so remove it to avoid confusion. 
all.counts3=all.counts1[!rownames(all.counts1) %in% c("snRNA-annotation","snoRNA-annotation",
                                                      "tRNA-annotation","Mt_tRNA-annotation",
                                                      "miRNA_primary_transcript","rRNA"),] 

all.counts2["snoRNA",]=all.counts2["snoRNA",]+all.counts1["snoRNA-annotation",colnames(all.counts2)]
all.counts2["snRNA",]=all.counts2["snRNA",]+all.counts1["snRNA-annotation",colnames(all.counts2)]
all.counts2["tRNA",]=all.counts2["tRNA",]+all.counts1["tRNA-annotation",colnames(all.counts2)]
all.counts2["mt_tRNA",]=all.counts2["mt_tRNA",]+all.counts1["Mt_tRNA-annotation",colnames(all.counts2)]
all.counts2["rRNA",]=all.counts2["rRNA",]+all.counts1["rRNA",colnames(all.counts2)]


all.counts<-rbind(all.counts3,all.counts2)

all.counts=rbind(all.counts,non_annotation=all.info[colnames(all.counts),]$bamSize-colSums(all.counts))

saveRDS(all.counts,paste0("/home/emily/data/EJ25/data/counts/",format(Sys.Date(),"%y%m%d"),".region.counts.refseq.rds"))
```
