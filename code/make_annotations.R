######## annotation for gene counts
annotation.file="/home/emily/data/EJ25/code/publish/annotation/GRCm39.gff"
annotation.miRNA.file="/home/emily/data/EJ25/code/publish/annotation/GRCm39.miRNA.gff3"

annotation.gr<-import.gff(annotation.file)
annotation.gene<-annotation.gr[annotation.gr$type=="gene" | annotation.gr$type=="pseudogene"]
annotation.biotypes=annotation.gene$gene_biotype
names(annotation.biotypes)=annotation.gene$gene

annotation.exons<-annotation.gr[annotation.gr$type=="exon"]
annotation.exons$gene_biotype=annotation.biotypes[annotation.exons$gene]

# remove miRNA and tRNA - miRNA will be added from miRbase and tRNA were included as extra scaffolds
annotation.exons<-annotation.exons[!annotation.exons$gene_biotype %in% c("miRNA","tRNA")]

mcols(annotation.exons)<-mcols(annotation.exons)[,c("gene_biotype","gene")]

miRNA.annotation.gr<-import.gff(annotation.miRNA.file)

# remove primary transcripts
miRNA.annotation.gr<-miRNA.annotation.gr[miRNA.annotation.gr$type=="miRNA"]
mcols(miRNA.annotation.gr)<-mcols(miRNA.annotation.gr)[,c("type","Name")]
colnames(mcols(miRNA.annotation.gr))<-c("gene_biotype","gene")

# combine main annotation and miRNA
all.gr<-c(annotation.exons,miRNA.annotation.gr)

saveRDS(all.gr,"/home/emily/data/EJ25/code/publish/annotation/GRCm39.genes.rds")


######

annotation.exons2<-annotation.gr[annotation.gr$type=="exon"]

annotation.exons2$gene_biotype=annotation.biotypes[annotation.exons2$gene]

annotation.exons2$gene_biotype<-gsub("snRNA","snRNA-annotation",annotation.exons2$gene_biotype)
annotation.exons2$gene_biotype<-gsub("snoRNA","snoRNA-annotation",annotation.exons2$gene_biotype)
annotation.exons2$gene_biotype<-gsub("tRNA","tRNA-annotation",annotation.exons2$gene_biotype)

annotation.exons2<-annotation.exons2[annotation.exons2$gene_biotype!="miRNA"]

mcols(annotation.exons2)<-mcols(annotation.exons2)[,c("gene_biotype","gene")]

# now fix the biotypes so that all the Mt_tRNA are marked properly
annotation.exons2[as.vector(seqnames(annotation.exons2))=="chrM" & annotation.exons2$gene_biotype=="tRNA-annotation"]$gene_biotype="Mt_tRNA-annotation"


# #also just get rid of any transcripts in the gencode that overlaps with an annotated miRNA from miRbase - this seems necessary, for example see Gm37327 (chr12 113432809-113432837), it's annotated as Ig_D gene but overlaps with miRNA. 
annotation.miRNA.overlaps=findOverlaps(annotation.exons2,miRNA.annotation.gr)
annotation.exons2<-annotation.exons2[setdiff(1:length(annotation.exons2),queryHits(annotation.miRNA.overlaps))]

all.gr<-c(annotation.exons2,miRNA.annotation.gr)

# split by biotype and then overlap them - find ambiguous regions
annotation.exons2.by_biotype=split(all.gr,all.gr$gene_biotype)

annotation.exons2.by_biotype.reduced=sapply(annotation.exons2.by_biotype,function(x){
  mcols(x)<-NULL
  return(reduce(x))
})

all.reduced=unlist(GRangesList(annotation.exons2.by_biotype.reduced))
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


ambiguous.pc_regions<-subsetByOverlaps(ambiguous.regions,annotation.exons2.by_biotype.reduced$protein_coding,ignore.strand=F)

annotation.exons2.by_biotype.reduced.unique=sapply(annotation.exons2.by_biotype.reduced,function(x){
  return(setdiff(x,ambiguous.regions,ignore.strand=F))
})

regions.biotype=c(annotation.exons2.by_biotype.reduced.unique,
                  ambiguous_non_coding=subsetByOverlaps(ambiguous.regions,annotation.exons2.by_biotype.reduced$protein_coding,ignore.strand=F,invert=T),
                  ambiguous_coding=subsetByOverlaps(ambiguous.regions,annotation.exons2.by_biotype.reduced$protein_coding,ignore.strand=F,invert=F))

regions.biotype=GRangesList(regions.biotype)

saveRDS(regions.biotype,"/home/emily/data/EJ25/code/publish/annotation/GRCm39.regions.rds")
