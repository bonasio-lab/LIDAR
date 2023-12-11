
suppressPackageStartupMessages(library(GenomicRanges))
suppressPackageStartupMessages(library(GenomicAlignments))
suppressPackageStartupMessages(library(rtracklayer))
suppressPackageStartupMessages(library(Rsamtools))
suppressPackageStartupMessages(library(parallel))
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(Biostrings))

args <- commandArgs(trailingOnly = TRUE)

bam_dir=args[1]
annotation_dir=args[2]
species=args[3]
single=as.numeric(args[4])




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



if(species=="mmus"){
  tRNA.map=readRDS(paste0(annotation_dir,"/tRNA/GRCm39.tRNA_map.rds"))
  tRNA.seqs=readDNAStringSet(paste0(annotation_dir,"/tRNA/GRCm39.tRNA.fa"))
}
if(species=="hsap"){
  tRNA.map=readRDS(paste0(annotation_dir,"/tRNA/GRCh38.tRNA_map.rds"))
  tRNA.seqs=readDNAStringSet(paste0(annotation_dir,"/tRNA/GRCh38.tRNA.fa"))
}

tRNA.widths=width(tRNA.seqs)
names(tRNA.widths)=sapply(names(tRNA.seqs),get_split," ",1)

tRNA.codons=sapply(names(tRNA.seqs),get_split,"-",2)
names(tRNA.codons)=sapply(names(tRNA.seqs),get_split," ",1)
tRNA.codons<-gsub("iMet","Met",tRNA.codons)

AAs=c("Ala","Arg","Asn","Asp","Cys","Gln","Glu","Gly","His","Ile","Leu","Lys","Met","Phe","Pro","Ser","Thr","Trp","Tyr","Val","SeC")

## for each tRNA, determine the position that will define fragments and halves
pos.tRF3_start=apply(tRNA.map,1,function(x){
  # position 48 = 70th column
  x=x[70:length(x)]
  x=x[!is.na(x)]
  return(min(x))
})

pos.half3_start=apply(tRNA.map,1,function(x){
  # position 35 = 38th column
  x=x[38:length(x)]
  x=x[!is.na(x)]
  return(min(x))
})

pos.tRF5_end=apply(tRNA.map,1,function(x){
  # position 31 = 34th column
  x=x[1:34]
  x=x[!is.na(x)]
  return(max(x))  
})

pos.half5_end=apply(tRNA.map,1,function(x){
  # position 35 = 38th column
  x=x[1:38]
  x=x[!is.na(x)]
  return(max(x))
})

single.files<-list.files(paste0(bam_dir,"tRNA/"),pattern="single.tRNA.bam",full.names=T)
names(single.files)<-gsub(".single.tRNA.bam","",sapply(single.files,get_split,"//",2))
single0<-sapply(single.files,readGAlignments,use.names=T,param=ScanBamParam(tag="RG",what=c("flag")))

if(single==0){
  paired.files<-list.files(paste0(bam_dir,"tRNA/"),pattern="paired.tRNA.bam",full.names=T)
  names(paired.files)<-gsub(".paired.tRNA.bam","",sapply(paired.files,get_split,"//",2))
  paired0<-sapply(paired.files,readGAlignments,use.names=T,param=ScanBamParam(tag=c("RG"),what=c("flag","mrnm","mpos")))
}




single.bam<-sapply(names(single0),function(bam0){
  #print(bam0)
  bam=single0[[bam0]]
  
  # only keep reads that are at least 15 nt
  if(length(bam)==0){
    return(bam)
  }
  bam=bam[mcols(bam)$flag %in% c(0,1024)]
  bam=bam[mcols(bam)$RG=="tRNA_GSNAP"]
  
  
  mcols(bam)$tRNA.width<-tRNA.widths[as.vector(seqnames(bam))]
  mcols(bam)$tRF3_start=pos.tRF3_start[as.vector(seqnames(bam))]
  mcols(bam)$half3_start=pos.half3_start[as.vector(seqnames(bam))]
  
  mcols(bam)$tRF5_end=pos.tRF5_end[as.vector(seqnames(bam))]
  mcols(bam)$half5_end=pos.half5_end[as.vector(seqnames(bam))]
  
  
  #####
  
  bam.start_frag=rep(F,length(bam))
  bam.end_frag=rep(F,length(bam))
  
  bam.full=rep(F,length(bam))
  bam.trf3=rep(F,length(bam))
  bam.trf5=rep(F,length(bam))
  bam.half3=rep(F,length(bam))
  bam.half5=rep(F,length(bam))
  bam.internal=rep(F,length(bam))
  bam.full_close=rep(F,length(bam))
  bam.frag3=rep(F,length(bam))
  bam.frag5=rep(F,length(bam))
  
  bam.close_end_frag=rep(F,length(bam))
  
  names(bam.start_frag)<-names(bam)
  names(bam.end_frag)<-names(bam)
  names(bam.full)<-names(bam)
  names(bam.close_end_frag)=names(bam)
  
  names(bam.trf3)<-names(bam)
  names(bam.trf5)<-names(bam)
  names(bam.half3)<-names(bam)
  names(bam.half5)<-names(bam)
  names(bam.internal)<-names(bam)
  
  names(bam.frag3)<-names(bam)
  names(bam.frag5)<-names(bam)
  names(bam.full_close)<-names(bam)
  
  
  bam.start_frag[start(bam)==1]=T
  bam.end_frag[(mcols(bam)$tRNA.width==end(bam))]=T
  
  end_dist=mcols(bam)$tRNA.width-end(bam)
  bam.close_end_frag[end_dist<=5]=T
  
  bam.full[bam.start_frag & (bam.end_frag)]=T
  
  bam.trf3[bam.close_end_frag==T & start(bam)>=mcols(bam)$tRF3_start]=T
  bam.half3[bam.close_end_frag==T & start(bam)>=mcols(bam)$half3_start]=T
  
  bam.trf5[bam.start_frag==T & end(bam)<=mcols(bam)$tRF5_end]=T
  bam.half5[bam.start_frag==T & end(bam)<=mcols(bam)$half5_end]=T
  
  # anything that is a tRF can't be a half
  bam.half3[bam.trf3]=F
  bam.half5[bam.trf5]=F
  
  bam.internal[start(bam)>=10 & (mcols(bam)$tRNA.width-end(bam))>=10]=T
  
  bam.frag3[bam.trf3 | bam.half3]=T
  bam.frag5[bam.trf5 | bam.half5]=T
  
  bam.full_close[bam.close_end_frag & bam.start_frag]=T
  
  mcols(bam)$start_frag=bam.start_frag
  mcols(bam)$end_frag=bam.end_frag
  mcols(bam)$full=bam.full
  mcols(bam)$close_end_frag=bam.close_end_frag
  
  mcols(bam)$half5=bam.half5
  mcols(bam)$trf5=bam.trf5
  
  mcols(bam)$half3=bam.half3
  mcols(bam)$trf3=bam.trf3
  
  mcols(bam)$internal=bam.internal
  mcols(bam)$frag5=bam.frag5
  mcols(bam)$frag3=bam.frag3
  mcols(bam)$full_close=bam.full_close
  
  
  return(bam)
})


if(single==0){
  paired.bam<-sapply(names(paired0),function(bam0){
    #print(bam0)
    bam=paired0[[bam0]]
    bam=bam[mcols(bam)$flag %in% c(99,147,1123,1171)]
    bam=bam[mcols(bam)$RG=="tRNA_GSNAP"]
    if(length(bam)==0){
      return(bam)
    }
    
    mcols(bam)$tRNA.width=tRNA.widths[as.vector(seqnames(bam))]
    
    mcols(bam)$tRNA.width<-tRNA.widths[as.vector(seqnames(bam))]
    mcols(bam)$tRF3_start=pos.tRF3_start[as.vector(seqnames(bam))]
    mcols(bam)$half3_start=pos.half3_start[as.vector(seqnames(bam))]
    
    mcols(bam)$tRF5_end=pos.tRF5_end[as.vector(seqnames(bam))]
    mcols(bam)$half5_end=pos.half5_end[as.vector(seqnames(bam))]
    
    
    bam.start_frag=rep(F,length(bam)/2)
    bam.end_frag=rep(F,length(bam)/2)
    bam.full=rep(F,length(bam)/2)
    bam.close_end_frag=rep(F,length(bam)/2)
    
    bam.full=rep(F,length(bam)/2)
    
    bam.trf3=rep(F,length(bam)/2)
    bam.trf5=rep(F,length(bam)/2)
    bam.half3=rep(F,length(bam)/2)
    bam.half5=rep(F,length(bam)/2)
    bam.internal=rep(F,length(bam)/2)
    
    bam.full_close=rep(F,length(bam)/2)
    bam.frag3=rep(F,length(bam)/2)
    bam.frag5=rep(F,length(bam)/2)
    
    
    names(bam.start_frag)<-names(bam[mcols(bam)$flag==99])
    names(bam.end_frag)<-names(bam[mcols(bam)$flag==99])
    names(bam.full)<-names(bam[mcols(bam)$flag==99])
    names(bam.close_end_frag)=names(bam[mcols(bam)$flag==99])
    
    
    names(bam.trf3)<-names(bam[mcols(bam)$flag==99])
    names(bam.trf5)<-names(bam[mcols(bam)$flag==99])
    names(bam.half3)<-names(bam[mcols(bam)$flag==99])
    names(bam.half5)<-names(bam[mcols(bam)$flag==99])
    names(bam.internal)<-names(bam[mcols(bam)$flag==99])
    
    names(bam.full_close)<-names(bam[mcols(bam)$flag==99])
    names(bam.frag3)<-names(bam[mcols(bam)$flag==99])
    names(bam.frag5)<-names(bam[mcols(bam)$flag==99])
    
    
    bam.start_frag[names(bam[start(bam)==1 & mcols(bam)$flag==99])]=T
    bam.end_frag[names(bam[(mcols(bam)$tRNA.width==end(bam)) & mcols(bam)$flag==147])]=T
    
    bam.close_end_frag[names(bam[((mcols(bam)$tRNA.width-end(bam)))<=5 & mcols(bam)$flag==147])]=T
    
    bam.full[bam.start_frag & (bam.end_frag)]=T
    
    
    bam.trf3.names=intersect(names(bam[start(bam)>=mcols(bam)$tRF3_start & mcols(bam)$flag==99]),
                             names(bam.close_end_frag[bam.close_end_frag]))
    
    bam.half3.names=intersect(names(bam[start(bam)>=mcols(bam)$half3_start & mcols(bam)$flag==99]),
                              names(bam.close_end_frag[bam.close_end_frag]))
    
    bam.frag3.names=c(bam.trf3.names,bam.half3.names)
    
    bam.trf5.names=intersect(names(bam[end(bam)<=mcols(bam)$tRF5_end & mcols(bam)$flag==147]),
                             names(bam.start_frag[bam.start_frag]))
    
    bam.half5.names=intersect(names(bam[end(bam)<=mcols(bam)$half5_end & mcols(bam)$flag==147]),
                              names(bam.start_frag[bam.start_frag]))
    bam.frag5.names=c(bam.trf5.names,bam.half5.names)
    
    bam.full_close.names=intersect(names(bam.close_end_frag[bam.close_end_frag]),names(bam.start_frag[bam.start_frag]))
    
    bam.trf3[bam.trf3.names]=T
    bam.half3[setdiff(bam.half3.names,bam.trf3.names)]=T
    bam.frag3[bam.frag3.names]=T
    
    bam.trf5[bam.trf5.names]=T
    bam.half5[setdiff(bam.half5.names,bam.trf5.names)]=T
    bam.frag5[bam.frag5.names]=T
    
    bam.internal.R1=names(bam[start(bam)>=10 & mcols(bam)$flag==99])
    bam.internal.R2=names(bam[mcols(bam)$tRNA.width-end(bam)>=10 & mcols(bam)$flag==147])
    
    bam.internal[intersect(bam.internal.R1,bam.internal.R2)]=T
    bam.full_close[bam.full_close.names]=T
    
    
    
    bam.pair=makeGAlignmentPairs(bam,use.names=T,use.mcols = T)
    
    mcols(bam.pair)$start_frag=bam.start_frag[names(bam.pair)]
    mcols(bam.pair)$end_frag=bam.end_frag[names(bam.pair)]
    mcols(bam.pair)$full=bam.full[names(bam.pair)]
    mcols(bam.pair)$close_end_frag=bam.close_end_frag[names(bam.pair)]
    mcols(bam.pair)$half5=bam.half5[names(bam.pair)]
    mcols(bam.pair)$trf5=bam.trf5[names(bam.pair)]
    mcols(bam.pair)$half3=bam.half3[names(bam.pair)]
    mcols(bam.pair)$trf3=bam.trf3[names(bam.pair)]
    mcols(bam.pair)$internal=bam.internal[names(bam.pair)]
    mcols(bam.pair)$frag3=bam.frag3[names(bam.pair)]
    mcols(bam.pair)$frag5=bam.frag5[names(bam.pair)]
    mcols(bam.pair)$full_close=bam.full_close[names(bam.pair)]
    
    return(bam.pair)
  })
}



tRNA.split.counts<-lapply(names(single.bam),function(x){
  #print(x)
  x.single.bam=single.bam[[x]]
  x.single.bam.split<-list(trf5=x.single.bam[mcols(x.single.bam)$trf5],
                           half5=x.single.bam[mcols(x.single.bam)$half5],
                           trf3=x.single.bam[mcols(x.single.bam)$trf3],
                           half3=x.single.bam[mcols(x.single.bam)$half3],
                           internal=x.single.bam[mcols(x.single.bam)$internal],
                           full_close=x.single.bam[mcols(x.single.bam)$full_close])
  
  x.single.bam.seqnames=sapply(x.single.bam.split,function(y){return(as.vector(seqnames(y)))})
  x.single.counts=sapply(x.single.bam.seqnames,function(y){
    sapply(names(tRNA.codons),function(z){
      return(sum(y==z))
    })
  })
  
  if(single==0){
    x.paired.bam=paired.bam[[x]]
    x.paired.bam.split<-list(trf5=x.paired.bam[mcols(x.paired.bam)$trf5],
                             half5=x.paired.bam[mcols(x.paired.bam)$half5],
                             trf3=x.paired.bam[mcols(x.paired.bam)$trf3],
                             half3=x.paired.bam[mcols(x.paired.bam)$half3],
                             internal=x.paired.bam[mcols(x.paired.bam)$internal],
                             full_close=x.paired.bam[mcols(x.paired.bam)$full_close])
    
    x.paired.bam.seqnames=sapply(x.paired.bam.split,function(y){return(as.vector(seqnames(y)))})
    x.paired.counts=sapply(x.paired.bam.seqnames,function(y){
      sapply(names(tRNA.codons),function(z){
        return(sum(y==z))
      })
    })
    x.counts=x.single.counts+x.paired.counts
  }else{
    x.counts=x.single.counts
  }
  return(x.counts)
  
})

names(tRNA.split.counts)<-names(single.bam)



tRNA.fragment.freq.cyto.frame<-do.call(rbind,lapply(names(tRNA.split.counts),function(x0){
  x=tRNA.split.counts[[x0]]
  x=x[!grepl("mito",rownames(x)),]
  x.sums=colSums(x[,colnames(x)!="full_close"])
  x.frame<-data.frame(sample=x0,
                      fragment=factor(names(x.sums),levels=c("half5","trf5","internal","half3","trf3")),
                      count=x.sums,
                      freq=100*(x.sums/sum(x.sums)))
  return(x.frame)
}))



tRNA.fragment.freq.mito.frame<-do.call(rbind,lapply(names(tRNA.split.counts),function(x0){
  x=tRNA.split.counts[[x0]]
  x=x[grepl("mito",rownames(x)),]
  x.sums=colSums(x[,colnames(x)!="full_close"])
  x.frame<-data.frame(sample=x0,
                      fragment=factor(names(x.sums),levels=c("half5","trf5","internal","half3","trf3")),
                      count=x.sums,
                      freq=100*(x.sums/sum(x.sums)))
  return(x.frame)
}))



num_samples=length(single.bam)
if(num_samples<10){
  fig.width=5
}
if(num_samples>=10 & num_samples<20){
  fig.width=10
}
if(num_samples>=20){
  fig.width=25
}


pdf(paste0(bam_dir,"../figures/",format(Sys.Date(),"%y%m%d"),".tRNA_fragment_distribution_cyto.pdf"),height=10,width=fig.width)
ggplot(tRNA.fragment.freq.cyto.frame,aes(x=sample,y=freq,fill=fragment))+geom_bar(stat="identity")+scale_fill_manual(values=c("gold1","gold3","gray","pink1","pink3"))+theme(axis.text.x=element_text(angle=45,hjust=1))+ylab("% of fragment reads")+xlab("")
invisible(dev.off())

pdf(paste0(bam_dir,"../figures/",format(Sys.Date(),"%y%m%d"),".tRNA_fragment_distribution_mito.pdf"),height=10,width=fig.width)
ggplot(tRNA.fragment.freq.mito.frame,aes(x=sample,y=freq,fill=fragment))+geom_bar(stat="identity")+scale_fill_manual(values=c("gold1","gold3","gray","pink1","pink3"))+theme(axis.text.x=element_text(angle=45,hjust=1))+ylab("% of fragment reads")+xlab("")
invisible(dev.off())





##############


vals=15:95

tRNA.fragment.width.cyto.frame<-do.call(rbind,lapply(names(single.bam),function(x){
  x.single.bam=single.bam[[x]]
  x.single.bam=x.single.bam[!grepl("mito",as.vector(seqnames(x.single.bam)))]

  if(single==0){
    x.paired.bam=GRanges(paired.bam[[x]])
    x.paired.bam=x.paired.bam[!grepl("mito",as.vector(seqnames(x.paired.bam)))]
    x.widths=c(width(x.single.bam),width(x.paired.bam))
  }else{
    x.widths=width(x.single.bam)
  }
  
  x.widths.count=sapply(vals,function(y){return(sum(x.widths==y))})
  x.frame=data.frame(sample=x,width=vals,count=x.widths.count,freq=100*(x.widths.count/sum(x.widths.count)))
  return(x.frame)
}))



######


num_samples=length(single.bam)
if(num_samples<10){
  fig.width=8
  fig.height=5
}
if(num_samples>=10 & num_samples<20){
  fig.width=8
  fig.height=10
}
if(num_samples>=20){
  fig.width=8
  fig.height=20
}


tRNA.fragment.width.mito.frame<-do.call(rbind,lapply(names(single.bam),function(x){
  x.single.bam=single.bam[[x]]
  x.single.bam=x.single.bam[grepl("mito",as.vector(seqnames(x.single.bam)))]
  
  if(single==0){
    x.paired.bam=GRanges(paired.bam[[x]])
    x.paired.bam=x.paired.bam[grepl("mito",as.vector(seqnames(x.paired.bam)))]
    x.widths=c(width(x.single.bam),width(x.paired.bam))
  }else{
    x.widths=width(x.single.bam)
  }
  
  x.widths.count=sapply(vals,function(y){return(sum(x.widths==y))})
  x.frame=data.frame(sample=x,width=vals,count=x.widths.count,freq=100*(x.widths.count/sum(x.widths.count)))
  return(x.frame)
}))


pdf(paste0(bam_dir,"../figures/",format(Sys.Date(),"%y%m%d"),".tRNA_cyto_fragment_lengths.pdf"),height=fig.height,width=fig.width)
ggplot(tRNA.fragment.width.cyto.frame,aes(x=width,y=freq,color=sample,group=sample))+geom_line()+facet_wrap(~sample,ncol=2)+scale_x_continuous("Fragment length",breaks=c(48,63,70,78))+ggtitle("Fragment lengths of reads aligning to tRNA")+ylab("% of reads")+coord_cartesian(xlim=c(15,95))
invisible(dev.off())


pdf(paste0(bam_dir,"../figures/",format(Sys.Date(),"%y%m%d"),".tRNA_mito_fragment_lengths.pdf"),height=fig.height,width=fig.width)
ggplot(tRNA.fragment.width.mito.frame,aes(x=width,y=freq,color=sample,group=sample))+geom_line()+facet_wrap(~sample,ncol=2)+scale_x_continuous("Fragment length",breaks=c(48,63,70,78))+ggtitle("Fragment lengths of reads aligning to tRNA")+ylab("% of reads")+coord_cartesian(xlim=c(15,95))
invisible(dev.off())



###########

stats.files=list.files(paste0(bam_dir,"../logs/"),pattern="stats$",full.names = T)
stats=lapply(stats.files,read.delim,header=F)
libSizes<-sapply(stats,function(x){
  x.sum=sum(x[x$V2=="map",]$V5)
  names(x.sum)=x$V1[1]
  return(x.sum)
})

canonical_coverage=lapply(names(single.bam),function(bam0){
  # print(bam0)
  x.single.bam=single.bam[[bam0]]
  x.single.coverage=sapply(coverage(x.single.bam)[names(tRNA.codons)],as.vector)
  
  
  if(single==0){
    x.paired.bam=GRanges(paired.bam[[bam0]])
    x.paired.coverage=sapply(coverage(x.paired.bam)[names(tRNA.codons)],as.vector)
  }

  coverage.map=t(sapply(names(x.single.coverage),function(x0){
    x.single=x.single.coverage[[x0]]
    
    if(single==0){
      x.paired=x.paired.coverage[[x0]]
      x=x.single+x.paired
    }else{
      x=x.single
    }

    names(x)=1:length(x)
    
    x.map=tRNA.map[x0,]
    
    x.cov.map=x[x.map]
    names(x.cov.map)=names(x.map)
    
    ## smoothing
    # for NA, find the value to the left and to the right that are not NA
    x.cov.map.NA=which(is.na(x.cov.map))
    x.cov.map.nonNA=which(!is.na(x.cov.map))
    
    x.cov.map.fix=x.cov.map
    for(i in 1:length(x.cov.map.NA)){
      i.select=x.cov.map.NA[i]
      i.left0=x.cov.map.nonNA[x.cov.map.nonNA<i.select]
      if(length(i.left0)==0){
        i.left=-1
      }else{
        i.left=max(i.left0)
      }
      i.right0=x.cov.map.nonNA[x.cov.map.nonNA>i.select]
      if(length(i.right0)==0){
        i.right=-1
      }else{
        i.right=min(i.right0)
      }
      if(i.left!=-1 & i.right!=-1){
        i.val=mean(x.cov.map[i.left],x.cov.map[i.right])
      }
      if(i.left!=-1 & i.right==-1){
        i.val=x.cov.map[i.left]
      }
      if(i.left==-1 & i.right!=-1){
        i.val=x.cov.map[i.right]
      }
      x.cov.map.fix[names(i.select)]=i.val
    }
    
    return(x.cov.map.fix)
  }))
  coverage.map[is.na(coverage.map)]=0
  
  # positions before var region
  var.start=min(which(grepl("e",colnames(tRNA.map))==T))
  var.end=max(which(grepl("e",colnames(tRNA.map))==T))
  
  
  coverage.map1=coverage.map[,1:(var.start-1)]
  coverage.map2=coverage.map[,var.start:var.end]
  coverage.map3=coverage.map[,(var.end+1):ncol(tRNA.map)]
  
  coverage.map.collapse=cbind(coverage.map1,var=rowMeans(coverage.map2),coverage.map3)
  
  coverage.map.collapse.norm=coverage.map.collapse/(libSizes[bam0]/1e6)
  
  return(coverage.map.collapse.norm)
  
})
names(canonical_coverage)=names(single.bam)

canonical_coverage.by_codon.cyto=lapply(canonical_coverage,function(x){
  x.by_codon=t(sapply(AAs,function(y){
    select.tRNA=names(tRNA.codons)[grepl(y,names(tRNA.codons))]
    select.tRNA=select.tRNA[!grepl("mito",select.tRNA)]
    x.select=x[select.tRNA,]
    if(length(select.tRNA)==1){
      return(x.select)
    }
    return(colSums(x.select))
  }))
  return(x.by_codon)
})

canonical_coverage.by_codon.mito=lapply(canonical_coverage,function(x){
  x.by_codon=t(sapply(AAs,function(y){
    select.tRNA=names(tRNA.codons)[grepl(y,names(tRNA.codons))]
    select.tRNA=select.tRNA[grepl("mito",select.tRNA)]
    x.select=x[select.tRNA,]
    if(length(select.tRNA)==1){
      return(x.select)
    }
    return(colSums(x.select))
  }))
  return(x.by_codon)
})

names(canonical_coverage.by_codon.cyto)=names(canonical_coverage)
names(canonical_coverage.by_codon.mito)=names(canonical_coverage)

coverage.by_codon.cyto.frame<-do.call(rbind,lapply(names(canonical_coverage.by_codon.cyto),function(x0){
  x=canonical_coverage.by_codon.cyto[[x0]]
  x.frame=make_data_frame(x,x_label="position",y_label="AA",value_label="RPM")
  x.frame$sample=x0
  return(x.frame)
}))

coverage.by_codon.mito.frame<-do.call(rbind,lapply(names(canonical_coverage.by_codon.mito),function(x0){
  x=canonical_coverage.by_codon.mito[[x0]]
  x.frame=make_data_frame(x,x_label="position",y_label="AA",value_label="RPM")
  x.frame$sample=x0
  return(x.frame)
}))

colors=c("#30123BFF","#3E378FFF","#455BCDFF","#477CF3FF","#3E9BFEFF","#28BBECFF","#18D6CBFF","#20EAABFF","#46F884FF","#78FE5AFF","#A2FC3CFF","#C4F134FF","#E1DD37FF","#F6C33AFF","#FEA632FF","#FB8022FF","#F05B12FF","#DD3D08FF","#C42503FF","#A31301FF","#7A0403FF")

pdf(paste0(bam_dir,"../figures/",format(Sys.Date(),"%y%m%d"),".tRNA_coverage_by_codon_cyto.pdf"),height=fig.height,width=fig.width)
ggplot(coverage.by_codon.cyto.frame,aes(x=position,y=RPM,fill=AA,group=AA))+geom_bar(stat="identity",width=1)+facet_wrap(~sample,ncol=2)+guides(guide_legend(fill=guide_legend(ncol=2,title="")))+theme(axis.text.x = element_blank(),axis.ticks.x = element_blank())+scale_fill_manual(values=colors)+ylab("% of reads mapping to tRNA")+xlab("Position along tRNA")+ggtitle("Codon distribution by position along cyto tRNA")+theme(legend.title=element_blank())
invisible(dev.off())

pdf(paste0(bam_dir,"../figures/",format(Sys.Date(),"%y%m%d"),".tRNA_coverage_by_codon_mito.pdf"),height=fig.height,width=fig.width)
ggplot(coverage.by_codon.mito.frame,aes(x=position,y=RPM,fill=AA,group=AA))+geom_bar(stat="identity",width=1)+facet_wrap(~sample,ncol=2)+guides(guide_legend(fill=guide_legend(ncol=2,title="")))+theme(axis.text.x = element_blank(),axis.ticks.x = element_blank())+scale_fill_manual(values=colors)+ylab("% of reads mapping to tRNA")+xlab("Position along tRNA")+ggtitle("Codon distribution by position along mito tRNA")+theme(legend.title=element_blank())
invisible(dev.off())




