# library(ggplot2)
library(RColorBrewer)
library("reshape2")
#library("Gviz")
#library("ggbio")
library("GenomicRanges")
library(rtracklayer)
library(biovizBase)
library(gridExtra)
#library(vcfR)
library(VariantAnnotation)

library(GenomicFeatures)
makeTxDbFromGFF(gff_virus)
#ref_fa=FaFile(paste0(basedir,"/References/viruses_short.fasta"))
txdb=makeTxDbFromGFF(paste0(basedir,"/References/viruses_short.gff"),format="auto",dataSource="ensemble",
                     chrominfo =  Seqinfo(seqnames = c("L","S"),seqlength=c(7229,3377), isCircular=c(FALSE, FALSE)),taxonomyId=11623)

# with VariantAnnotation
#coding <- predictCoding(vcf, txdb ,seqSource = ref_fa)



get_afs_annos  <- function(vcf,region_idx,sample) {
  af_list=geno(vcf)[["AF"]][region_idx,sample]
  annos=vcfInfo(vcf)[["ANN"]][region_idx]
  xpos=1:length(af_list)
  afs = matrix(unlist(sapply(xpos, 
                             function (x) unlist( 
                               sapply( 1:length(af_list[[ x ]]), 
                                       function (y) return(c(x,y,af_list[[ x ]][y]))  )) 
  )
  )
  , ncol=3,byrow=TRUE)
  colnames(afs)=c("xpos","ALL","AF")
  afs=as.data.frame(afs)
  afs=afs[! is.na(afs[,3]),]
  if (nrow(afs) == 0) { 
    return (list())
    }
  afs$ANNO=sapply(1:nrow(afs), function (x) annos[[ afs$xpos[x] ]][ afs$ALL[x] ] )
  afs$ANNO[is.na(afs$ANNO)]="-"
  afs$ANNO[grep("synonymous_variant",afs$ANNO)]="synonymous_variant"
  afs$ANNO[grep("missense_variant",afs$ANNO)]="missense_variant"
  afs$ANNO[grep("frameshift_variant",afs$ANNO)]="frameshift_variant"
  afs$ANNO[grep("stop_gained",afs$ANNO)]="stop_gained"
  afs$ANNO[grep("stop_lost",afs$ANNO)]="stop_lost"
  afs$ANNO[grep("start_lost",afs$ANNO)]="start_lost"
  afs$ANNO[grep("start_gained",afs$ANNO)]="start_gained"
  afs$ANNO=factor(afs$ANNO,levels=c("-","synonymous_variant","missense_variant","stop_lost","stop_gained","start_lost","inframe_deletion","frameshift_variant","disruptive_inframe_deletion"))
  afs$ANNO[is.na(afs$ANNO)]="-"
  return(afs)
}

plot_AFs_depths <- function(vcf,chrom,sample,annos=TRUE,depths=TRUE, chrom_len = FALSE,
                            cols=setNames(c("black","green","red","orangered","orange","orangered4","pink","purple","purple4"),
                                          c("-","synonymous_variant","missense_variant","stop_lost","stop_gained","start_lost","inframe_deletion","frameshift_variant","disruptive_inframe_deletion")
                            ),symb=c(3,4,8,20,0,5,6),ylims=c(0.05,1.0)){
  df=as.data.frame(rowRanges(vcf))[,c("seqnames","start")]
  region_idx = which(df$seqnames == chrom)
  df$DP=geno(vcf)$DP[,sample]
  df_region=df[region_idx,]
  df_region$xpos=seq(1,length(df_region[,1]))
  if (chrom_len) {
    xlimit=c(1,chrom_len)
    xcol = "start"
  }
  else {
    xlimit=c(df_region$xpos[1]-0.25,df_region$xpos[length(df_region$xpos)]+0.15)
    xcol = "xpos"
  }
  par(mar = c(0, 4, 2, 4))
  plot(df_region[,xcol], rep.int(0.5,length(df_region[,xcol])), pch=NULL, col="white", main=sample,log="y",
       ylab="AF(%)",xlim=xlimit, xlab=NA,xaxt='n', yaxt="n",
       ylim=c(0.005,1))
  axis(2,at=c(0.01,0.05,0.1,0.5,1.0),labels = FALSE )
  mtext(c("1","5","10","50","100"),side=2,at=c(0.01,0.05,0.1,0.5,1.0),line=1,outer=F,cex=0.60,las=1)
  abline(h = c(0.01,0.1,0.5,1.0), lty=3,col="lightgrey")
  afs=get_afs_annos(vcf,region_idx,sample)
  if ( length(afs) > 0 ) {
    afs$start = sapply(afs$xpos,function (x)  df_region$start[df_region$xpos == x ])
    points(afs[,xcol],afs$AF,col=cols[afs$ANNO],cex=1,pch=symb[afs$ALL]) 
  }
  par(new = T)
  plot(df_region[,xcol], df_region$DP, type="l", axes=F, log="y", xlab=NA, ylab=NA,col="darkgrey")
  axis(side = 4)
  mtext(side = 4, line = 2.75, 'Depth',cex=0.75)
  #par(mar = c(0, 4.1, 4.1, 2.1))
  return(df_region[,c("start","xpos")])
  
}


#for drawing genome features 
get.xpos.from.coords <- function(df,coords) {
  stretch <-  df$xpos[ coords[1] <= df$start & df$start <= coords[2]]
  return(c(min(stretch)-0.25,max(stretch)+0.25))
}

load_vcf_file <- function(vcf_file,scan_form=c("DP","AF"), scan_inf="ANN") {
  svp <- ScanVcfParam(info=scan_inf,geno=scan_form)
  vcf <- readVcf( vcf_file,"viruses",svp )
  annos=vcfInfo(vcf)[["ANN"]]       
  annos=lapply(annos, function (x) if( length(x) == 0 ) {return (NA)} else {return(unlist(strsplit(x,",",fixed=TRUE)))})
  idx = which(sapply(1:length(annos), function(x) ! is.na(annos[[x]][1]) ) )
  annos[idx]=sapply(idx, function (x) sapply(strsplit(annos[[x]],"|",fixed=TRUE), "[",2 ))  
  vcfInfo(vcf)[["ANN"]] <- annos
  return(vcf)
}

plot_element <- function (gR_entry, chrom_len, rh = 0.5,xoff=0.01,
                          cols=c("grey","darkgrey","lightblue","lightblue","lightgrey","darkgrey"),label=TRUE) {
  xstart=start(gR_entry)
  xstop=end(gR_entry)
  ym = 1
  xoff=xoff*chrom_len
  if (as.character(strand(gR_entry)) == "-") { ym = -1 }
  if ( gR_entry$type == "gene") {
    xoff=100
    xs = c(xstart,xstart,xstart,xstop,xstop,xstop)
    ys = c(0,ym*0.5,ym,ym,ym*0.5,0)*rh
    if (ym < 0) { xs[c(1,3)] = xs[c(1,3)] + xoff } else {xs[c(4,6)] = xs[c(4,6)] - xoff }
    polygon(xs, ys, col = cols[gR_entry$type])
    
  }
  else {
    rect(xstart, 0, xstop, ym*rh,col=cols[gR_entry$type])
  }
  if (label) { 
    text((xstart+xstop)/2,ym*0.5*rh,labels=c(gR_entry$gene),cex=0.65,srt=90)
  }
} 



##################
## Code for drawing
#################


basedir="~/LCMV_Data/"
gff_virus <- import.gff(paste0(basedir,"/References/viruses_short.gff"))

# read vcf file
vcf_file=paste0(basedir,"/Run_0355/lofreq2_all_samp_bed_norm_0.05_snpeff.vcf")
#vcf <- read.vcfR( vcf_file, verbose = FALSE )
svp <- ScanVcfParam(info="ANN",geno=c("DP","AF"))
vcf <- readVcf( vcf_file,"viruses",svp )
# only get the second field of each functional annotation for each allele and put it into the ANN field, NA if none
annos=vcfInfo(vcf)[["ANN"]]       
annos=lapply(annos, function (x) if( length(x) == 0 ) {return (NA)} else {return(unlist(strsplit(x,",",fixed=TRUE)))})
idx = which(sapply(1:length(annos), function(x) ! is.na(annos[[x]][1]) ) )
annos[idx]=sapply(idx, function (x) sapply(strsplit(annos[[x]],"|",fixed=TRUE), "[",2 ))  
vcfInfo(vcf)[["ANN"]] <- annos

vir_clean = gff_virus[ ! (gff_virus$type %in% c("CDS","databank_entry")) ]
mcols(vir_clean)$gene=c("5pUTR","GP","IGR","NP","3pUTR","5pUTR","Z","IGR","L","3pUTR")

samples <- rownames(colData(vcf))

layout(matrix(1:5, ncol = 1), widths = 1, respect = FALSE)
chrom_lens=c(3377,7229)
names(chrom_lens) = c("S","L")
chrom="S"
for (i in 1:4) {
  print(paste("plotting sample",samples[i]))
  reg <- plot_AFs_depths(vcf,chrom,samples[i],chrom_len = chrom_lens[chrom])
}


# idxst=seq(1,48,by = 5)
smpls=split(samples,ceiling(seq_along(samples)/5))
for (smps in smpls){
  for (chrom in c("S","L")){
    layout(matrix(1:(length(smps)+2), ncol = 1), widths = 1, respect = FALSE)
    for (i in smps) {
      print(paste("plotting sample",i))
      reg <- plot_AFs_depths(vcf,chrom,i,chrom_len = chrom_lens[chrom])
    }
    par(mar = c(3, 4, 0, 4))
    plot(1,10,pch=NULL,xlim=c(1,chrom_lens[chrom]),ylim=c(-1,1),ylab=NA,axes=FALSE,xlab=chrom,bty="n")
    axis(side=1,at=0:(chrom_lens[chrom]/1000)*1000)
    mtext(chrom, side=1, line=2)
    for( i in 1:length(vir_clean[seqnames(vir_clean) == chrom])){
      plot_element(vir_clean[seqnames(vir_clean) == chrom][i],chrom_lens[chrom],label=FALSE)
    }
    plot.new()
    legend("bottom",c("-","synonymous_variant","missense_variant","stop_lost","stop_gained","start_lost","inframe_deletion","frameshift_variant","disruptive_inframe_deletion","alt. allele 1", "alt. allele 2","alt. allele 3"),
           col=c("black","green","red","orangered","orange","orangered4","pink","purple","purple",rep("black",3)),
           pch=c(rep(3,9),3,4,8),cex=0.75,ncol=4)
    dev.copy2pdf(file=paste0("plot_S",smps[1],"_",smps[length(smps)],"_",chrom,".pdf"))
  }
}


q(save="no")

vcf_file="/Volumes/vetlinux01/LCMV/Run_0355/VarDict_2/all_samps_vardict_filt_norm_0.001_snpeff_snp_only.vcf"
scan_form=c("DP","AF","AD")
scan_inf="ANN"
svp <- ScanVcfParam(info=scan_inf,geno=scan_form)
vcf <- readVcf( vcf_file,"viruses",svp )
vcf2 <- expand(vcf,row.names=T)
geno(vcf)[["AD"]][4,c("S01")]
head(geno(vcf2)[["AD"]])
head(geno(vcf2)[["DP"]])
library(adegenet)

vcf2genpop <- function(vcf_file){
  # should create a genpop object from a vcf of pooled indv.
  # uses only AD field as of now, read depths becomes number of individuals
  # if NA/. and DP set use all reference as allele, else NA
  
  
}


