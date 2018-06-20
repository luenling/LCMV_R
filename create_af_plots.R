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
basedir="~/Data/LCMV_project/"
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

# get sample sheet from google drive

library(gsheet)
sample_sheet="https://docs.google.com/spreadsheets/d/1L2u3CZV2v75bsRhfepGThshCqUjNhA8bzo453bDohSQ/edit?usp=sharing"
smp_sh = gsheet2tbl(sample_sheet)
smp_sh <- as.data.frame(lapply(smp_sh, factor)) 
smp_sh = smp_sh[,1:4]
smp_sh$type = as.factor(gsub('\\.\\d+$',"",smp_sh$Description,perl=T))
smp_sh$txo = as.factor(paste(smp_sh$type,smp_sh$Origin,sep="x"))
smp_sh$txoxd = as.factor(paste(smp_sh$type,smp_sh$Origin,smp_sh$Days,sep="x"))


# load vcf file, should only contain SNPs and fully be annotated with AD2 and DPs
vcf_file="/Volumes/vetlinux01/LCMV/Run_0355/VarDict_2/all_samps_vardict_filt_norm_0.01_snpeff_snp_only.vcf"
#vcf_file="all_samps_vardict_filt_norm_0.01_snpeff_snp_only.vcf"
scan_form=c("DP","RD","AF","AD","AD2")
scan_inf="ANN"
svp <- ScanVcfParam(info=scan_inf,geno=scan_form)
vcf <- readVcf( vcf_file,"viruses",svp,row.names=T )
vcf2 <- expand(vcf,row.names=T)
geno(vcf)[["AD"]][4,c("S01")]
head(geno(vcf2)[["AD"]])
head(geno(vcf2)[["DP"]])
library(adegenet)

vcf2genpop <- function(vcf) {
  # should create a genpop object from a vcf of pooled indv.
  # uses preferentially AD field fallback to AFs, read depths becomes number of individuals
  # if NA/. and DP set use all reference as allele, else NA
  samps = rownames(colData(vcf))
  vcf_tab_t = matrix(NA, ncol = length(samps), nrow = 0)
  geno_names=names(geno(vcf))
  for (loc in rownames(vcf) ) {
    l_base <- gsub("_[ACTG/]+$","",loc,perl=TRUE)
    l_base <- gsub(":","_",l_base,fixed=TRUE)
    a.num <- length(fixed(vcf[loc,])$ALT[[1]]) + 1
    if ("AD" %in% geno_names) {
      # vardict or other caller using AD field 
      rvecs <- vapply(geno(vcf)[["AD"]][loc,],'[',1:a.num,1:a.num)
      # check if reference depth is set use this to fill depths of some NAs
      if ("RD" %in% geno_names) {
        rvecs[1,is.na(rvecs[1,])] = rowSums(geno(vcf)$RD[loc,is.na(rvecs[1,]),])
      }
      # use the depths calculated by samtools mpileup (in AD2) or finally set to DP
      rvecs[1,is.na(rvecs[1,])] = if ( length(geno(vcf)$AD2[loc,1][[1]]) > 0 ) vapply(geno(vcf)$AD2[loc,is.na(rvecs[1,])],'[[',1,1) else 100
      # not reference allele fill NA with 0
      rvecs[2:a.num,][which(is.na(rvecs[2:a.num,]))] = 0
    } else if ("AF" %in% geno_names) {
      # use allele frequencies and DP to get approximate read counts
      rvecs <- vapply(geno(vcf)[["AF"]][loc,],'[',as.double(1:(a.num-1)) )
      rvecs[is.na(rvecs)] <- 0
      
      rvecs <- if(a.num < 3)  rbind(1-rvecs,rvecs)  else rbind(1-colSums(rvecs),rvecs)
      rvecs <- t(apply(rvecs,1,function(x) round(x * geno(vcf)$DP[loc,])))
    }
    # add rows to transposed tab and name them
    vcf_tab_t = rbind(vcf_tab_t,rvecs)
    rownames(vcf_tab_t)[(nrow(vcf_tab_t)-a.num+1):nrow(vcf_tab_t)] = paste0(l_base,".",1:a.num )
  }
  return(t(vcf_tab_t))
}

vcf_tab<-vcf2genpop(vcf)
vcf_pop <- as.genpop(vcf_tab,ploidy=1,type="codom")
d_nei <- dist.genpop(vcf_pop)
ca1 <- dudi.coa(tab(vcf_pop),scannf=FALSE,nf=3)
barplot(ca1$eig,main="Correspondance Analysis eigenvalues",
        col=heat.colors(length(ca1$eig)))
s.label(ca1$li, sub="CA 1-2",csub=2)
add.scatter.eig(ca1$eig,nf=3,xax=1,yax=2,posi="bottomright")
s.label(ca1$li, sub="CA 1-2",csub=2)
add.scatter.eig(ca1$eig,nf=3,xax=1,yax=2,posi="bottomright")

s.label(ca1$li,xax=2,yax=3,lab=popNames(vcf_pop),sub="CA 1-3",csub=2)
add.scatter.eig(ca1$eig,nf=3,xax=2,yax=3,posi="topleft")

s.class(ca1$li,fac=smp_sh$txoxd,xax=2,yax=3,label=NULL,
        col=fac2col(levels(smp_sh$txoxd)),sub="CA 1-2",csub=1)

b6rag2 = factor(smp_sh$Sample[smp_sh$type ==  "B6-RAG2-/-LY5" & smp_sh$Sample != "S15"] )
sm_b6rag2 <- subset(smp_sh, Sample %in% b6rag2)
sm_b6rag2 <- lapply(sm_b6rag2, function(x) if(is.factor(x)) factor(x) else x)
vcf_pop_b6rag2 <- vcf_pop[b6rag2,]

library(corrplot)
corrplot(as.matrix(dist.genpop(vcf_pop_b6rag2,method=2)), type = "full",is.corr=F)

ca1 <- dudi.coa(tab(vcf_pop_b6rag2),scannf=F,nf=4)
ca2 <- dudi.pco(dist.genpop(vcf_pop_b6rag2,method=2), nf=4,scannf=F)
plot(ca2$li[,1],ca2$li[,2],pch=19,col=fac2col(sm_b6rag2$Origin),xlab="PC1",ylab="PC2")
abline(h=0)
abline(v=0)
legend("bottomright",legend=levels(sm_b6rag2$Origin),fill=fac2col(levels(sm_b6rag2$Origin)),ncol=2)
text(ca2$li[,1],ca2$li[,2],labels = rownames(ca2$li), cex=0.75)
add.scatter.eig(ca2$eig,nf=4,xax=1,yax=2,posi="bottomright",ratio=0.2)

# find clusters
find.clusters(tab(vcf_pop_b6rag2),max.n=6)

# hclust tree
hc_b6rag2 = hclust(dist.genpop(vcf_pop_b6rag2,method=1), method= "complete" )
plot(hc_b6rag2,labels=sm_b6rag2$Sample,col=fac2col(sm_b6rag2$Origin))

# MDS
mds_b6rag2 <- cmdscale(dist.genpop(vcf_pop_b6rag2,method=2),eig=TRUE,k=4,x.ret=TRUE)
plot(mds_b6rag2$points[,1],mds_b6rag2$points[,2])


s.class(ca1$li,fac=sm_b6rag2$Origin,xax=1,yax=2,label=NULL,col=fac2col(levels(sm_b6rag2$Origin)),sub="CA 1-2",csub=1)
add.scatter.eig(ca1$eig,nf=3,xax=1,yax=2,posi="bottomright",ratio=0.2)
plot(ca1$li[,1],ca1$li[,2],pch=19,col=fac2col(sm_b6rag2$Origin),xlab="PC1",ylab="PC2")
abline(h=0)
abline(v=0)
legend("bottomright",legend=levels(sm_b6rag2$Origin),fill=fac2col(levels(sm_b6rag2$Origin)),ncol=2)
add.scatter.eig(ca1$eig,nf=3,xax=1,yax=2,posi="topright",ratio=0.15)

s.class(ca1$li,fac=sm_b6rag2$Origin,xax=2,yax=3,label=NULL,col=fac2col(levels(sm_b6rag2$Origin),col.pal=funky),sub="CA 2-3",csub=1)
add.scatter.eig(ca1$eig,nf=3,xax=2,yax=3,posi="bottomright",ratio=0.2)


#lofreq
lofreq5_vcf="/Volumes/vetlinux01/LCMV/Run_0355/BQSR/lofreq2_all_samp_bed_norm_0.05_snpeff_snp_only.vcf"
lofreq_vcf="/Volumes/vetlinux01/LCMV/Run_0355/BQSR/lofreq2_all_samp_bed_norm_0.01_snpeff_snp_only.vcf"
scan_form=c("DP","AF")
scan_inf="ANN"
svp <- ScanVcfParam(info=scan_inf,geno=scan_form)
lof_vcf <- readVcf( lofreq_vcf,"viruses",svp,row.names=T )
lof5_vcf <- readVcf( lofreq5_vcf,"viruses",svp,row.names=T )

lof_tab<-vcf2genpop(lof_vcf)
lof_pop <- as.genpop(lof_tab,ploidy=1,type="codom")
lof5_tab<-vcf2genpop(lof5_vcf)
lof5_pop <- as.genpop(lof5_tab,ploidy=1,type="codom")
lof5_d_nei <- dist.genpop(lof5_pop)


b6rag2 = factor(smp_sh$Sample[smp_sh$type ==  "B6-RAG2-/-LY5"] )
sm_b6rag2 <- subset(smp_sh, Sample %in% b6rag2)
sm_b6rag2 <- as.data.frame(lapply(sm_b6rag2, function(x) if(is.factor(x)) factor(x) else x))
rownames(sm_b6rag2) <- sm_b6rag2$Sample
lof_pop_b6rag2 <- lof_pop[b6rag2,]
lof5_pop_b6rag2 <- lof5_pop[b6rag2,]


# 5% MAF Lofreq
ca1 <- dudi.coa(tab(lof5_pop_b6rag2),scannf=FALSE,nf=3)
barplot(ca1$eig,main="Correspondance Analysis eigenvalues",
        col=heat.colors(length(ca1$eig)))
s.class(ca1$li,fac=sm_b6rag2$txoxd,xax=1,yax=2,label=NULL,
        col=fac2col(levels(sm_b6rag2$txoxd)),sub="CA 1-2",csub=1)
legend("topright",legend=levels(sm_b6rag2$txoxd),fill=fac2col(levels(sm_b6rag2$txoxd)),ncol=2,cex=0.5)

plot(ca1$li[,1],ca1$li[,2],pch=19,col=fac2col(sm_b6rag2$Origin),xlab="PC1",ylab="PC2",main=paste("Cor. Anal, MAF5%:",sm_b6rag2$type[1],sm_b6rag2$Days[1],sep="  ") )
#plot(ca1$li[,1],ca1$li[,3],pch=19,col=fac2col(sm_b6rag2$Origin),xlab="PC1",ylab="PC3")
abline(h=0)
abline(v=0)
legend("topright",legend=levels(sm_b6rag2$Origin),fill=fac2col(levels(sm_b6rag2$Origin)),ncol=2)
text(ca1$li[,1],ca1$li[,2],labels = rownames(ca1$li), cex=0.75)
# looking at contributions
inertia.ca1 <- inertia.dudi(ca1, row.inertia = TRUE,col.inertia = TRUE)
ca1.pc1.maxall = rownames(inertia.ca1$col.abs)[order(inertia.ca1$col.abs[,1],decreasing = T)][1:10]
colSums(tab(lof5_pop_b6rag2)[,ca1.pc1.maxall] > 0.0)
# makes no sense - different depths on different sites
inertia.ca1$col.abs/100
ca1$co[rowSums(ca1$co)>0.1,]


b6rag2 = factor(smp_sh$Sample[smp_sh$type ==  "B6-RAG2-/-LY5" & ! smp_sh$Sample %in% c("S15")] )
sm_b6rag2 <- subset(smp_sh, Sample %in% b6rag2)
sm_b6rag2 <- as.data.frame(lapply(sm_b6rag2, function(x) if(is.factor(x)) factor(x) else x))
rownames(sm_b6rag2) <- sm_b6rag2$Sample
lof5_pop_b6rag2 <- lof5_pop[b6rag2,]

ca2 <- dudi.pco(dist.genpop(lof5_pop_b6rag2,method=2), nf=4,scannf=F)
plot(ca2$li[,1],ca2$li[,2],pch=19,col=fac2col(sm_b6rag2$Origin),xlab="PC1",ylab="PC2",
     main=paste("PCoA MAF5% Reynolds distance:",sm_b6rag2$type[1],sm_b6rag2$Days[1],sep=" "))
#plot(ca2$li[,1],ca2$li[,3],pch=19,col=fac2col(sm_b6rag2$Origin),xlab="PC1",ylab="PC3")
abline(h=0)
abline(v=0)
legend("topright",legend=levels(sm_b6rag2$Origin),fill=fac2col(levels(sm_b6rag2$Origin)),ncol=2)
text(ca2$li[,1],ca2$li[,2],labels = rownames(ca2$li), cex=0.75)



b6rag2 = factor(smp_sh$Sample[smp_sh$type ==  "B6-RAG2-/-LY5" & ! smp_sh$Sample %in% c("S15","S13")] )
sm_b6rag2 <- subset(smp_sh, Sample %in% b6rag2)
sm_b6rag2 <- as.data.frame(lapply(sm_b6rag2, function(x) if(is.factor(x)) factor(x) else x))
rownames(sm_b6rag2) <- sm_b6rag2$Sample
lof5_pop_b6rag2 <- lof5_pop[b6rag2,]

tabs = tab(lof5_pop_b6rag2,freq=T)
tabs = tabs[, - grep('\\.1$',colnames(tabs),perl=T)]
ca3 <- dudi.pca(tabs,scale = FALSE, scannf = FALSE, nf = 3)
barplot(ca3$eig)
plot(ca3$li[,1],ca3$li[,2],pch=19,col=fac2col(sm_b6rag2$Origin),xlab="PC1",ylab="PC2",
     main=paste("PCA MAF5%:",sm_b6rag2$type[1],sm_b6rag2$Days[1],sep=" "))
#plot(ca3$li[,2],ca3$li[,3],pch=19,col=fac2col(sm_b6rag2$Origin),xlab="PC2",ylab="PC3")
abline(h=0)
abline(v=0)
legend("topright",legend=levels(sm_b6rag2$Origin),fill=fac2col(levels(sm_b6rag2$Origin)),ncol=2)
text(ca3$li[,1],ca3$li[,2],labels = rownames(ca3$li), cex=0.75)
# looking at contributions
inertia.ca3 <- inertia.dudi(ca3, row.inertia = TRUE,col.inertia = TRUE)
ca3.pc1.maxall = rownames(inertia.ca3$col.abs)[order(inertia.ca3$col.abs[,1],decreasing = T)][1:9]
colSums(tab(lof5_pop_b6rag2,freq = T)[,ca3.pc1.maxall] > 0.0)
round(tab(lof5_pop_b6rag2,freq = T)[,ca3.pc1.maxall],digits = 3)
summary(inertia.ca3)
# only 4 loci really contribute to first axis
inertia.ca1$col.abs/100
ca1$co[rowSums(ca1$co)>0.1,]

# drawing some trees
library(ape)
library(dendextend)
hc1 = as.dendrogram(hclust(dist.genpop(lof5_pop_b6rag2,method=4),method="complete"))
hc2 = as.dendrogram(hclust(dist.genpop(lof5_pop_b6rag2,method=4),method="single"))
hc3 = as.dendrogram(hclust(dist.genpop(lof5_pop_b6rag2,method=4),method="ward.D"))
lof5_pop_b6rag2 %>% dist.genpop(method=4) %>% hclust(method="complete") %>% 
  as.dendrogram() %>% set("labels_col", fac2col(sm_b6rag2[labels(.),"Origin"])) %>% 
  plot()


tanglegram(hc1,hc3)
#plot(hc, col=fac2col(sm_b6rag2$Origin))
#labels(hc)
hc1 %>%  set("labels_col", fac2col(sm_b6rag2[labels(hc),"Origin"])) %>% plot()
legend("topright",legend=levels(sm_b6rag2$Origin),fill=fac2col(levels(sm_b6rag2$Origin)),ncol=2)


nj_lof5_b6rag2 = bionj(dist.genpop(lof5_pop_b6rag2,method=2))
plot(nj_lof5_b6rag2,show.tip.label=F)
tiplabels(nj_lof5_b6rag2$tip.label, adj = c(0.0, 0.5), frame = "none" , bg = NULL, col = fac2col(sm_b6rag2$Origin))
legend("topright",legend=levels(sm_b6rag2$Origin),fill=fac2col(levels(sm_b6rag2$Origin)),ncol=2)


b6rag2 = factor(smp_sh$Sample[smp_sh$type ==  "B6-RAG2-/-LY5" & smp_sh$Sample != "S15"] )
sm_b6rag2 <- subset(smp_sh, Sample %in% b6rag2)
sm_b6rag2 <- lapply(sm_b6rag2, function(x) if(is.factor(x)) factor(x) else x)
lof_pop_b6rag2 <- lof_pop[b6rag2,]

#library(corrplot)
corrplot(as.matrix(dist.genpop(lof_pop_b6rag2,method=4)), type = "full",is.corr=F)

ca1 <- dudi.coa(tab(lof_pop_b6rag2),scannf=F,nf=4)
plot(ca1$li[,1],ca1$li[,2],pch=19,col=fac2col(sm_b6rag2$Origin),xlab="PC1",ylab="PC2",main=paste("Cor. Anal:",sm_b6rag2$type[1],sm_b6rag2$Days[1],sep="  ") )
#plot(ca1$li[,1],ca1$li[,3],pch=19,col=fac2col(sm_b6rag2$Origin),xlab="PC1",ylab="PC3")
abline(h=0)
abline(v=0)
legend("topleft",legend=levels(sm_b6rag2$Origin),fill=fac2col(levels(sm_b6rag2$Origin)),ncol=2)
text(ca1$li[,1],ca1$li[,2],labels = rownames(ca1$li), cex=0.75)


ca2 <- dudi.pco(dist.genpop(lof_pop_b6rag2,method=2), nf=4,scannf=F)
plot(ca2$li[,1],ca2$li[,2],pch=19,col=fac2col(sm_b6rag2$Origin),xlab="PC1",ylab="PC2")
plot(ca2$li[,1],ca2$li[,3],pch=19,col=fac2col(sm_b6rag2$Origin),xlab="PC1",ylab="PC3")
abline(h=0)
abline(v=0)
legend("bottomright",legend=levels(sm_b6rag2$Origin),fill=fac2col(levels(sm_b6rag2$Origin)),ncol=2)
text(ca2$li[,1],ca2$li[,2],labels = rownames(ca2$li), cex=0.75)

add.scatter.eig(ca2$eig,nf=4,xax=1,yax=2,posi="bottomright",ratio=0.2)

c57bl_6 = factor(smp_sh$Sample[smp_sh$type ==  "C57BL/6" & smp_sh$Days =="40 dpi" ] )
sm_c57bl_6 <- subset(smp_sh, Sample %in% c57bl_6)
sm_c57bl_6 <- lapply(sm_c57bl_6, function(x) if(is.factor(x)) factor(x) else x)
lof_pop_c57bl_6 <- lof_pop[c57bl_6,]

corrplot(as.matrix(dist.genpop(lof_pop_c57bl_6,method=2)), type = "full",is.corr=F)

ca1 <- dudi.coa(tab(lof_pop_c57bl_6),scannf=F,nf=4)
ca2 <- dudi.pco(dist.genpop(lof_pop_c57bl_6,method=2), nf=4,scannf=F)
plot(ca2$li[,1],ca2$li[,2],pch=19,col=fac2col(sm_c57bl_6$Origin),xlab="PC1",ylab="PC2")
plot(ca2$li[,2],ca2$li[,3],pch=19,col=fac2col(sm_c57bl_6$Origin),xlab="PC2",ylab="PC3")



#library(corrplot)
corrplot(as.matrix(dist.genpop(lof_pop_b6rag2,method=2)), type = "full",is.corr=F)



pairwise.fst(vcf_pop)
dudi.pco(d_nei)
library(adegenet)
library(ade4)
