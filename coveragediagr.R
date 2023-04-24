setwd("/Volumes/Temp/LCMV/Data/")
#setwd("/Volumes//vetgrid01/LCMV_project//Data/BSF_0176_H325VBBXX_5_samples/Coverages/")
setwd("/Volumes//vetgrid01/LCMV_project/Run_0277/Mappings/Coverages")
virus.covs=read.table("bsf_0277_all.coverages",header=T)
#levels(virus.covs$CHROM)=c("S","L")
primers=read.table("/Volumes//vetgrid01/LCMV_project/References/primers.bed")
colnames(primers)=c("CHROM","Start","End","Primer","MQ","Strand","Flag","CIGAR","MateCHR")
# bed format 0 based open intervals:
primers$Start = primers$Start + 1
levels(primers$CHROM)
levels(primers$CHROM)=c("S","L")

l.covs=virus.covs[virus.covs$CHROM=="L",2:ncol(virus.covs)]
s.covs=virus.covs[virus.covs$CHROM=="S",2:ncol(virus.covs)]
l.covs$windows=cut(l.covs$BPS,seq(0,max(l.covs$BPS),by=50))
s.covs$windows=cut(s.covs$BPS,seq(0,max(s.covs$BPS),by=50))

l.covs.50bp=aggregate(l.covs[,1:(ncol(l.covs)-1)],list(l.covs$win),mean)
s.covs.50bp=aggregate(s.covs[,1:(ncol(s.covs)-1)],list(s.covs$win),mean)
matplot(l.covs.50bp$BPS,l.covs.50bp[,3:5],type="l",xlab="segment L",ylab="coverage (50bp average)")
matplot(l.covs.50bp$BPS,l.covs.50bp[,3:ncol(l.covs.50bp)],type="l",xlab="segment L",ylab="coverage (50bp average)")

pdf("l.cov.pdf",width=8,height=6)
matplot(l.covs.50bp$BPS,l.covs.50bp[,4:ncol(l.covs.50bp)]+0.01,type="l",xlab="base position",ylab="coverage (50bp average)", log="y", main="L segment")
dev.off()
pdf("s.cov.pdf",width=8,height=6)
matplot(s.covs.50bp$BPS,s.covs.50bp[,4:ncol(s.covs.50bp)]+0.01,type="l",xlab="base position",ylab="coverage (50bp average)", log="y", main="S segment")
dev.off()

colSums(l.covs[,2:(ncol(l.covs)-1)] > 100)/nrow(l.covs)
colSums(s.covs[,2:(ncol(s.covs)-1)] > 100)/nrow(s.covs)

library(ggplot2)
library(RColorBrewer)
library("reshape2")
l.covs.50bp.melted <- melt(l.covs.50bp[,c(2,4:ncol(l.covs.50bp))],id.vars="BPS",variable.name = "sample", value.name="coverage")
s.covs.50bp.melted <- melt(s.covs.50bp[,c(2,4:ncol(s.covs.50bp))],id.vars="BPS",variable.name = "sample", value.name="coverage")

bp = ggplot(data=l.covs.50bp.melted, aes(x=BPS, y=coverage, group = sample, linetype=sample, colour = sample)) +
  geom_line() + scale_y_log10() + scale_linetype_manual(values = c(rep("solid", 12), rep("dashed", 11))) +
  scale_color_manual(values = c(brewer.pal(9, "Set1"), brewer.pal(6,"Accent"),brewer.pal(8, "Dark2"))) + theme_bw() + 
  theme(legend.direction ="horizontal",legend.position = "bottom") + guides(col=guide_legend(nrow=3)) + ggtitle("L Segment") 

primers$Y=10
primers$YE=10

bp + geom_vline(xintercept = primers$Start[primers$CHR == "L" ],colour="black", linetype = "dotted", linewidth=5) +
 annotate("text", x =primers$Start[primers$CHR == "L"], y=9*10^6, label=gsub("RT_","",primers$Primer[primers$CHR == "L"]), size=3)
# geom_text(data= primers[primers$CHR == "L",], aes(x = Start, y = Start, label= Primer))

bpS = ggplot(data=s.covs.50bp.melted, aes(x=BPS, y=coverage, group = sample, linetype=sample,colour = sample)) +
  geom_line() + scale_y_log10() + scale_linetype_manual(values = c(rep("solid", 9), rep("dashed", 7))) +
  scale_color_manual(values = c(brewer.pal(9, "Set1"), brewer.pal(8, "Dark2"))) + theme_bw() +
  guides(linetype=guide_legend(ncol=2)) + ggtitle("S Segment") + ggsave("ssegment.pdf",width=8,height=6)
bpS + ggsave("ssegment.pdf",width=8,height=6)
bpS +  geom_vline(xintercept = primers$Start[primers$CHR == "S" ],colour="black", linetype = "dotted", linewidth=5) +
  annotate("text", x =primers$Start[primers$CHR == "S"], y=primers$Start[primers$CHR == "S"]*0+c(9*10^6,7*10^6), label=gsub("RT_","",primers$Primer[primers$CHR == "S"]), size=3) + ggsave("ssegment_primers.pdf",width=8,height=6)


bpL = ggplot(data=l.covs.50bp.melted, aes(x=BPS, y=coverage, group = sample, linetype=sample, colour = sample)) +
  geom_line() + scale_y_log10() + scale_linetype_manual(values = c(rep("solid", 9), rep("dashed", 7))) +
  scale_color_manual(values = c(brewer.pal(9, "Set1"), brewer.pal(8, "Dark2"))) + theme_bw() + 
  guides(linetype=guide_legend(ncol =2)) + ggtitle("L Segment") 
bpL + ggsave("lsegment.pdf",width=8,height=6)
bpL +  geom_vline(xintercept = primers$Start[primers$CHR == "L" ],colour="black", linetype = "dotted") +
  annotate("text", x =primers$Start[primers$CHR == "L"], y=primers$Start[primers$CHR == "L"]*0+c(9*10^6,7*10^6), label=gsub("RT_","",primers$Primer[primers$CHR == "L"]), size=3) + ggsave("Lsegment_primers.pdf",width=8,height=6)

low_samps =c("S_1","S_2","S_3","S_4","S_5","S_8","S_9","S_10")

bpS = ggplot(data=s.covs.50bp.melted[s.covs.50bp.melted$sample %in% low_samps,], aes(x=BPS, y=coverage, group = sample, linetype=sample,colour = sample)) +
  geom_line() +  scale_linetype_manual(values = c(rep("solid", 9), rep("dashed", 7))) +
  scale_color_manual(values = c(brewer.pal(9, "Set1"), brewer.pal(8, "Dark2"))) + theme_bw() +
  guides(linetype=guide_legend(ncol=2)) + ggtitle("S Segment") + ggsave("ssegment_low_covs.pdf",width=8,height=6)
bpS + ggsave("ssegment_low_covs.pdf",width=8,height=6)
bpS +  geom_vline(xintercept = primers$Start[primers$CHR == "S" ],colour="black", linetype = "dotted", linewidth=5) +
  annotate("text", x =primers$Start[primers$CHR == "S"], y=primers$Start[primers$CHR == "S"]*0+c(9*10^6,7*10^6), label=gsub("RT_","",primers$Primer[primers$CHR == "S"]), size=3) + ggsave("ssegment_primers_low.pdf",width=8,height=6)








setwd("/Volumes/vetgrid01/LCMV_project/Variants/Varseq2/")
dat = dcast(m2data, type + index ~ variable, value.var="position")





variants=read.table("ds1000_varscan_perc_head.tab",header=T,na.strings = c("."))
dos = read.table("ds1000_varscan_DP_head.tab",header=T,na.strings = c("."))

all.var=variants[ ! is.na(rowSums(variants[,5:27])),]
fit <- princomp(sqrt(all.var[,5:27]/100), cor=TRUE)
summary(fit) # print variance accounted for 
loadings(fit) # pc loadings 
plot(fit,type="lines") # scree plot 
fit$scores # the principal components
biplot(fit)
plot(fit)
library(psych)
fit <- factor.pa(all.var[,5:27], nfactors=2, rotation="varimax")
fit # print results
str(fit)
library("FactoMineR", lib.loc="/Library/Frameworks/R.framework/Versions/3.1/Resources/library")
result <- PCA(all.var[,5:27])
library("PopGenome")
GENOME.class <- readData("ds1000_varscan.vcf", format="VCF", include.unknown=TRUE)

fsts.varscan=read.table("../ds1000_varscan_snps.fst.tab",header=TRUE)
fsts.median=sapply(fsts.varscan[,3:ncol(fsts.varscan)], median)
fsts.mean=sapply(fsts.varscan[,3:ncol(fsts.varscan)], median)

dm.fst.mean=matrix(0,nrow=23,ncol=23)
dimnames(dm.fst.mean)=list(colnames(all.var[,5:27]),colnames(all.var[,5:27]))
a=c()
for(i in 1:22){
  for(j in (i+1):23){
    a=append(a,( (22 - i / 2 ) * (i-1) + j -1 ) )
    dm.fst.mean[j,i] =  fsts.mean[ ( (22 - i / 2 ) * (i-1) + j -1 ) ]
  }
} 
dm.fst.mean=as.dist(dm.fst.mean,diag = FALSE, upper = FALSE)
hc <- hclust(dm.fst.mean)
plot(hc)

setwd("/Volumes/vetgrid01/LCMV_project/NewData/bwa_mem_mm10/Coverages/")
virus.covs=read.table("ND_bwa_mem_mm10_all.coverages",header=T)
levels(virus.covs$CHROM)=c("S","L")
primers=read.table("/Volumes//vetgrid01/LCMV_project/References/primers.bed")
colnames(primers)=c("CHROM","Start","End","Primer","MQ","Strand","Flag","CIGAR","MateCHR")
primers$Start = primers$Start + 1
levels(primers$CHROM)
levels(primers$CHROM)=c("S","L")
l.covs=virus.covs[virus.covs$CHROM=="L",2:ncol(virus.covs)]
s.covs=virus.covs[virus.covs$CHROM=="S",2:ncol(virus.covs)]
l.covs$windows=cut(l.covs$BPS,seq(0,max(l.covs$BPS),by=50))
s.covs$windows=cut(s.covs$BPS,seq(0,max(s.covs$BPS),by=50))

l.covs.50bp=aggregate(l.covs[,1:(ncol(l.covs)-1)],list(l.covs$win),mean)
s.covs.50bp=aggregate(s.covs[,1:(ncol(s.covs)-1)],list(s.covs$win),mean)

l.covs$windows=cut(l.covs$BPS,seq(0,max(l.covs$BPS),by=10))
s.covs$windows=cut(s.covs$BPS,seq(0,max(s.covs$BPS),by=10))

l.covs.10bp=aggregate(l.covs[,1:(ncol(l.covs)-1)],list(l.covs$win),mean)
s.covs.10bp=aggregate(s.covs[,1:(ncol(s.covs)-1)],list(s.covs$win),mean)

matplot(l.covs.50bp$BPS,l.covs.50bp[,3:5],type="l",xlab="segment L",ylab="coverage (50bp average)")
matplot(l.covs.50bp$BPS,l.covs.50bp[,3:ncol(l.covs.50bp)],type="l",xlab="segment L",ylab="coverage (50bp average)")


library(ggplot2)
library(RColorBrewer)
library("reshape2")
library("Gviz")
library("ggbio")
library("GenomicRanges")
library(GenomicAlignments)

l.covs.50bp.melted <- melt(l.covs.50bp[,c(2,3:ncol(l.covs.50bp))],id.vars="BPS",variable.name = "sample", value.name="coverage")
s.covs.50bp.melted <- melt(s.covs.50bp[,c(2,3:ncol(s.covs.50bp))],id.vars="BPS",variable.name = "sample", value.name="coverage")
l.covs.10bp.melted <- melt(l.covs.10bp[,c(2,3:ncol(l.covs.10bp))],id.vars="BPS",variable.name = "sample", value.name="coverage")
s.covs.10bp.melted <- melt(s.covs.10bp[,c(2,3:ncol(s.covs.10bp))],id.vars="BPS",variable.name = "sample", value.name="coverage")
colors = colorRampPalette(c("black", "yellow", "blue", "green", "red"))( 14 ) 
colors=c(rev(brewer.pal(6,"Set1")),"darkgrey")

cl13p=grep("Cl13.*plasmid",levels(l.covs.50bp.melted$sample))
cl13bh=grep("Cl13.*BHK",levels(l.covs.50bp.melted$sample))
brain=grep("Brain",levels(l.covs.50bp.melted$sample))
first=seq(1,6)
serum=grep("Serum",levels(l.covs.50bp.melted$sample))

displ=brain
bp = ggplot(data=l.covs.10bp.melted[l.covs.10bp.melted$sample %in% levels(l.covs.10bp.melted$sample)[displ],], aes(x=BPS, y=coverage, group = sample, linetype=sample, colour = sample)) +
  geom_line() + scale_y_log10() + scale_linetype_manual(values = c(rep("solid", 14), rep("dashed", 14))) +
  scale_color_manual(values = c(colors,colors)) + theme_bw() + 
  theme(legend.direction ="horizontal",legend.position = "bottom") + guides(col=guide_legend(nrow=3)) + ggtitle("L Segment") 
bp + geom_vline(xintercept = primers$Start[primers$CHR == "L" ], aes(colour="black", linetype = "dotted", linewidth=5)) +
  annotate("text", x =primers$Start[primers$CHR == "L"], y=9*10^6, label=gsub("RT_","",primers$Primer[primers$CHR == "L"]), size=3) + ggsave("lsegment_brains.pdf",width=8,height=6)

bps = ggplot(data=s.covs.10bp.melted[s.covs.10bp.melted$sample %in% levels(s.covs.10bp.melted$sample)[displ],], aes(x=BPS, y=coverage, group = sample, linetype=sample, colour = sample)) +
  geom_line() + scale_y_log10() + scale_linetype_manual(values = c(rep("solid", 14), rep("dashed", 14))) +
  scale_color_manual(values = c(colors,colors)) + theme_bw() + 
  theme(legend.direction ="horizontal",legend.position = "bottom") + guides(col=guide_legend(nrow=3)) + ggtitle("S Segment") 
bps + geom_vline(xintercept = primers$Start[primers$CHR == "S" ], aes(colour="black", linetype = "dotted", linewidth=5)) +
  annotate("text", x =primers$Start[primers$CHR == "S"], y=9*10^6, label=gsub("RT_","",primers$Primer[primers$CHR == "S"]), size=3) + ggsave("ssegment_brains.pdf",width=8,height=6)

displ=cl13p
bp = ggplot(data=l.covs.10bp.melted[l.covs.10bp.melted$sample %in% levels(l.covs.10bp.melted$sample)[displ],], aes(x=BPS, y=coverage, group = sample, linetype=sample, colour = sample)) +
  geom_line() + scale_y_log10() + scale_linetype_manual(values = c(rep("solid", 7), rep("dashed", 7))) +
  scale_color_manual(values = c(colors,colors)) + theme_bw() + 
  theme(legend.direction ="horizontal",legend.position = "bottom") + guides(col=guide_legend(nrow=3)) + ggtitle("L Segment") 
bp + geom_vline(xintercept = primers$Start[primers$CHR == "L" ], aes(colour="black", linetype = "dotted", linewidth=5)) +
  annotate("text", x =primers$Start[primers$CHR == "L"], y=9*10^6, label=gsub("RT_","",primers$Primer[primers$CHR == "L"]), size=3) + ggsave("lsegment_cl13_plasmid.pdf",width=8,height=6)

bps = ggplot(data=s.covs.10bp.melted[s.covs.10bp.melted$sample %in% levels(s.covs.10bp.melted$sample)[displ],], aes(x=BPS, y=coverage, group = sample, linetype=sample, colour = sample)) +
  geom_line() + scale_y_log10() + scale_linetype_manual(values = c(rep("solid", 7), rep("dashed", 7))) +
  scale_color_manual(values = c(colors,colors)) + theme_bw() + 
  theme(legend.direction ="horizontal",legend.position = "bottom") + guides(col=guide_legend(nrow=3)) + ggtitle("S Segment") 
bps + geom_vline(xintercept = primers$Start[primers$CHR == "S" ], aes(colour="black", linetype = "dotted", linewidth=5)) +
  annotate("text", x =primers$Start[primers$CHR == "S"], y=9*10^6, label=gsub("RT_","",primers$Primer[primers$CHR == "S"]), size=3) + ggsave("ssegment_cl13_plasmid.pdf",width=8,height=6)


displ=cl13bh
bp = ggplot(data=l.covs.10bp.melted[l.covs.10bp.melted$sample %in% levels(l.covs.10bp.melted$sample)[displ],], aes(x=BPS, y=coverage, group = sample, linetype=sample, colour = sample)) +
  geom_line() + scale_y_log10() + scale_linetype_manual(values = c(rep("solid", 7), rep("dashed", 7))) +
  scale_color_manual(values = c(colors,colors)) + theme_bw() + 
  theme(legend.direction ="horizontal",legend.position = "bottom") + guides(col=guide_legend(nrow=3)) + ggtitle("L Segment") 
bp + geom_vline(xintercept = primers$Start[primers$CHR == "L" ], aes(colour="black", linetype = "dotted", linewidth=5)) +
  annotate("text", x =primers$Start[primers$CHR == "L"], y=9*10^6, label=gsub("RT_","",primers$Primer[primers$CHR == "L"]), size=3) + ggsave("lsegment_cl13_bhk21.pdf",width=8,height=6)

bps = ggplot(data=s.covs.10bp.melted[s.covs.10bp.melted$sample %in% levels(s.covs.10bp.melted$sample)[displ],], aes(x=BPS, y=coverage, group = sample, linetype=sample, colour = sample)) +
  geom_line() + scale_y_log10() + scale_linetype_manual(values = c(rep("solid", 7), rep("dashed", 7))) +
  scale_color_manual(values = c(colors,colors)) + theme_bw() + 
  theme(legend.direction ="horizontal",legend.position = "bottom") + guides(col=guide_legend(nrow=3)) + ggtitle("S Segment") 
bps + geom_vline(xintercept = primers$Start[primers$CHR == "S" ], aes(colour="black", linetype = "dotted", linewidth=5)) +
  annotate("text", x =primers$Start[primers$CHR == "S"], y=9*10^6, label=gsub("RT_","",primers$Primer[primers$CHR == "S"]), size=3) + ggsave("ssegment_cl13_bhk21.pdf",width=8,height=6)

displ=first

bp = ggplot(data=l.covs.10bp.melted[l.covs.10bp.melted$sample %in% levels(l.covs.10bp.melted$sample)[displ],], aes(x=BPS, y=coverage, group = sample, linetype=sample, colour = sample)) +
  geom_line() + scale_y_log10() + scale_linetype_manual(values = c(rep("solid", 7), rep("dashed", 7))) +
  scale_color_manual(values = c(colors,colors)) + theme_bw() + 
  theme(legend.direction ="horizontal",legend.position = "bottom") + guides(col=guide_legend(nrow=3)) + ggtitle("L Segment") 
bp + geom_vline(xintercept = primers$Start[primers$CHR == "L" ], aes(colour="black", linetype = "dotted", linewidth=5)) +
  annotate("text", x =primers$Start[primers$CHR == "L"], y=9*10^6, label=gsub("RT_","",primers$Primer[primers$CHR == "L"]), size=3) + ggsave("lsegment_0_ARM_B16.pdf",width=8,height=6)

bps = ggplot(data=s.covs.10bp.melted[s.covs.10bp.melted$sample %in% levels(s.covs.10bp.melted$sample)[displ],], aes(x=BPS, y=coverage, group = sample, linetype=sample, colour = sample)) +
  geom_line() + scale_y_log10() + scale_linetype_manual(values = c(rep("solid", 7), rep("dashed", 7))) +
  scale_color_manual(values = c(colors,colors)) + theme_bw() + 
  theme(legend.direction ="horizontal",legend.position = "bottom") + guides(col=guide_legend(nrow=3)) + ggtitle("S Segment") 
bps + geom_vline(xintercept = primers$Start[primers$CHR == "S" ], aes(colour="black", linetype = "dotted", linewidth=5)) +
  annotate("text", x =primers$Start[primers$CHR == "S"], y=9*10^6, label=gsub("RT_","",primers$Primer[primers$CHR == "S"]), size=3) + ggsave("ssegment_0_ARM_B16.pdf",width=8,height=6)

displ=serum

bp = ggplot(data=l.covs.10bp.melted[l.covs.10bp.melted$sample %in% levels(l.covs.10bp.melted$sample)[displ],], aes(x=BPS, y=coverage, group = sample, linetype=sample, colour = sample)) +
  geom_line() + scale_y_log10() + scale_linetype_manual(values = c(rep("solid", 7), rep("dashed", 7))) +
  scale_color_manual(values = c(colors,colors)) + theme_bw() + 
  theme(legend.direction ="horizontal",legend.position = "bottom") + guides(col=guide_legend(nrow=3)) + ggtitle("L Segment") 
bp + geom_vline(xintercept = primers$Start[primers$CHR == "L" ], aes(colour="black", linetype = "dotted", linewidth=5)) +
  annotate("text", x =primers$Start[primers$CHR == "L"], y=9*10^6, label=gsub("RT_","",primers$Primer[primers$CHR == "L"]), size=3) + ggsave("lsegment_serum.pdf",width=8,height=6)

bps = ggplot(data=s.covs.10bp.melted[s.covs.10bp.melted$sample %in% levels(s.covs.10bp.melted$sample)[displ],], aes(x=BPS, y=coverage, group = sample, linetype=sample, colour = sample)) +
  geom_line() + scale_y_log10() + scale_linetype_manual(values = c(rep("solid", 7), rep("dashed", 7))) +
  scale_color_manual(values = c(colors,colors)) + theme_bw() + 
  theme(legend.direction ="horizontal",legend.position = "bottom") + guides(col=guide_legend(nrow=3)) + ggtitle("S Segment") 
bps + geom_vline(xintercept = primers$Start[primers$CHR == "S" ], aes(colour="black", linetype = "dotted", linewidth=5)) +
  annotate("text", x =primers$Start[primers$CHR == "S"], y=9*10^6, label=gsub("RT_","",primers$Primer[primers$CHR == "S"]), size=3) + ggsave("ssegment_serum.pdf",width=8,height=6)



primers$Y=10
primers$YE=10

# geom_text(data= primers[primers$CHR == "L",], aes(x = Start, y = Start, label= Primer))

bpS = ggplot(data=s.covs.50bp.melted, aes(x=BPS, y=coverage, group = sample, linetype=sample,colour = sample)) +
  geom_line() + scale_y_log10() + scale_linetype_manual(values = c(rep("solid", 6), rep("dotted", 6))) +
  scale_color_manual(values = c(brewer.pal(6, "Set1"), brewer.pal(6,"Accent"),brewer.pal(8, "Dark2"))) + theme_bw() + 
  guides(linetype=guide_legend(ncol=2)) + ggtitle("S Segment") + ggsave("ssegment.pdf",width=8,height=6)

library(ggplot2)
library(RColorBrewer)
library("reshape2")
#library("Gviz")
library("ggbio")
library("GenomicRanges")
library(rtracklayer)
library(biovizBase)
library(gridExtra)
setwd("/Volumes//vetgrid01/LCMV_project/Run_0277/Mappings/RF_coverages/")
rf.covs=read.table("rf_all.coverages",header=T)
gff_virus <- import.gff("/Volumes/vetgrid01/LCMV_project/References/viruses_short.gff")

length(gff_virus)
ranges(gff_virus)
seqnames(gff_virus)
str(gff_virus)
covs.long=melt(rf.covs,id.vars=c("CHROM","BPS"),value.name="COVERAGE",variable.name="SAMPLE_STRAND")
covs.long$SAMPLE=as.factor(sub("_[rf]+","",covs.long$SAMPLE_STRAND,perl=TRUE))
covs.long$STRAND=as.factor(sub("S_\\d+_","",covs.long$SAMPLE_STRAND,perl=TRUE))
str(covs.long)
gff_virus_short = gff_virus[! (a$type %in% c("databank_entry","CDS"))]
seqlengths(gff_virus_short) <- c(7229,3377)
levels(mcols(gff_virus_short)$type)=c(NA,"5pUTR","gene",NA,"IGR","3pUTR")
mcols(gff_virus_short)$ID=as.character(mcols(gff_virus_short)$type)
gnms=mcols(gff_virus_short)$gene
mcols(gff_virus_short)$ID[! is.na(gnms)] = gnms[! is.na(gnms)]

covs = ggplot(covs.long,aes(x=BPS, y=COVERAGE,colour=STRAND,group=STRAND) ) +  facet_grid( SAMPLE ~ CHROM , scales = "free")
covs = covs + geom_smooth(method="loess", span=0.05, se=FALSE)

covS=ggplot(covs.long[covs.long$CHROM=="S",],aes(x=BPS, y=COVERAGE,colour=STRAND,group=STRAND)) +
  geom_smooth(method="loess", span=0.05, se=FALSE) + facet_grid( SAMPLE ~ . , scales = "free")
gff_virus_short_S=gff_virus_short[seqnames(gff_virus_short) == "S"]
trackS=ggplot(gff_virus_short_S, aes(color = type, fill=type)) + 
  geom_alignment( label=FALSE, range.geom = "arrowrect",gap.geom = "chevron") + 
  annotate("text",x=(start(gff_virus_short_S)+end(gff_virus_short_S))/2,y=1,label=as.character(gff_virus_short_S$ID),size=3,angle = 25) +
  guides(fill=FALSE,color=FALSE)
tracks(covS,trackS,heights=c(13,1),xlim=c(0,seqlengths(gff_virus_short)["S"]),title="S")

library(ggrepel)
covL=ggplot(covs.long[covs.long$CHROM=="L",],aes(x=BPS, y=COVERAGE,colour=STRAND,group=STRAND)) + 
  geom_smooth(method="loess", span=0.05, se=FALSE)+ facet_grid( SAMPLE ~ . , scales = "free")
gff_virus_short_L=gff_virus_short[seqnames(gff_virus_short) == "L"]
trackL=ggplot(gff_virus_short_L, aes(color = type, fill=type)) +
  geom_alignment( label=FALSE, range.geom = "arrowrect",gap.geom = "chevron") + 
  annotate("text",x=(start(gff_virus_short_L)+end(gff_virus_short_L))/2,y=1,label=as.character(gff_virus_short_L$ID),size=3,angle = 25) +
  guides(fill=FALSE,color=FALSE)
tracks(covL,trackL,heights=c(13,1),xlim=c(0,seqlengths(gff_virus_short)["L"]),title="L")

geom_text_repel(data=as.data.frame(mcols(gff_virus_short_L)),x=(start(gff_virus_short_L)+end(gff_virus_short_L))/2,y=1,label=as.character(gff_virus_short_L$ID))



ggplot(gff_virus_short_S) + 
  geom_alignment(names.expr="type", label=FALSE, range.geom = "arrowrect",gap.geom = "chevron")
autoplot(gff_virus_short_S,names.expr="gene_name")
trackb=ggplot(gff_virus_short) + geom_alignment(aes(color = strand,fill=strand),  range.geom = "arrowrect",gap.geom = "chevron")+ facet_grid(~ seqnames,scales ="free_x")
tracks(covs,trackb,heights=c(12,1))
trackal=ggplot(gff_virus_short) + geom_alignment( range.geom = "arrowrect", gap.geom = "chevron")+ facet_grid(~ seqnames,scales ="free") 

#tracks(covs, trackal,heights=c(12,1),xlim=ggplot_build(covs)$layout$panel_ranges[[1]]$x.range)
tracks(covs, trackal,heights=c(12,1))
fixed(covs) = TRUE
#ggplot(covs.long,aes(x=BPS, y=COVERAGE,colour=SAMPLE,group=STRAND) ) +  geom_line() + facet_grid( . ~ CHROM , scales = "free")
ggplot_build(covs)$layout$panel_ranges

library(GenomicFeatures)
makeTxDbFromGFF(gff_virus)
library(vcfR)
library(VariantAnnotation)
ref_fa=FaFile("/Volumes/vetgrid01/LCMV_project/References/viruses_short.fasta")
txdb=makeTxDbFromGFF("/Volumes/vetgrid01/LCMV_project/References/viruses_short.gff",format="auto",dataSource="ensemble",
                     chrominfo =  Seqinfo(seqnames = c("L","S"),seqlength=c(7229,3377), isCircular=c(FALSE, FALSE)),taxonomyId=11623)
coding <- predictCoding(vcf, txdb ,seqSource = ref_fa)

vcf_file="/Volumes/vetgrid01/LCMV_project/Run_0277/Mappings/Maxcov10K/Freebayes/run_0277.maxcov10K_co_freebayes_nocmpl_0.01_anno.vcf"
vcf_file="/Volumes/vetgrid01/LCMV_project/Run_0277/Mappings/Maxcov10K/Lofreq2/ClipOverlaps/lofreq2_all_samp_bed_norm_0.001_all.vcf.gz"
vcf_file="/Volumes/vetgrid01/LCMV_project/Run_0277/Mappings/Maxcov10K/VPhaser2/all_out_dps_anno.vcf"
vcf <- read.vcfR( vcf_file, verbose = FALSE )
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
quartz(width=10,height=10)
layout(matrix(1:7, ncol = 1), widths = 1, respect = FALSE)
chrom_lens=c(3377,7229)
names(chrom_lens) = c("S","L")
for (i in 2:7) {
  print(paste("plotting sample",samples[i]))
  reg <- plot_AFs_depths(vcf,chrom,samples[i],chrom_len = chrom_lens[chrom])
}

par(mar = c(3, 4, 0, 4))
plot(1,10,pch=NULL,xlim=c(1,chrom_lens[chrom]),ylim=c(-1,1),ylab=NA,axes=FALSE,xlab=chrom,bty="n")
axis(side=1,at=0:(chrom_lens[chrom]/1000)*1000)
mtext(chrom, side=1, line=2)
for( i in 1:length(vir_clean[seqnames(vir_clean) == chrom])){plot_element(vir_clean[seqnames(vir_clean) == chrom][i],chrom_lens[chrom])}


vcf_fb="/Volumes/vetgrid01/LCMV_project/Run_0277/Mappings/Maxcov10K/Freebayes/run_0277.maxcov10K_co_freebayes_nocmpl_0.01_anno_un.vcf"
vcf_lf="/Volumes/vetgrid01/LCMV_project/Run_0277/Mappings/Maxcov10K/Lofreq2/ClipOverlaps/lofreq2_all_samp_bed_norm_0.001_all.vcf.gz"
vcf_vp="/Volumes/vetgrid01/LCMV_project/Run_0277/Mappings/Maxcov10K/VPhaser2/all_out_dps_anno.vcf"
vcf_vs="/Volumes/vetgrid01/LCMV_project/Run_0277/Mappings/Maxcov10K/Varscan2/bam_varscan_afs_0.01_anno.vcf"
fb_vcf = load_vcf_file(vcf_fb)
vs_vcf = load_vcf_file(vcf_vs)


samps_1=c("S_1","S_2","S_3","S_4","S_5","S_6","S_7","S_8","S_9" )
samps_1=c("S_10","S_11","S_12","S_13","S_14","S_15")
samps = list("s1_9" = c("S_1","S_2","S_3","S_4","S_5","S_6","S_7","S_8","S_9" ),
             "s10_15" = c("S_10","S_11","S_12","S_13","S_14","S_15"))
chrom="S"
vcf = fb_vcf
bn="freebayes"
vcf = vs_vcf
bn="varscan2"
samp=c("s1_9")
setwd("/Volumes/vetgrid01/LCMV_project/Run_0277/Mappings/Maxcov10K/Lofreq2/")
for ( chrom in c("L","S")) {
  for (samp in names(samps)) {
    samps_1=samps[[samp]]
    layout(matrix(1:(length(samps_1)+1), ncol = 1), widths = 1, respect = FALSE)
    for (i in samps_1) {
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
    pdfn=paste(bn,samp,"chrom",chrom,sep="_")
    dev.copy2pdf(file=paste0(pdfn,".pdf"))
    
  }
}




layout(matrix(1:(length(samps_1)+1), ncol = 1), widths = 1, respect = FALSE)
for (i in samps_1) {
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

plot_AFs_depths(vcf,chrom,samples[i])
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
  afs$start = sapply(afs$xpos,function (x)  df_region$start[df_region$xpos == x ])
  points(afs[,xcol],afs$AF,col=cols[afs$ANNO],cex=1,pch=symb[afs$ALL]) 
  par(new = T)
  plot(df_region[,xcol], df_region$DP, type="l", axes=F, log="y", xlab=NA, ylab=NA,col="darkgrey")
  axis(side = 4)
  mtext(side = 4, line = 2.75, 'Depth',cex=0.75)
  #par(mar = c(0, 4.1, 4.1, 2.1))
  return(df_region[,c("start","xpos")])
  
}

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
  afs=afs[! is.na(afs[,3]),]
  afs=data.frame(afs)
  colnames(afs)=c("xpos","ALL","AF")
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



