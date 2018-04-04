# get file name 
args <- commandArgs(TRUE)
# this script requires the libraries ggplot2, grid, hexbin, and grid.extra to be installed
# you should be able to do this by:
# install.packages(c("ggplot2","gridExtra","hexbin"),dependencies="Depends")
# only takes one argument, the tab file with the metrics for each SNP
# the tab delimited metrics file should have the following fields:
# Chromosome Position Alleles FS SB RTD ATD aRPB Ref
file_name <- args[1]
base_name <- args[2]
# getting rid of ending
if(is.null(base_name)) {
	base_name <- sub("([.][^.]+)$",'',basename(file_name),perl=T)
}
value_tab=read.table(file_name,header=F,na.strings = c("NA","nan"))
if ( ncol(value_tab) == 10) { 
    colnames(value_tab)=c("CHR","BPS","ALL","FS","SB","RTD","ATD","TDB","aRPB","Ref")
} else {
    colnames(value_tab)=c("CHR","BPS","ALL","FS","SB","RTD","ATD","TDB","aRPB")
    value_tab$Ref=0
}
library(ggplot2)
### ggplot2 had a change of synthax around version 0.9.1 (or so), so we need to check the version and set the synthax accordingly:
opt.val="opts"
if (packageVersion("ggplot2") >= '0.9.2'){
opt.val="labs"
}

if(length(value_tab$Ref[value_tab$Ref==1]) == 0){
    # only create simple plots with all snps 
    mybreaks=c(1,10,50,100,250,500,1000,5000,10000)
    xlimit=c(0,600)
    #pdf(file=paste(base_name,"_FS_SB.pdf",sep=""))
    good=ggplot(data=value_tab) + xlim(xlimit) +
        stat_binhex(data=value_tab,aes(x=FS, y=SB),bins=40) +
            scale_fill_gradientn(colours=c("red","green","blue"),trans="log10", name = "count",breaks= mybreaks, labels=mybreaks, na.value=NA) + get(opt.val)(title="SNPs") #+ labs(title="SNPs")
    #dev.off()
    ggsave(good,file=paste(base_name,"_FS_SB.pdf",sep=""))
    #pdf(file=paste(base_name,"_aRPB_ATD.pdf",sep=""))
    good=ggplot(data=value_tab) + stat_binhex(data=value_tab,aes(x=aRPB, y=ATD),binwidth=c(0.05,1)) +
        scale_fill_gradientn(colours=c("red","green","blue"),trans="log", name = "count",breaks= mybreaks, labels=mybreaks, na.value=NA) + get(opt.val)(title="SNPs")
    #dev.off()
    ggsave(good,file=paste(base_name,"_aRPB_ATD.pdf",sep=""))
    #pdf(file=paste(base_name,"_aRPB_ATD.pdf",sep=""))
    xlimit=c(0,100)
    good=ggplot(data=value_tab) + xlim(xlimit) + stat_binhex(data=value_tab,aes(x=TDB, y=ATD),binwidth=c(2.5,1.0)) +
        scale_fill_gradientn(colours=c("red","green","blue"),trans="log", name = "count",breaks= mybreaks, labels=mybreaks, na.value=NA) + get(opt.val)(title="SNPs")
    #dev.off()
    ggsave(good,file=paste(base_name,"_TDB_ATD.pdf",sep=""))
    quit(save="no")
}

library(grid)
library(gridExtra)

good_snps = value_tab$Ref == 1
num_good=length(value_tab$Ref[good_snps])
bad_snps =  value_tab$Ref == 0
num_bad=length(value_tab$Ref[bad_snps])
mybreaks=c(1,10,50,100,250,500,1000,5000,10000)
xlimit=c(0,600)
pdf(file=paste(base_name,"_FS_SB_double.pdf",sep=""))
good=ggplot(data=value_tab) + xlim(xlimit) +
     stat_binhex(data=value_tab[good_snps,],aes(x=FS, y=SB),bins=40) +
    scale_fill_gradientn(colours=c("red","green","blue"),trans="log", name = "count",breaks= mybreaks, labels=mybreaks, na.value=NA) + get(opt.val)(title="SNPs in Reference")
bad=ggplot(data=value_tab) + xlim(xlimit) +
     stat_binhex(data=value_tab[bad_snps,],aes(x=FS, y=SB),bins=40) + scale_fill_gradientn(colours=c("red","green","blue"),trans="log", name = "count",breaks= mybreaks, labels=mybreaks, na.value=NA) + get(opt.val)(title="SNPs not in Reference")
grid.arrange(good,bad,ncol=1)  
dev.off()

pdf(file=paste(base_name,"_aRPB_ATD_double.pdf",sep=""))
good=ggplot(data=value_tab) +
     stat_binhex(data=value_tab[good_snps,],aes(x=aRPB, y=ATD),binwidth=c(0.5,1)) +
    scale_fill_gradientn(colours=c("red","green","blue"),trans="log", name = "count",breaks= mybreaks, labels=mybreaks, na.value=NA) + get(opt.val)(title="SNPs in Reference")
bad=ggplot(data=value_tab) +
     stat_binhex(data=value_tab[bad_snps,],aes(x=aRPB, y=ATD),binwidth=c(0.5,1)) + scale_fill_gradientn(colours=c("red","green","blue"),trans="log", name = "count",breaks= mybreaks, labels=mybreaks, na.value=NA) + get(opt.val)(title="SNPs not in Reference")
grid.arrange(good,bad,ncol=1)     
dev.off()

pdf(file=paste(base_name,"_TDB_ATD_double.pdf",sep=""))
good=ggplot(data=value_tab) +
     stat_binhex(data=value_tab[good_snps,],aes(x=TDB, y=ATD),binwidth=c(1.0,1)) +
    scale_fill_gradientn(colours=c("red","green","blue"),trans="log", name = "count",breaks= mybreaks, labels=mybreaks, na.value=NA) + get(opt.val)(title="SNPs in Reference")
bad=ggplot(data=value_tab) +
     stat_binhex(data=value_tab[bad_snps,],aes(x=TDB, y=ATD),binwidth=c(1.0,1)) + scale_fill_gradientn(colours=c("red","green","blue"),trans="log", name = "count",breaks= mybreaks, labels=mybreaks, na.value=NA) + get(opt.val)(title="SNPs not in Reference")
grid.arrange(good,bad,ncol=1)     
dev.off()


# just to be nice remove the crappy plot variables
rm("good","bad")
gc()

# create ROC like curves for retrieval of SNPs
# using thresholds for the PHRED scaled fisher test P value for strand imbalance (FS) and the fraction of the less common strand supporting the alternative allele (SB)
pdf(file=paste(base_name,"_FS_SB_SNPs_filtered.pdf",sep=""))
par(mar=c(5, 4, 4, 8) + 0.1)
fs_vals=c(5,10,20,30,60,100,150,200,250)
good_vals = sapply(fs_vals,function (x) length(value_tab$FS[good_snps & value_tab$FS >= x] ))
bad_vals =  sapply(fs_vals,function (x) length(value_tab$FS[bad_snps & value_tab$FS >= x] ))
plot(fs_vals, good_vals/num_good,type="l", col="blue", lty=1,ylim=c(0,0.2),xlab="FS",ylab="fraction filtered")
lines(fs_vals, bad_vals/num_bad,type="l", col="red", lty=1)
SB_vals=c(0.2,0.1,0.05,0.01)
ltypes=c(2,3,4,5)
for(i in 1:4){
good_vals = sapply(fs_vals,function (x) length(value_tab$FS[good_snps & value_tab$FS >= x & value_tab$SB <= SB_vals[i]] ))
bad_vals =  sapply(fs_vals,function (x) length(value_tab$FS[bad_snps & value_tab$FS >= x & value_tab$SB <=  SB_vals[i]] ))
lines(fs_vals, good_vals/num_good,type="l", col="blue", lty=ltypes[i])
lines(fs_vals, bad_vals/num_bad,type="l", col="red", lty=ltypes[i])
}
axis(4,col.axis="red",at=seq(0,0.2,by=0.05),labels=round(seq(0,0.2,by=0.05)*num_bad,digits=0))
axis(4,col.axis="blue",outer=T,at=seq(0,0.2,by=0.05),labels=round(seq(0,0.2,by=0.05)*num_good,digits=0),line=-3.5)
title(main="SNPs filtered out")
legend("topright",legend=c("in reference SNP set","not in reference SNP set (novel or false)","SB <= 1.0","SB <= 0.2","SB <= 0.1","SB <= 0.05","SB <= 0.01"),col=c("blue","red","grey10","grey10","grey10","grey10","grey10"),lty=c(1,1,1,2,3,4,5))
dev.off()
# thresholds for the absolute value of the Read posiotn bias (aRPB) and the average alternative tail distance (ATD)
pdf(file=paste(base_name,"_aRPB_ATD_SNPs_filtered.pdf",sep=""))
par(mar=c(5, 4, 4, 8) + 0.1)
aRPB_vals=c(0,1,2,3,3.5,4,4.5,5,7.5,10,15)
good_vals = sapply(aRPB_vals,function (x) length(value_tab$aRPB[good_snps & value_tab$aRPB >= x] ))
bad_vals =  sapply(aRPB_vals,function (x) length(value_tab$aRPB[bad_snps & value_tab$aRPB >= x] ))
plot(aRPB_vals, good_vals/num_good,type="l", col="blue", lty=1,ylim=c(0,0.25),xlab="aRPB",ylab="fraction filtered")
lines(aRPB_vals, bad_vals/num_bad,type="l", col="red", lty=1)
ATD_vals=c(15,10,7.5,5,2.5,1)
ltypes=c(2,3,4,5,6,2)
for(i in 1:4){
    good_vals = sapply(aRPB_vals,function (x) length(value_tab$aRPB[good_snps & value_tab$aRPB >= x & value_tab$ATD <= ATD_vals[i]] ))
    bad_vals =  sapply(aRPB_vals,function (x) length(value_tab$aRPB[bad_snps & value_tab$aRPB >= x & value_tab$ATD <=  ATD_vals[i]] ))
    lines(aRPB_vals, good_vals/num_good,type="l", col="blue", lty=ltypes[i])
    lines(aRPB_vals, bad_vals/num_bad,type="l", col="red", lty=ltypes[i])
}
axis(4,col.axis="red",at=seq(0,0.25,by=0.05),labels=round(seq(0,0.25,by=0.05)*num_bad,digits=0))
axis(4,col.axis="blue",outer=T,at=seq(0,0.25,by=0.05),labels=round(seq(0,0.25,by=0.05)*num_good,digits=0),line=-3.5)
title(main="SNPs filtered out")
legend("topright",legend=c("in reference SNP set","not in reference SNP set (novel or false)","ATD <= 15","ATD <= 10","ATD <= 7.5","ATD <= 5","ATD <= 2.5","ATD <= 1.0"),col=c("blue","red","grey10","grey10","grey10","grey10","grey10","grey10"),lty=c(1,1,2,3,4,5,6,2))
dev.off()

pdf(file=paste(base_name,"_TDB_ATD_SNPs_filtered.pdf",sep=""))
par(mar=c(5, 4, 4, 8) + 0.1)
TDB_vals=c(0,1,5,10,15,20,25,30,35,40,50,75,100,120)
good_vals = sapply(TDB_vals,function (x) length(value_tab$TDB[good_snps & value_tab$TDB >= x] ))
bad_vals =  sapply(TDB_vals,function (x) length(value_tab$TDB[bad_snps & value_tab$TDB >= x] ))
plot(TDB_vals, good_vals/num_good,type="l", col="blue", lty=1,ylim=c(0,0.25),xlab="TDB",ylab="fraction filtered")
lines(TDB_vals, bad_vals/num_bad,type="l", col="red", lty=1)
ATD_vals=c(15,10,7.5,5,2.5,1)
ltypes=c(2,3,4,5,6,2)
for(i in 1:4){
    good_vals = sapply(TDB_vals,function (x) length(value_tab$TDB[good_snps & value_tab$TDB >= x & value_tab$ATD <= ATD_vals[i]] ))
    bad_vals =  sapply(TDB_vals,function (x) length(value_tab$TDB[bad_snps & value_tab$TDB >= x & value_tab$ATD <=  ATD_vals[i]] ))
    lines(TDB_vals, good_vals/num_good,type="l", col="blue", lty=ltypes[i])
    lines(TDB_vals, bad_vals/num_bad,type="l", col="red", lty=ltypes[i])
}
axis(4,col.axis="red",at=seq(0,0.25,by=0.05),labels=round(seq(0,0.25,by=0.05)*num_bad,digits=0))
axis(4,col.axis="blue",outer=T,at=seq(0,0.25,by=0.05),labels=round(seq(0,0.25,by=0.05)*num_good,digits=0),line=-3.5)
title(main="SNPs filtered out")
legend("topright",legend=c("in reference SNP set","not in reference SNP set (novel or false)","ATD <= 15","ATD <= 10","ATD <= 7.5","ATD <= 5","ATD <= 2.5","ATD <= 1.0"),col=c("blue","red","grey10","grey10","grey10","grey10","grey10","grey10"),lty=c(1,1,2,3,4,5,6,2))
dev.off()

