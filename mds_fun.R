install.packages("rafalib")
install.packages("downloader")
install.packages("devtools")
library(rafalib)
library(downloader)

library(devtools)
install_github("genomicsclass/tissuesGeneExpression")
library(tissuesGeneExpression)
data(tissuesGeneExpression)
dim(e) ##e contains the expression data

imagemat(e[1:200,1:50])
colind <- tissue%in%c("kidney","colon","liver")
mat <- e[,colind]
group <- factor(tissue[colind])
group
s <- svd(mat-rowMeans(mat))
PC1 <- s$d[1]*s$v[,1]
PC2 <- s$d[2]*s$v[,2]

barplot(s$d/sum(s$d))
plot(1:length(s$d),cumsum(s$d)/sum(s$d))

library(edgeR)
a <- plotMDS(mat, ndim = 2)
b <- cmdscale(a$distance.matrix, eig=TRUE)
b <- cmdscale(a$distance.matrix, k=98,eig=TRUE)
barplot(sort(abs(b$eig),decreasing =T)/sum(abs(b$eig)))

data("eurodist")
str(eurodist)
euro.msd <- cmdscale(eurodist,eig=T)
plot(euro.msd$points[,2]~euro.msd$points[,1] )
text(euro.msd$points[,2]~euro.msd$points[,1], labels = labels(eurodist), pos = 3)
barplot(sort(abs(euro.msd$eig),decreasing =T)/sum(abs(euro.msd$eig)))


str(UScitiesD)
us.msd <- cmdscale(UScitiesD,eig=T)
plot(us.msd$points[,1]~us.msd$points[,2] )
text(us.msd$points[,1]~us.msd$points[,2], labels = labels(UScitiesD), pos = 3)
barplot(sort(abs(us.msd$eig),decreasing =T)/sum(abs(us.msd$eig)))

library(foreign)
cong <- read.dta("~/Downloads/hou113kh.dta")
str(cong)
votes <- cong[,10:ncol(cong)]
votes[votes == 9] = 0
votes[votes == 6] = -1
summary(votes)
cong.msd <- cmdscale(dist(votes),eig=T)
party <- ifelse(cong$party == 100, "blue", "red") 
plot(cong.msd$points[,1],cong.msd$points[,2], col = party)
barplot(sort(abs(cong.msd$eig[1:10]),decreasing =T)/sum(abs(cong.msd$eig)),xlab = "% variation 10 highest eigenvalues")
plotMDS(votes)
library(adegenet)
data(nancycats)
is.genind(nancycats)
data(microbov)
obj <- genind2genpop(microbov)
ca1 <- dudi.coa(tab(obj),scannf=FALSE,nf=3)
barplot(ca1$eig,main="Correspondance Analysis eigenvalues",
        col=heat.colors(length(ca1$eig)))
s.label(ca1$li, sub="CA 1-2",xax=1,yax=2,csub=2)
add.scatter.eig(ca1$eig,nf=3,xax=1,yax=2,posi="bottomright")

s.label(ca1$li,xax=2,yax=3, sub="CA 2-3",csub=2)
add.scatter.eig(ca1$eig,nf=3,xax=2,yax=3,posi="bottomright")
