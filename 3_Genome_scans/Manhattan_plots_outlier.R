library(dplyr)
library(ggplot2)
library(doBy)
library(wesanderson)
library(expss)
library(maditr)

# load pbs results
setwd("~/FST")
pbs_raw <- read.table('FST_North_Southres', header=T, stringsAsFactors=F)

#Select autosomes only
chrom_list <- read.table('Sca_chrom_list.txt', header=T, stringsAsFactors=F)

pbs<-subset(pbs_raw,(chr %in% chrom_list$chr))

df_pbs<-merge(pbs, chrom_list , by="chr")


# order by chr number
res = df_pbs[which(df_pbs[6] %in% 'sca1'),]
for(i in 1:34){
  res = rbind(res, df_pbs[which(df_pbs[,7] %in% paste0('sca',i)),])
}
res = cbind(res,01:nrow(res))


# find the middle points for each chr (for naming)
#run:

table(res[,7])


l =c(11477,11273,10424,10251,9590,9526,9171,8083,7697,7595,7281,7232,7208,6892,6659,6252,6175,6103,5979,5977,5866,5657,5443,5398,5190,4879,4455,4365,4233,4072,4069,3765)
ind = c(5738)
for (j in 2:length(l)){
  ind[j] = ind[j-1]+round(l[j-1]/2)+round(l[j]/2)
}

# plot
pdf('Moose_North_South_FST.pdf', height=7, width=14)
plot(NULL, pch=16, cex.axis=1.5, cex.lab = 2, ylab='Fst',ylim=c(0,1), xlim=c(0,nrow(res)), xlab='',xaxt='n', main='Top 0.1% Fst values')
axis(1, at=ind, labels=unique(res[,7]), cex.axis=0.85, las=2)
points(res[which(res[,7] %in% c('sca1','sca3','sca5','sca7','sca9','sca11','sca13','sca15','sca17','sca19','sca22','sca25','sca27','sca29','sca31','sca33')),8], res[which(res[,7] %in% c('sca1','sca3','sca5','sca7','sca9','sca11','sca13','sca15','sca17','sca19','sca22','sca25','sca27','sca29','sca31','sca33')),5], pch=1, cex=0.4, col='black',cex.axis=0.8)
points(res[which(res[,7] %in% c('sca2','sca4','sca6','sca8','sca10','sca12','sca14','sca16','sca18','sca20','sca23','sca26','sca28','sca30','sca32','sca34')),8], res[which(res[,7] %in% c('sca2','sca4','sca6','sca8','sca10','sca12','sca14','sca16','sca18','sca20','sca23','sca26','sca28','sca30','sca32','sca34')),5], pch=1, cex=0.4, col='grey',cex.axis=0.8)

# select top hit level
q = quantile(as.numeric(res[,5]),0.999)

#points(res[res[,9]>=q,12], res[res[,9]>=q,9], pch=16, cex=0.3, col='blue')
points(res[res[,5]>=q,8], res[res[,5]>=q,5], pch=16, cex=0.8, col = wes_palette("Zissou1",n=1))
abline(h=q, col="red",lwd=1, lty=2)
dev.off()


#Z transformation of FST
#https://github.rcac.purdue.edu/MarkRChristieGroup/steelhead-poolse

#ZFst = (Fst – μ Fst)/ σ Fst, where Fst is the Fst in a window, μ Fst is an average Fst over all windows, and σ Fst
#is a standard deviation of Fst values of all windows 

mean_FST<-mean(df_pbs$Fst)
SD_FST<-sd(df_pbs$Fst)
FST<-df_pbs$Fst
ZFST<- (FST-mean_FST)/SD_FST

min(ZFST)
max(ZFST)
res2<-data.frame(res,ZFST)


# plot
pdf('Moose_North_South_Z_FST.pdf', height=7, width=14)
plot(NULL, pch=16, cex.axis=1.5, cex.lab = 2, ylab='Z(Fst)', xlim=c(0,nrow(res2)),ylim=c(0,12), xlab='',xaxt='n', main='')
axis(1, at=ind, labels=unique(res2[,7]), cex.axis=0.85, las=2)
points(res2[which(res2[,7] %in% c('sca1','sca3','sca5','sca7','sca9','sca11','sca13','sca15','sca17','sca19','sca22','sca25','sca27','sca29','sca31','sca33')),8], res2[which(res2[,7] %in% c('sca1','sca3','sca5','sca7','sca9','sca11','sca13','sca15','sca17','sca19','sca22','sca25','sca27','sca29','sca31','sca33')),9], pch=1, cex=0.4, col='black',cex.axis=0.8)
points(res2[which(res2[,7] %in% c('sca2','sca4','sca6','sca8','sca10','sca12','sca14','sca16','sca18','sca20','sca23','sca26','sca28','sca30','sca32','sca34')),8], res2[which(res2[,7] %in% c('sca2','sca4','sca6','sca8','sca10','sca12','sca14','sca16','sca18','sca20','sca23','sca26','sca28','sca30','sca32','sca34')),9], pch=1, cex=0.4, col='grey',cex.axis=0.8)
points(res[res2[,9]>=q,8], res2[res2[,9]>=q,9], pch=16, cex=0.8, col = wes_palette("Zissou1",n=1))

# select top hit level (i.e. >=5)
q=5
abline(h=q, col="red",lwd=1, lty=2)
dev.off()

##extract hits with q >=5

nrow(res2[res2[,9]>=q,])
#683

write.table(res2[res2[,9]>=q,], file='ZFST_5_North_South.bed', row.names=F, col.names=F, quote=F, sep='\t') 


# Based on this *bed file, one can then cross reference it with the annotation 
#(i.e. bedtools intersect -a MAM_SEL/sorted.Loxodonta_africana.loxAfr4.gff3.gz -b top01_MAM_ELE.bed > matches) 
#and you get a list of annotations in those top hits.
