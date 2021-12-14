install.packages("ks")
library(ks)
library(RColorBrewer)
library(zoo)


dir2 = "/Metapopulation_model/Figures/"
dir4 = '/Metapopulation_model/metapopulation_polio/'

oneVillage64000_TBP = read.csv(paste0(dir4,"time_between_pcases_64000reintRate_0.001000migRate_0.000000_paper.csv"),sep=",",header = FALSE,col.names = paste0("V",seq_len(9000)), fill = TRUE)
oneVillage64000_extInt = read.table(paste0(dir4,"extinction_interval_64000reintRate_0.001000migRate_0.000000_paper.csv"))
twoVillage32000_TBP = read.csv(paste0(dir4,"time_between_pcases_3200032000reintRate_0.001000migRate_0.000000_paper.csv"),sep=",",header = FALSE,col.names = paste0("V",seq_len(2000)), fill = TRUE)
twoVillage32000_extInt = read.table(paste0(dir4,"extinction_interval_3200032000reintRate_0.001000migRate_0.000000_paper.csv"))
fourVillage16000_TBP = read.csv(paste0(dir4,"time_between_pcases_16000160001600016000reintRate_0.001000migRate_0.000000_paper.csv"),sep=",",header = FALSE,col.names = paste0("V",seq_len(500)), fill = TRUE)
fourVillage16000_extInt = read.table(paste0(dir4,"extinction_interval_16000160001600016000reintRate_0.001000migRate_0.000000_paper.csv"))
eightVillage8000_TBP = read.csv(paste0(dir4,"time_between_pcases_80008000800080008000800080008000reintRate_0.001000migRate_0.000000_paper.csv"),sep=",",header = FALSE,col.names = paste0("V",seq_len(500)), fill = TRUE)
eightVillage8000_extInt = read.table(paste0(dir4,"extinction_interval_80008000800080008000800080008000reintRate_0.001000migRate_0.000000_paper.csv"))
sixteenVillage4000_TBP = read.csv(paste0(dir4,"time_between_pcases_4000400040004000400040004000400040004000400040004000400040004000reintRate_0.001000migRate_0.000000_paper.csv"),sep=",",header = FALSE,col.names = paste0("V",seq_len(500)), fill = TRUE)
sixteenVillage4000_extInt = read.table(paste0(dir4,"extinction_interval_4000400040004000400040004000400040004000400040004000400040004000reintRate_0.001000migRate_0.000000_paper.csv"))
thirtytwoVillage2000_TBP = read.csv(paste0(dir4,"time_between_pcases_20002000200020002000200020002000200020002000200020002000200020002000200020002000200020002000200020002000200020002000200020002000reintRate_0.001000migRate_0.000000_paper.csv"),sep=",",header = FALSE,col.names = paste0("V",seq_len(500)), fill = TRUE)
thirtytwoVillage2000_extInt = read.csv(paste0(dir4,"extinction_interval_32x2000reintRate_0.001000migRate_0.000000_paper.csv"))



#remove nas
new <- unlist(oneVillage64000_extInt)[!is.na(unlist(oneVillage64000_extInt))]
new1<-unlist(twoVillage32000_extInt)[!is.na(unlist(twoVillage32000_extInt))]
new2<-unlist(fourVillage16000_extInt)[!is.na(unlist(fourVillage16000_extInt))]
new3<-unlist(eightVillage8000_extInt)[!is.na(unlist(eightVillage8000_extInt))]
new4<-unlist(sixteenVillage4000_extInt)[!is.na(unlist(sixteenVillage4000_extInt))]
new5<-unlist(thirtytwoVillage2000_extInt)[!is.na(unlist(thirtytwoVillage2000_extInt))]

#density estimates of log of data
k = kde(x=as.matrix(new),positive=TRUE)
k1 = kde(x=as.matrix((new1)),positive=TRUE)
k2 = kde(x=as.matrix((new2)),positive=TRUE)
k3 = kde(x=as.matrix((new3)),positive=TRUE)
k4 = kde(x=as.matrix((new4)),positive=TRUE)
k5 = kde(x=as.matrix((new5)),positive=TRUE)

#plot fit to hist
png(paste0(dir2,"hist_gauss_curve_fit.png"), width=1400, height=800, res=150) 
hist(new,freq=FALSE,xlab="Extinction intervals (years)",main="Gaussian fit to N=64,000 pop size")
lines(k$eval.points,k$estimate,col='red')
dev.off()

#plot smoothed curve on top of hist
hist(new,freq=FALSE)
lines(k$eval.points,k$estimate,col='red')
hist(new1,freq=FALSE)
lines(k1$eval.points,k1$estimate,col='red')
hist(new2,freq=FALSE)
lines(k2$eval.points,k2$estimate,col='red')
hist(new3,freq=FALSE)
lines(k3$eval.points,k3$estimate,col='red')
hist(new4,freq=FALSE)
lines(k4$eval.points,k4$estimate,col='red')
hist(new5,freq=FALSE)
lines(k5$eval.points,k5$estimate,col='red')


marker=colorRampPalette(c('orange','red','purple','royalblue'))(7)

#plot smooth curves on top of each other
png(paste0(dir2,"extInt_dist_noMig_multipatch_paper_v1.png"), width=1400, height=800, res=150) 
par(fig=c(0,1,0,1),mar=c(3.3,3.6,1.1,1.1),oma=c(1,0,0,3),mgp=c(2.2,0.45,0),tcl=-0.4,cex=1.2)
plot((k$eval.points),k$estimate,col='black',type='l',ylim=c(0,1.4),xlim=c(0,5),xlab='',ylab='Density',lwd=2)
lines((k1$eval.points),k1$estimate,col=marker[[1]],lwd=2)
lines((k2$eval.points),k2$estimate,col=marker[[2]],lwd=2)
lines((k3$eval.points),k3$estimate,col=marker[[4]],lwd=2)
lines((k4$eval.points),k4$estimate,col=marker[[6]],lwd=2)
lines((k5$eval.points),k5$estimate,col=marker[[7]],lwd=2)
#mtext("Density",side=2,outer=TRUE,cex=1.2)
mtext('Length of extinction interval (Years)',side=1,outer=TRUE,cex=1.2)
legend(x=1.5,y=1.2,legend=rev(c('1x64k','2x32k','4x16k','8x8k','16x4k','32x2k')),col=rev(c('black',marker[1:2],marker[4],marker[6:7])),lwd=2,bty='n')
par(new=TRUE)
par(mar=c(0,0,0,0))
par(fig=c(0.6,0.95,0.5,0.85))
plot((k$eval.points),k$estimate,col='black',type='l',ylim=c(0,0.04),xlim=c(3,5),xlab='',ylab='',lwd=2)
lines((k1$eval.points),k1$estimate,col=marker[[1]],lwd=2)
lines((k2$eval.points),k2$estimate,col=marker[[2]],lwd=2)
lines((k3$eval.points),k3$estimate,col=marker[[4]],lwd=2)
lines((k4$eval.points),k4$estimate,col=marker[[6]],lwd=2)
lines((k5$eval.points),k5$estimate,col=marker[[7]],lwd=2)
dev.off()



new <- unlist(oneVillage64000_TBP)[!is.na(unlist(oneVillage64000_TBP))]
new1<-unlist(twoVillage32000_TBP)[!is.na(unlist(twoVillage32000_TBP))]
new2<-unlist(fourVillage16000_TBP)[!is.na(unlist(fourVillage16000_TBP))]
new3<-unlist(eightVillage8000_TBP)[!is.na(unlist(eightVillage8000_TBP))]
new4<-unlist(sixteenVillage4000_TBP)[!is.na(unlist(sixteenVillage4000_TBP))]
new5<-unlist(thirtytwoVillage2000_TBP)[!is.na(unlist(thirtytwoVillage2000_TBP))]

#density estimates of log of data
k = kde(x=as.matrix(new),positive=TRUE)
k1 = kde(x=as.matrix((new1)),positive=TRUE)
k2 = kde(x=as.matrix((new2)),positive=TRUE)
k3 = kde(x=as.matrix((new3)),positive=TRUE)
k4 = kde(x=as.matrix((new4)),positive=TRUE)
k5 = kde(x=as.matrix((new5)),positive=TRUE)

#plot smoothed curve on top of hist
hist(new,freq=FALSE)
lines(k$eval.points,k$estimate,col='red')
hist(new1,freq=FALSE)
lines(k1$eval.points,k1$estimate,col='red')
hist(new2,freq=FALSE)
lines(k2$eval.points,k2$estimate,col='red')
hist(new3,freq=FALSE)
lines(k3$eval.points,k3$estimate,col='red')
hist(new4,freq=FALSE)
lines(k4$eval.points,k4$estimate,col='red')
hist(new5,freq=FALSE)
lines(k5$eval.points,k5$estimate,col='red')


#plot smooth curves on top of each other
png(paste0(dir2,"interInt_dist_noMig_multipatch_paper_v1.png"), width=1400, height=800, res=150) 
par(fig=c(0,1,0,1),mar=c(3.3,3.6,1.1,1.1),oma=c(1,0,0,3),mgp=c(2.2,0.45,0),tcl=-0.4,cex=1.2)
plot((k$eval.points),k$estimate,col='black',type='l',ylim=c(0,1),xlim=c(-.3,5),xlab='',ylab='Density',lwd=2)
lines((k1$eval.points),k1$estimate,col=marker[[1]],lwd=2)
lines((k2$eval.points),k2$estimate,col=marker[[2]],lwd=2)
lines((k3$eval.points),k3$estimate,col=marker[[4]],lwd=2)
lines((k4$eval.points),k4$estimate,col=marker[[6]],lwd=2)
lines((k5$eval.points),k5$estimate,col=marker[[7]],lwd=2)
#mtext("Density",side=2,outer=TRUE,cex=1.2)
mtext('Length of intercase interval (Years)',side=1,outer=TRUE,cex=1.2)
legend(x=1,y=0.8,legend=rev(c('1x64k','2x32k','4x16k','8x8k','16x4k','32x2k')),col=rev(c('black',marker[1:2],marker[4],marker[6:7])),lwd=2,bty='n')
par(new=TRUE)
par(mar=c(0,0,0,0))
par(fig=c(0.6,0.95,0.5,0.85))
plot((k$eval.points),k$estimate,col='black',type='l',ylim=c(0,0.01),xlim=c(3,5),xlab='',ylab='',lwd=2)
lines((k1$eval.points),k1$estimate,col=marker[[1]],lwd=2)
lines((k2$eval.points),k2$estimate,col=marker[[2]],lwd=2)
lines((k3$eval.points),k3$estimate,col=marker[[4]],lwd=2)
lines((k4$eval.points),k4$estimate,col=marker[[6]],lwd=2)
lines((k5$eval.points),k5$estimate,col=marker[[7]],lwd=2)
dev.off()


twoVillage32000_TBP_migRate_0.1 = read.csv(paste0(dir4,"time_between_pcases_3200032000reintRate_0.001000migRate_0.100000_paper.csv"),sep=",",header = FALSE,col.names = paste0("V",seq_len(5000)), fill = TRUE)
twoVillage32000_extInt_migRate_0.1 = read.table(paste0(dir4,"extinction_interval_3200032000reintRate_0.001000migRate_0.100000_paper.csv"))
fourVillage16000_TBP_migRate_0.1 = read.csv(paste0(dir4,"time_between_pcases_16000160001600016000reintRate_0.001000migRate_0.100000_paper.csv"),sep=",",header = FALSE,col.names = paste0("V",seq_len(2000)), fill = TRUE)
fourVillage16000_extInt_migRate_0.1 = read.table(paste0(dir4,"extinction_interval_16000160001600016000reintRate_0.001000migRate_0.100000_paper.csv"))
eightVillage8000_TBP_migRate_0.1 = read.csv(paste0(dir4,"time_between_pcases_80008000800080008000800080008000reintRate_0.001000migRate_0.100000_paper.csv"),sep=",",header = FALSE,col.names = paste0("V",seq_len(1000)), fill = TRUE)
eightVillage8000_extInt_migRate_0.1 = read.table(paste0(dir4,"extinction_interval_80008000800080008000800080008000reintRate_0.001000migRate_0.100000_paper.csv"))
sixteenVillage4000_TBP_migRate_0.1 = read.csv(paste0(dir4,"time_between_pcases_4000400040004000400040004000400040004000400040004000400040004000reintRate_0.001000migRate_0.100000_paper.csv"),sep=",",header = FALSE,col.names = paste0("V",seq_len(1000)), fill = TRUE)
sixteenVillage4000_extInt_migRate_0.1 = read.table(paste0(dir4,"extinction_interval_4000400040004000400040004000400040004000400040004000400040004000reintRate_0.001000migRate_0.100000_paper.csv"))
thirtytwoVillage2000_TBP_migRate_0.1 = read.csv(paste0(dir4,"time_between_pcases_20002000200020002000200020002000200020002000200020002000200020002000200020002000200020002000200020002000200020002000200020002000reintRate_0.001000migRate_0.100000_paper.csv"),sep=",",header = FALSE,col.names = paste0("V",seq_len(500)), fill = TRUE)
thirtytwoVillage2000_extInt_migRate_0.1 = read.csv(paste0(dir4,"extinction_interval_20002000200020002000200020002000200020002000200020002000200020002000200020002000200020002000200020002000200020002000200020002000reintRate_0.001000migRate_0.100000_paper.csv"))


#remove nas
new <- unlist(oneVillage64000_extInt)[!is.na(unlist(oneVillage64000_extInt))]
new1<-unlist(twoVillage32000_extInt_migRate_0.1)[!is.na(unlist(twoVillage32000_extInt_migRate_0.1))]
new2<-unlist(fourVillage16000_extInt_migRate_0.1)[!is.na(unlist(fourVillage16000_extInt_migRate_0.1))]
new3<-unlist(eightVillage8000_extInt_migRate_0.1)[!is.na(unlist(eightVillage8000_extInt_migRate_0.1))]
new4<-unlist(sixteenVillage4000_extInt_migRate_0.1)[!is.na(unlist(sixteenVillage4000_extInt_migRate_0.1))]
new5<-unlist(thirtytwoVillage2000_extInt_migRate_0.1)[!is.na(unlist(thirtytwoVillage2000_extInt_migRate_0.1))]

#density estimates of log of data
k = kde(x=as.matrix(new),positive=TRUE)
k1 = kde(x=as.matrix((new1)),positive=TRUE)
k2 = kde(x=as.matrix((new2)),positive=TRUE)
k3 = kde(x=as.matrix((new3)),positive=TRUE)
k4 = kde(x=as.matrix((new4)),positive=TRUE)
k5 = kde(x=as.matrix((new5)),positive=TRUE)

#plot smoothed curve on top of hist
hist(new,freq=FALSE)
lines(k$eval.points,k$estimate,col='red')
hist(new1,freq=FALSE)
lines(k1$eval.points,k1$estimate,col='red')
hist(new2,freq=FALSE)
lines(k2$eval.points,k2$estimate,col='red')
hist(new3,freq=FALSE)
lines(k3$eval.points,k3$estimate,col='red')
hist(new4,freq=FALSE)
lines(k4$eval.points,k4$estimate,col='red')
hist(new5,freq=FALSE)
lines(k5$eval.points,k5$estimate,col='red')


png(paste0(dir2,"extInt_dist_migRate0point1_multipatch_paper_v1.png"), width=1400, height=800, res=150) 
par(fig=c(0,1,0,1),mar=c(3.3,3.6,1.1,1.1),oma=c(1,0,0,3),mgp=c(2.2,0.45,0),tcl=-0.4,cex=1.2)
plot((k$eval.points),k$estimate,col='black',type='l',ylim=c(0,1),xlim=c(0,5),xlab='',ylab='Density',lwd=2)
lines((k1$eval.points),k1$estimate,col=marker[[1]],lwd=2)
lines((k2$eval.points),k2$estimate,col=marker[[2]],lwd=2)
lines((k3$eval.points),k3$estimate,col=marker[[4]],lwd=2)
lines((k4$eval.points),k4$estimate,col=marker[[6]],lwd=2)
lines((k5$eval.points),k5$estimate,col=marker[[7]],lwd=2)
#mtext("Density",side=2,outer=TRUE,cex=1.2)
mtext('Length of extinction interval (Years)',side=1,outer=TRUE,cex=1.2)
legend(x=1.5,y=.9,legend=rev(c('1x64k','2x32k','4x16k','8x8k','16x4k','32x2k')),col=rev(c('black',marker[1:2],marker[4],marker[6:7])),lwd=2,bty='n')
par(new=TRUE)
par(mar=c(0,0,0,0))
par(fig=c(0.6,0.95,0.5,0.85))
plot((k$eval.points),k$estimate,col='black',type='l',ylim=c(0,0.02),xlim=c(3,5),xlab='',ylab='',lwd=2)
lines((k1$eval.points),k1$estimate,col=marker[[1]],lwd=2)
lines((k2$eval.points),k2$estimate,col=marker[[2]],lwd=2)
lines((k3$eval.points),k3$estimate,col=marker[[4]],lwd=2)
lines((k4$eval.points),k4$estimate,col=marker[[6]],lwd=2)
lines((k5$eval.points),k5$estimate,col=marker[[7]],lwd=2)
dev.off()


new <- unlist(oneVillage64000_TBP)[!is.na(unlist(oneVillage64000_TBP))]
new1<-unlist(twoVillage32000_TBP_migRate_0.1)[!is.na(unlist(twoVillage32000_TBP_migRate_0.1))]
new2<-unlist(fourVillage16000_TBP_migRate_0.1)[!is.na(unlist(fourVillage16000_TBP_migRate_0.1))]
new3<-unlist(eightVillage8000_TBP_migRate_0.1)[!is.na(unlist(eightVillage8000_TBP_migRate_0.1))]
new4<-unlist(sixteenVillage4000_TBP_migRate_0.1)[!is.na(unlist(sixteenVillage4000_TBP_migRate_0.1))]
new5<-unlist(thirtytwoVillage2000_TBP_migRate_0.1)[!is.na(unlist(thirtytwoVillage2000_TBP_migRate_0.1))]

#density estimates of log of data
k = kde(x=as.matrix(new),positive=TRUE)
k1 = kde(x=as.matrix((new1)),positive=TRUE)
k2 = kde(x=as.matrix((new2)),positive=TRUE)
k3 = kde(x=as.matrix((new3)),positive=TRUE)
k4 = kde(x=as.matrix((new4)),positive=TRUE)
k5 = kde(x=as.matrix((new5)),positive=TRUE)


#plot smoothed curve on top of hist
hist(new,freq=FALSE)
lines(k$eval.points,k$estimate,col='red')
hist(new1,freq=FALSE)
lines(k1$eval.points,k1$estimate,col='red')
hist(new2,freq=FALSE)
lines(k2$eval.points,k2$estimate,col='red')
hist(new3,freq=FALSE)
lines(k3$eval.points,k3$estimate,col='red')
hist(new4,freq=FALSE)
lines(k4$eval.points,k4$estimate,col='red')
hist(new5,freq=FALSE)
lines(k5$eval.points,k5$estimate,col='red')

png(paste0(dir2,"interInt_dist_migRate0point1_multipatch_paper_v1.png"), width=1400, height=800, res=150) 
par(fig=c(0,1,0,1),mar=c(3.3,3.6,1.1,1.1),oma=c(1,0,0,3),mgp=c(2.2,0.45,0),tcl=-0.4,cex=1.2)
plot((k$eval.points),k$estimate,col='black',type='l',ylim=c(0,1),xlim=c(-.3,5),xlab='',ylab='Density',lwd=2)
lines((k1$eval.points),k1$estimate,col=marker[[1]],lwd=2)
lines((k2$eval.points),k2$estimate,col=marker[[2]],lwd=2)
lines((k3$eval.points),k3$estimate,col=marker[[4]],lwd=2)
lines((k4$eval.points),k4$estimate,col=marker[[6]],lwd=2)
lines((k5$eval.points),k5$estimate,col=marker[[7]],lwd=2)
#mtext("Density",side=2,outer=TRUE,cex=1.2)
mtext('Length of intercase interval (Years)',side=1,outer=TRUE,cex=1.2)
legend(x=1,y=0.8,legend=rev(c('1x64k','2x32k','4x16k','8x8k','16x4k','32x2k')),col=rev(c('black',marker[1:2],marker[4],marker[6:7])),lwd=2,bty='n')
par(new=TRUE)
par(mar=c(0,0,0,0))
par(fig=c(0.6,0.95,0.5,0.85))
plot((k$eval.points),k$estimate,col='black',type='l',ylim=c(0,0.0015),xlim=c(3,5),xlab='',ylab='',lwd=2)
lines((k1$eval.points),k1$estimate,col=marker[[1]],lwd=2)
lines((k2$eval.points),k2$estimate,col=marker[[2]],lwd=2)
lines((k3$eval.points),k3$estimate,col=marker[[4]],lwd=2)
lines((k4$eval.points),k4$estimate,col=marker[[6]],lwd=2)
lines((k5$eval.points),k5$estimate,col=marker[[7]],lwd=2)
dev.off()



##single pop
dir4 = '/Users/Celeste/Desktop/multipatch_model/metapopulation_polio/'


oneVillage64000_TBP_v2 = read.csv(paste0(dir4,"time_between_pcases_64000reintRate_0.001000migRate_0.000000_paper.csv"),sep=",",header = FALSE,col.names = paste0("V",seq_len(9000)), fill = TRUE)
oneVillage64000_extInt_v2 = read.table(paste0(dir4,"extinction_interval_64000reintRate_0.001000migRate_0.000000_paper.csv"))
oneVillage32000_TBP_v2 = read.csv(paste0(dir4,"time_between_pcases_32000reintRate_0.001000migRate_0.000000_paper.csv"),sep=",",header = FALSE,col.names = paste0("V",seq_len(1000)), fill = TRUE)
oneVillage32000_extInt_v2 = read.table(paste0(dir4,"extinction_interval_32000reintRate_0.001000migRate_0.000000_paper.csv"))
oneVillage16000_TBP_v2 = read.csv(paste0(dir4,"time_between_pcases_16000reintRate_0.001000migRate_0.000000_paper.csv"),sep=",",header = FALSE,col.names = paste0("V",seq_len(180)), fill = TRUE)
oneVillage16000_extInt_v2 = read.table(paste0(dir4,"extinction_interval_16000reintRate_0.001000migRate_0.000000_paper.csv"))
oneVillage8000_TBP_v2 = read.csv(paste0(dir4,"time_between_pcases_8000reintRate_0.001000migRate_0.000000_paper.csv"),sep=",",header = FALSE,col.names = paste0("V",seq_len(180)), fill = TRUE)
oneVillage8000_extInt_v2 = read.table(paste0(dir4,"extinction_interval_8000reintRate_0.001000migRate_0.000000_paper.csv"))
oneVillage4000_TBP_v2 = read.csv(paste0(dir4,"time_between_pcases_4000reintRate_0.001000migRate_0.000000_paper.csv"),sep=",",header = FALSE,col.names = paste0("V",seq_len(180)), fill = TRUE)
oneVillage4000_extInt_v2 = read.table(paste0(dir4,"extinction_interval_4000reintRate_0.001000migRate_0.000000_paper.csv"))
oneVillage2000_TBP_v2 = read.csv(paste0(dir4,"time_between_pcases_2000reintRate_0.001000migRate_0.000000_paper.csv"),sep=",",header = FALSE,col.names = paste0("V",seq_len(180)), fill = TRUE)
oneVillage2000_extInt_v2 = read.table(paste0(dir4,"extinction_interval_2000reintRate_0.001000migRate_0.000000_paper.csv"))

new <- unlist(oneVillage64000_extInt_v2)[!is.na(unlist(oneVillage64000_extInt_v2))]
new1<-unlist(oneVillage32000_extInt_v2)[!is.na(unlist(oneVillage32000_extInt_v2))]
new2<-unlist(oneVillage16000_extInt_v2)[!is.na(unlist(oneVillage16000_extInt_v2))]
new3<-unlist(oneVillage8000_extInt_v2)[!is.na(unlist(oneVillage8000_extInt_v2))]
new4<-unlist(oneVillage4000_extInt_v2)[!is.na(unlist(oneVillage4000_extInt_v2))]
new5<-unlist(oneVillage2000_extInt_v2)[!is.na(unlist(oneVillage2000_extInt_v2))]

#density estimates of log of data
k = kde(x=as.matrix(new),positive=TRUE)
k1 = kde(x=as.matrix((new1)),positive=TRUE)
k2 = kde(x=as.matrix((new2)),positive=TRUE)
k3 = kde(x=as.matrix((new3)),positive=TRUE)
k4 = kde(x=as.matrix((new4)),positive=TRUE)
k5 = kde(x=as.matrix((new5)),positive=TRUE)

#example of fit
png(paste0(dir2,"hist_gauss_curve_fit_64000.png"), width=1400, height=800, res=150) 
hist(new,freq=FALSE,xlab="Extinction intervals (years)",main="")
lines(k$eval.points,k$estimate,col='red')
dev.off()


png(paste0(dir2,"extInt_dist_singlePop_paper_v1.png"), width=1400, height=800, res=150) 
par(fig=c(0,1,0,1),mar=c(3.3,3.6,1.1,1.1),oma=c(1,0,0,3),mgp=c(2.2,0.45,0),tcl=-0.4,cex=1.2)
plot((k$eval.points),k$estimate,col='black',type='l',ylim=c(0,2.1),xlim=c(0,5),xlab='',ylab='Density',lwd=2)
lines((k1$eval.points),k1$estimate,col=marker[[1]],lwd=2)
lines((k2$eval.points),k2$estimate,col=marker[[2]],lwd=2)
lines((k3$eval.points),k3$estimate,col=marker[[4]],lwd=2)
lines((k4$eval.points),k4$estimate,col=marker[[6]],lwd=2)
lines((k5$eval.points),k5$estimate,col=marker[[7]],lwd=2)
#mtext("Density",side=2,outer=TRUE,cex=1.2)
mtext('Length of extinction interval (Years)',side=1,outer=TRUE,cex=1.2)
legend(x=1,y=2,legend=rev(c('1x64k','1x32k','1x16k','1x8k','1x4k','1x2k')),col=rev(c('black',marker[1:2],marker[4],marker[6:7])),lwd=2,bty='n')
par(new=TRUE)
par(mar=c(0,0,0,0))
par(fig=c(0.6,0.95,0.5,0.85))
plot((k$eval.points),k$estimate,col='black',type='l',ylim=c(0,0.02),xlim=c(3,5),xlab='',ylab='',lwd=2)
lines((k1$eval.points),k1$estimate,col=marker[[1]],lwd=2)
lines((k2$eval.points),k2$estimate,col=marker[[2]],lwd=2)
lines((k3$eval.points),k3$estimate,col=marker[[4]],lwd=2)
lines((k4$eval.points),k4$estimate,col=marker[[6]],lwd=2)
lines((k5$eval.points),k5$estimate,col=marker[[7]],lwd=2)
dev.off()

new <- unlist(oneVillage64000_TBP_v2)[!is.na(unlist(oneVillage64000_TBP_v2))]
new1<-unlist(oneVillage32000_TBP_v2)[!is.na(unlist(oneVillage32000_TBP_v2))]
new2<-unlist(oneVillage16000_TBP_v2)[!is.na(unlist(oneVillage16000_TBP_v2))]
new3<-unlist(oneVillage8000_TBP_v2)[!is.na(unlist(oneVillage8000_TBP_v2))]
new4<-unlist(oneVillage4000_TBP_v2)[!is.na(unlist(oneVillage4000_TBP_v2))]
new5<-unlist(oneVillage2000_TBP_v2)[!is.na(unlist(oneVillage2000_TBP_v2))]
#density estimates of log of data
k = kde(x=as.matrix(new),positive=TRUE)
k1 = kde(x=as.matrix((new1)),positive=TRUE)
k2 = kde(x=as.matrix((new2)),positive=TRUE)
k3 = kde(x=as.matrix((new3)),positive=TRUE)
k4 = kde(x=as.matrix((new4)),positive=TRUE)
k5 = kde(x=as.matrix((new5)),positive=TRUE)

png(paste0(dir2,"interInt_dist_singlePop_paper_v1.png"), width=1400, height=800, res=150) 
par(fig=c(0,1,0,1),mar=c(3.3,3.6,1.1,1.1),oma=c(1,0,0,3),mgp=c(2.2,0.45,0),tcl=-0.4,cex=1.2)
plot((k$eval.points),k$estimate,col='black',type='l',ylim=c(0,1),xlim=c(-.3,5),xlab='',ylab='Density',lwd=2)
lines((k1$eval.points),k1$estimate,col=marker[[1]],lwd=2)
lines((k2$eval.points),k2$estimate,col=marker[[2]],lwd=2)
lines((k3$eval.points),k3$estimate,col=marker[[4]],lwd=2)
lines((k4$eval.points),k4$estimate,col=marker[[6]],lwd=2)
lines((k5$eval.points),k5$estimate,col=marker[[7]],lwd=2)
#mtext("Density",side=2,outer=TRUE,cex=1.2)
mtext('Length of intercase interval (Years)',side=1,outer=TRUE,cex=1.2)
legend(x=1,y=0.8,legend=rev(c('1x64k','1x32k','1x16k','1x8k','1x4k','1x2k')),col=rev(c('black',marker[1:2],marker[4],marker[6:7])),lwd=2,bty='n')
par(new=TRUE)
par(mar=c(0,0,0,0))
par(fig=c(0.6,0.95,0.5,0.85))
plot((k$eval.points),k$estimate,col='black',type='l',ylim=c(0,0.015),xlim=c(3,5),xlab='',ylab='',lwd=2)
lines((k1$eval.points),k1$estimate,col=marker[[1]],lwd=2)
lines((k2$eval.points),k2$estimate,col=marker[[2]],lwd=2)
lines((k3$eval.points),k3$estimate,col=marker[[4]],lwd=2)
lines((k4$eval.points),k4$estimate,col=marker[[6]],lwd=2)
lines((k5$eval.points),k5$estimate,col=marker[[7]],lwd=2)
dev.off()

oneVillage64000_TBP_v2 = read.csv(paste0(dir4,"time_between_pcases_64000reintRate_0.001000migRate_0.000000_paper.csv"),sep=",",header = FALSE,col.names = paste0("V",seq_len(9000)), fill = TRUE)
oneVillage64000_extInt_v2 = read.table(paste0(dir4,"extinction_interval_64000reintRate_0.001000migRate_0.000000_paper.csv"))
oneVillage32000_TBP_v2_mig01 = read.csv(paste0(dir4,"time_between_pcases_3200032000reintRate_0.001000migRate_0.100000_paper.csv"),sep=",",header = FALSE,col.names = paste0("V",seq_len(5000)), fill = TRUE)
oneVillage32000_extInt_v2_mig01 = read.table(paste0(dir4,"extinction_interval_3200032000reintRate_0.001000migRate_0.100000_paper.csv"))
oneVillage16000_TBP_v2_mig01 = read.csv(paste0(dir4,"time_between_pcases_16000160001600016000reintRate_0.001000migRate_0.100000_paper.csv"),sep=",",header = FALSE,col.names = paste0("V",seq_len(5000)), fill = TRUE)
oneVillage16000_extInt_v2_mig01 = read.table(paste0(dir4,"extinction_interval_16000160001600016000reintRate_0.001000migRate_0.100000_paper.csv"))
oneVillage8000_TBP_v2_mig01 = read.csv(paste0(dir4,"time_between_pcases_80008000800080008000800080008000reintRate_0.001000migRate_0.100000_paper.csv"),sep=",",header = FALSE,col.names = paste0("V",seq_len(2000)), fill = TRUE)
oneVillage8000_extInt_v2_mig01 = read.table(paste0(dir4,"extinction_interval_80008000800080008000800080008000reintRate_0.001000migRate_0.100000_paper.csv"))
oneVillage4000_TBP_v2_mig01 = read.csv(paste0(dir4,"time_between_pcases_4000400040004000400040004000400040004000400040004000400040004000reintRate_0.001000migRate_0.100000_paper.csv"),sep=",",header = FALSE,col.names = paste0("V",seq_len(180)), fill = TRUE)
oneVillage4000_extInt_v2_mig01 = read.table(paste0(dir4,"extinction_interval_4000400040004000400040004000400040004000400040004000400040004000reintRate_0.001000migRate_0.100000_paper.csv"))
oneVillage2000_TBP_v2_mig01 = read.csv(paste0(dir4,"time_between_pcases_20002000200020002000200020002000200020002000200020002000200020002000200020002000200020002000200020002000200020002000200020002000reintRate_0.001000migRate_0.100000_paper.csv"),sep=",",header = FALSE,col.names = paste0("V",seq_len(180)), fill = TRUE)
oneVillage2000_extInt_v2_mig01 = read.table(paste0(dir4,"extinction_interval_20002000200020002000200020002000200020002000200020002000200020002000200020002000200020002000200020002000200020002000200020002000reintRate_0.001000migRate_0.100000_paper.csv"))

new <- unlist(oneVillage64000_extInt_v2)[!is.na(unlist(oneVillage64000_extInt_v2))]
new1<-unlist(oneVillage32000_extInt_v2_mig01)[!is.na(unlist(oneVillage32000_extInt_v2_mig01))]
new2<-unlist(oneVillage16000_extInt_v2_mig01)[!is.na(unlist(oneVillage16000_extInt_v2_mig01))]
new3<-unlist(oneVillage8000_extInt_v2_mig01)[!is.na(unlist(oneVillage8000_extInt_v2_mig01))]
new4<-unlist(oneVillage4000_extInt_v2_mig01)[!is.na(unlist(oneVillage4000_extInt_v2_mig01))]
new5<-unlist(oneVillage2000_extInt_v2_mig01)[!is.na(unlist(oneVillage2000_extInt_v2_mig01))]

#density estimates of log of data
k = kde(x=as.matrix(new),positive=TRUE)
k1 = kde(x=as.matrix((new1)),positive=TRUE)
k2 = kde(x=as.matrix((new2)),positive=TRUE)
k3 = kde(x=as.matrix((new3)),positive=TRUE)
k4 = kde(x=as.matrix((new4)),positive=TRUE)
k5 = kde(x=as.matrix((new5)),positive=TRUE)

#example of fit
png(paste0(dir2,"hist_gauss_curve_fit_64000.png"), width=1400, height=800, res=150) 
hist(new,freq=FALSE,xlab="Extinction intervals (years)",main="")
lines(k$eval.points,k$estimate,col='red')
dev.off()


png(paste0(dir2,"extInt_dist_multi_pop_mig01_paper.png"), width=1400, height=800, res=150) 
par(mfrow=c(1,2),oma=c(1,1,1,1),cex=1.2,mar=c(2.5,2.5,1,1))
plot(k$eval.points,k$estimate,col='black',type='l',ylim=c(0,1),xlim=c(0,5),xlab='',ylab='',lwd=2)
lines((k1$eval.points),k1$estimate,col=marker$color[4],lwd=2)
lines((k2$eval.points),k2$estimate,col=marker$color[5],lwd=2)
lines((k3$eval.points),k3$estimate,col=marker$color[6],lwd=2)
lines((k4$eval.points),k4$estimate,col=marker$color[7],lwd=2)
lines((k5$eval.points),k5$estimate,col=marker$color[8],lwd=2)
mtext("Density",side=2,outer=TRUE,cex=1.2)
legend('topright',legend=c('1x64000','2x32000','4x16000','8x8000','16x4000','32x2000'),col=c('black',marker$color[4:8]),lwd=2,bty='n')
plot(k$eval.points,k$estimate,col='black',type='l',ylim=c(0,0.018),xlim=c(3,5),xlab='',ylab='',lwd=2)
lines((k1$eval.points),k1$estimate,col=marker$color[4],lwd=2)
lines((k2$eval.points),k2$estimate,col=marker$color[5],lwd=2)
lines((k3$eval.points),k3$estimate,col=marker$color[6],lwd=2)
lines((k4$eval.points),k4$estimate,col=marker$color[7],lwd=2)
lines((k5$eval.points),k5$estimate,col=marker$color[8],lwd=2)
mtext('Length of extinction interval (years)',side=1,outer=TRUE,cex=1.2)
dev.off()

new <- unlist(oneVillage64000_TBP_v2)[!is.na(unlist(oneVillage64000_TBP_v2))]
new1<-unlist(oneVillage32000_TBP_v2_mig01)[!is.na(unlist(oneVillage32000_TBP_v2_mig01))]
new2<-unlist(oneVillage16000_TBP_v2_mig01)[!is.na(unlist(oneVillage16000_TBP_v2_mig01))]
new3<-unlist(oneVillage8000_TBP_v2_mig01)[!is.na(unlist(oneVillage8000_TBP_v2_mig01))]
new4<-unlist(oneVillage4000_TBP_v2_mig01)[!is.na(unlist(oneVillage4000_TBP_v2_mig01))]
new5<-unlist(oneVillage2000_TBP_v2_mig01)[!is.na(unlist(oneVillage2000_TBP_v2_mig01))]
#density estimates of log of data
k = kde(x=as.matrix(new),positive=TRUE)
k1 = kde(x=as.matrix((new1)),positive=TRUE)
k2 = kde(x=as.matrix((new2)),positive=TRUE)
k3 = kde(x=as.matrix((new3)),positive=TRUE)
k4 = kde(x=as.matrix((new4)),positive=TRUE)
k5 = kde(x=as.matrix((new5)),positive=TRUE)

png(paste0(dir2,"interInt_dist_multi_pop_mig01_paper.png"), width=1400, height=800, res=150) 
par(mfrow=c(1,2),oma=c(1,1,1,1),cex=1.2,mar=c(2.5,2.5,1,1))
plot(k$eval.points,k$estimate,col='black',type='l',ylim=c(0,1),xlim=c(0,5),xlab='',ylab='',lwd=2)
lines((k1$eval.points),k1$estimate,col=marker$color[4],lwd=2)
lines((k2$eval.points),k2$estimate,col=marker$color[5],lwd=2)
lines((k3$eval.points),k3$estimate,col=marker$color[6],lwd=2)
lines((k4$eval.points),k4$estimate,col=marker$color[7],lwd=2)
lines((k5$eval.points),k5$estimate,col=marker$color[8],lwd=2)
mtext("Density",side=2,outer=TRUE,cex=1.2)
legend('topright',legend=c('1x64000','2x32000','4x16000','8x8000','16x4000','32x2000'),col=c('black',marker$color[4:8]),lwd=2,bty='n')
plot(k$eval.points,k$estimate,col='black',type='l',ylim=c(0,0.0015),xlim=c(3,5),xlab='',ylab='',lwd=2)
lines((k1$eval.points),k1$estimate,col=marker$color[4],lwd=2)
lines((k2$eval.points),k2$estimate,col=marker$color[5],lwd=2)
lines((k3$eval.points),k3$estimate,col=marker$color[6],lwd=2)
lines((k4$eval.points),k4$estimate,col=marker$color[7],lwd=2)
lines((k5$eval.points),k5$estimate,col=marker$color[8],lwd=2)
mtext('Length of intercase interval (years)',side=1,outer=TRUE,cex=1.2)
dev.off()












         
         
