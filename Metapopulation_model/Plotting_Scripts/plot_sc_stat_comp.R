# plots affect of population partitioning and movement rate on sc statistic

dir = '/Users/Celeste/Desktop/multipatch_model/SC_statistic/'
dir1 = '/Users/Celeste/Desktop/multipatch_model/'
dir2 = '/Metapopulation_model/Figures/'

library(RColorBrewer)


#plot comparison of 0pcase (assuming case at start of sim) vs 1pcase (waiting for simulated case)
sixteenVillage_4000_0pcase = read.table(paste0(dir,'0pcase_16x4000_migRate_0_paper_092120.csv'), col.names=c('time', 'E&D'))
sixteenVillage_4000_1pcase = read.table(paste0(dir,'1pcase_16x4000_migRate_0_paper.csv'), col.names=c('time', 'E&D'))


png(paste0(dir2,"SC_comp_0pcase_vs_1pcase_16x4000_migRate_0_092120.png"), width=1400, height=800, res=150)
plot(sixteenVillage_4000_1pcase,type='l',xlab="Years since last detected paralytic case",ylab="Probability of silent circulation")
lines(sixteenVillage_4000_0pcase,col='red')
legend('topright',legend=c('NICA','ICA'),col=c('black','red'),lwd=2,bty='n')
dev.off()


# single pop
oneVillage_64000_seasonality1yrperiod = read.table(paste0(dir,'1pcase_1x64000_migRate_0_paper.csv'), col.names=c('time', 'E&D'))
oneVillage_16000_seasonality1yrperiod = read.table(paste0(dir,'1pcase_1x16000_migRate_0_paper.csv'), col.names=c('time', 'E&D'))
oneVillage_32000_seasonality1yrperiod = read.table(paste0(dir,'1pcase_1x32000_migRate_0_paper.csv'), col.names=c('time', 'E&D'))
oneVillage_8000_seasonality1yrperiod = read.table(paste0(dir,'1pcase_1x8000_migRate_0_paper.csv'), col.names=c('time', 'E&D'))
oneVillage_4000_seasonality1yrperiod = read.table(paste0(dir,'1pcase_1x4000_migRate_0_paper.csv'), col.names=c('time', 'E&D'))
oneVillage_2000_seasonality1yrperiod = read.table(paste0(dir,'1pcase_1x2000_migRate_0_paper.csv'), col.names=c('time', 'E&D'))

a0 = approx(oneVillage_64000_seasonality1yrperiod$time,oneVillage_64000_seasonality1yrperiod$E.D,xout=seq(0,10,0.01))
a = approx(oneVillage_16000_seasonality1yrperiod$time,oneVillage_16000_seasonality1yrperiod$E.D,xout=seq(0,10,0.01))
a1 = approx(oneVillage_32000_seasonality1yrperiod$time,oneVillage_32000_seasonality1yrperiod$E.D,xout=seq(0,10,0.01))
a2 = approx(oneVillage_8000_seasonality1yrperiod$time,oneVillage_8000_seasonality1yrperiod$E.D,xout=seq(0,10,0.01))
a3 = approx(oneVillage_4000_seasonality1yrperiod$time,oneVillage_4000_seasonality1yrperiod$E.D,xout=seq(0,10,0.01))
a4 = approx(oneVillage_2000_seasonality1yrperiod$time,oneVillage_2000_seasonality1yrperiod$E.D,xout=seq(0,10,0.01))


png(paste0(dir2,"SC_comp_village_subdivision_1x64000_single_pop_zoom.png"), width=1400, height=800, res=150) 
par(mfrow=c(1,2),oma=c(1,1,1,1),cex=1.2,mar=c(2.5,2.5,1,1))
plot(a0$x,(a0$y-a1$y),type='l',col=marker$color[4],xlab="",ylab="",xlim=c(0,5),lwd=2)
lines(a0$x,(a0$y-a$y),col=marker$color[5],lwd=2)
lines(a0$x,(a0$y-a2$y),col=marker$color[6],lwd=2)
lines(a0$x,(a0$y-a3$y),col=marker$color[7],lwd=2)
lines(a0$x,(a0$y-a4$y),col=marker$color[8],lwd=2)
legend('topright',legend=c('1x64000','1x32000','1x16000','1x8000','1x4000','1x2000'),col=marker$color[4:8],lwd=2,bty='n')
mtext("Probability of silent circulation differential",side=2,outer=TRUE,cex=1.2)
#abline(v=3,col='blue',lty=3)
#legend('topright',legend=c('1 village (1x64000)','2 villages (2x32000)','4 villages (4x16000)','8 villages (8x8000)','16 villages (16x4000)','5 villages (1x32000,4x8000)'),col=c('black',marker$color[4:7],'darkgreen'),lwd=c(2,2,2,2,2,4),bty='n')
plot(a0$x,(a0$y-a1$y),type='l',col=marker$color[4],xlab="",ylab="",xlim=c(3,5),lwd=2)
lines(a0$x,(a0$y-a$y),col=marker$color[5],lwd=2)
lines(a0$x,(a0$y-a2$y),col=marker$color[6],lwd=2)
lines(a0$x,(a0$y-a3$y),col=marker$color[7],lwd=2)
lines(a0$x,(a0$y-a4$y),col=marker$color[8],lwd=2)
mtext("Time since last detected paralytic case (years)",side=1,outer=TRUE,cex=1.2)
dev.off()


#compare single pop to multipop by adding patches
fourVillage_16000_seasonality1yrperiod = read.table(paste0(dir,'1pcase_4x16000_migRate_0_paper.csv'), col.names=c('time', 'E&D'))
twoVillage_32000_seasonality1yrperiod = read.table(paste0(dir,'1pcase_2x32000_migRate_0_paper.csv'), col.names=c('time', 'E&D'))
eightVillage_8000_seasonality1yrperiod = read.table(paste0(dir,'1pcase_8x8000_migRate_0_paper.csv'), col.names=c('time', 'E&D'))
sixteenVillage_4000_seasonality1yrperiod = read.table(paste0(dir,'1pcase_16x4000_migRate_0_paper.csv'), col.names=c('time', 'E&D'))
thirtytwoVillage_2000_seasonality1yrperiod = read.table(paste0(dir,'1pcase_32x2000_migRate_0_paper.csv'), col.names=c('time', 'E&D'))

#first approx functions for equal comparison
a0 = approx(oneVillage_64000_seasonality1yrperiod$time,oneVillage_64000_seasonality1yrperiod$E.D,xout=seq(0,10,0.01))
aa = approx(twoVillage_32000_seasonality1yrperiod$time,twoVillage_32000_seasonality1yrperiod$E.D,xout=seq(0,10,0.01))
aa1 = approx(fourVillage_16000_seasonality1yrperiod$time,fourVillage_16000_seasonality1yrperiod$E.D,xout=seq(0,10,0.01))
aa2 = approx(eightVillage_8000_seasonality1yrperiod$time,eightVillage_8000_seasonality1yrperiod$E.D,xout=seq(0,10,0.01))
aa3 = approx(sixteenVillage_4000_seasonality1yrperiod$time,sixteenVillage_4000_seasonality1yrperiod$E.D,xout=seq(0,10,0.01))
aa4 = approx(thirtytwoVillage_2000_seasonality1yrperiod$time,thirtytwoVillage_2000_seasonality1yrperiod$E.D,xout=seq(0,10,0.01))

a0$y[is.na(a0$y)]<-0
aa$y[is.na(aa$y)]<-0
aa1$y[is.na(aa1$y)]<-0
aa2$y[is.na(aa2$y)]<-0
aa3$y[is.na(aa3$y)]<-0
aa4$y[is.na(aa4$y)]<-0

marker=colorRampPalette(c('orange','red','purple','royalblue'))(7)


png(paste0(dir2,"psc_differential_metapop_paper_v1.png"), width=1400, height=800, res=150) 
par(fig=c(0,1,0,1),mar=c(3.3,3.6,1.1,1.1),oma=c(1,1,0,3),mgp=c(2.2,0.45,0),tcl=-0.4,cex=1.2)
plot(a0$x,(a0$y-aa$y),type='l',col=marker[[1]],lwd=2,ylim=c(-.15,.8),xlim=c(-.05,5),xlab="",ylab="Probability of silent circualation differential")
lines(a0$x,(a0$y-aa1$y),col=marker[[2]],lwd=2)
lines(a0$x,(a0$y-aa2$y),col=marker[[4]],lwd=2)
lines(a0$x,(a0$y-aa3$y),col=marker[[6]],lwd=2)
lines(a0$x,(a0$y-aa4$y),col=marker[[7]],lwd=2)
legend('topleft',legend=rev(c('2x32k','4x16k','8x8k','16x4k','32x2k')),col=rev(c(marker[1:2],marker[4],marker[6:7])),lwd=2,bty='n')
mtext(side=1,'Years since paralytic case',outer=TRUE,cex=1.2)
abline(h=0,lty=3)
mtext("1x64k\n has\n higher\n probability",side=4,at=.6,las=2,line=-1,outer=TRUE)
mtext("1x64k\n has\n lower\n probability",side=4,at=.22,las=2,line=-1,outer=TRUE)
par(new=TRUE)
par(mar=c(0,0,0,0))
par(fig=c(0.58,0.92,0.5,0.85))
plot(a0$x,(a0$y-aa$y),type='l',lwd=2,col=marker[[1]],ylim=c(-.06,.01),xlim=c(3,5),xlab="",ylab="")
lines(a0$x,(a0$y-aa1$y),col=marker[[2]],lwd=2)
lines(a0$x,(a0$y-aa2$y),col=marker[[4]],lwd=2)
lines(a0$x,(a0$y-aa3$y),col=marker[[6]],lwd=2)
lines(a0$x,(a0$y-aa4$y),col=marker[[7]],lwd=2)
abline(h=0,lty=3)
dev.off()

#add migration
fourVillage_16000_seasonality1yrperiod_mig01 = read.table(paste0(dir,'1pcase_4x16000_migRate_0.1_paper.csv'), col.names=c('time', 'E&D'))
twoVillage_32000_seasonality1yrperiod_mig01 = read.table(paste0(dir,'1pcase_2x32000_migRate_0.1_paper.csv'), col.names=c('time', 'E&D'))
eightVillage_8000_seasonality1yrperiod_mig01 = read.table(paste0(dir,'1pcase_8x8000_migRate_0.1_paper.csv'), col.names=c('time', 'E&D'))
sixteenVillage_4000_seasonality1yrperiod_mig01 = read.table(paste0(dir,'1pcase_16x4000_migRate_0.1_paper.csv'), col.names=c('time', 'E&D'))
thirtytwoVillage_2000_seasonality1yrperiod_mig01 = read.table(paste0(dir,'1pcase_32x2000_migRate_0.1_paper.csv'), col.names=c('time', 'E&D'))

#first approx functions for equal comparison
a0 = approx(oneVillage_64000_seasonality1yrperiod$time,oneVillage_64000_seasonality1yrperiod$E.D,xout=seq(0,10,0.01))
aa = approx(twoVillage_32000_seasonality1yrperiod_mig01$time,twoVillage_32000_seasonality1yrperiod_mig01$E.D,xout=seq(0,10,0.01))
aa1 = approx(fourVillage_16000_seasonality1yrperiod_mig01$time,fourVillage_16000_seasonality1yrperiod_mig01$E.D,xout=seq(0,10,0.01))
aa2 = approx(eightVillage_8000_seasonality1yrperiod_mig01$time,eightVillage_8000_seasonality1yrperiod_mig01$E.D,xout=seq(0,10,0.01))
aa3 = approx(sixteenVillage_4000_seasonality1yrperiod_mig01$time,sixteenVillage_4000_seasonality1yrperiod_mig01$E.D,xout=seq(0,10,0.01))
aa4 = approx(thirtytwoVillage_2000_seasonality1yrperiod_mig01$time,thirtytwoVillage_2000_seasonality1yrperiod_mig01$E.D,xout=seq(0,10,0.01))

a0$y[is.na(a0$y)]<-0
aa$y[is.na(aa$y)]<-0
aa1$y[is.na(aa1$y)]<-0
aa2$y[is.na(aa2$y)]<-0
aa3$y[is.na(aa3$y)]<-0
aa4$y[is.na(aa4$y)]<-0

marker=colorRampPalette(c('orange','red','purple','royalblue'))(7)


png(paste0(dir2,"psc_differential_metapop_paper_mig01_v1.png"), width=1400, height=800, res=150) 
par(fig=c(0,1,0,1),mar=c(3.3,3.6,1.1,1.1),oma=c(1,1,0,3),mgp=c(2.2,0.45,0),tcl=-0.4,cex=1.2)
plot(a0$x,(a0$y-aa$y),type='l',col=marker[[1]],lwd=2,ylim=c(-.15,.45),xlim=c(-.05,5),xlab="",ylab="Probability of silent circualation differential")
lines(a0$x,(a0$y-aa1$y),col=marker[[2]],lwd=2)
lines(a0$x,(a0$y-aa2$y),col=marker[[4]],lwd=2)
lines(a0$x,(a0$y-aa3$y),col=marker[[6]],lwd=2)
lines(a0$x,(a0$y-aa4$y),col=marker[[7]],lwd=2)
legend('topleft',legend=rev(c('2x32k','4x16k','8x8k','16x4k','32x2k')),col=rev(c(marker[1:2],marker[4],marker[6:7])),lwd=2,bty='n')
mtext(side=1,'Years since paralytic case',outer=TRUE,cex=1.2)
abline(h=0,lty=3)
mtext("1x64k\n has\n higher\n probability",side=4,at=.6,las=2,line=-1,outer=TRUE)
mtext("1x64k\n has\n lower\n probability",side=4,at=.25,las=2,line=-1,outer=TRUE)
par(new=TRUE)
par(mar=c(0,0,0,0))
par(fig=c(0.58,0.92,0.5,0.85))
plot(a0$x,(a0$y-aa$y),type='l',lwd=2,col=marker[[1]],ylim=c(-.035,.005),xlim=c(3,5),xlab="",ylab="")
lines(a0$x,(a0$y-aa1$y),col=marker[[2]],lwd=2)
lines(a0$x,(a0$y-aa2$y),col=marker[[4]],lwd=2)
lines(a0$x,(a0$y-aa3$y),col=marker[[6]],lwd=2)
lines(a0$x,(a0$y-aa4$y),col=marker[[7]],lwd=2)
abline(h=0,lty=3)
dev.off()




##with migration
fourVillage_16000_seasonality1yrperiod_mig01 = read.table(paste0(dir,'1pcase_4x16000_migRate_0.1_paper.csv'), col.names=c('time', 'E&D'))
twoVillage_32000_seasonality1yrperiod_mig01 = read.table(paste0(dir,'1pcase_2x32000_migRate_0.1_paper.csv'), col.names=c('time', 'E&D'))
eightVillage_8000_seasonality1yrperiod_mig01 = read.table(paste0(dir,'1pcase_8x8000_migRate_0.1_paper.csv'), col.names=c('time', 'E&D'))
sixteenVillage_4000_seasonality1yrperiod_mig01 = read.table(paste0(dir,'1pcase_16x4000_migRate_0.1_paper.csv'), col.names=c('time', 'E&D'))
mixed_1x32000_4x8000_seasonality1yrperiod_mig01 = read.table(paste0(dir,'1pcase_1x32000_4x8000_migRate_0.1_paper.csv'), col.names=c('time', 'E&D'))
mixed_1x32000_8x4000_seasonality1yrperiod_mig01 = read.table(paste0(dir,'1pcase_1x32000_8x4000_migRate_0.1_paper.csv'), col.names=c('time', 'E&D'))
thirtytwoVillage_2000_seasonality1yrperiod_mig01 = read.table(paste0(dir,'1pcase_32x2000_migRate_0.1_paper.csv'), col.names=c('time', 'E&D'))

a0 = approx(oneVillage_64000_seasonality1yrperiod$time,oneVillage_64000_seasonality1yrperiod$E.D,xout=seq(0,10,0.01))
a1 = approx(fourVillage_16000_seasonality1yrperiod_mig01$time,fourVillage_16000_seasonality1yrperiod_mig01$E.D,xout=seq(0,10,0.01))
a = approx(twoVillage_32000_seasonality1yrperiod_mig01$time,twoVillage_32000_seasonality1yrperiod_mig01$E.D,xout=seq(0,10,0.01))
a2 = approx(eightVillage_8000_seasonality1yrperiod_mig01$time,eightVillage_8000_seasonality1yrperiod_mig01$E.D,xout=seq(0,10,0.01))
a3 = approx(sixteenVillage_4000_seasonality1yrperiod_mig01$time,sixteenVillage_4000_seasonality1yrperiod_mig01$E.D,xout=seq(0,10,0.01))
a6 = approx(thirtytwoVillage_2000_seasonality1yrperiod_mig01$time,thirtytwoVillage_2000_seasonality1yrperiod_mig01$E.D,xout=seq(0,10,0.01))
a4 = approx(mixed_1x32000_4x8000_seasonality1yrperiod_mig01$time,mixed_1x32000_4x8000_seasonality1yrperiod_mig01$E.D,xout=seq(0,10,0.01))
a5 = approx(mixed_1x32000_8x4000_seasonality1yrperiod_mig01$time,mixed_1x32000_8x4000_seasonality1yrperiod_mig01$E.D,xout=seq(0,10,0.01))

a0$y[is.na(a0$y)] = 0
a1$y[is.na(a1$y)]= 0
a$y[is.na(a$y)]=0
a2$y[is.na(a2$y)]=0
a3$y[is.na(a3$y)]=0
a4$y[is.na(a4$y)]=0
a5$y[is.na(a5$y)]=0
a6$y[is.na(a6$y)]=0

marker=colorRampPalette(c('orange','red','purple','royalblue'))(7)


png(paste0(dir2,"SC_mig_comp_mixed_pop_v4.png"), width=1400, height=800, res=150) 
par(fig=c(0,1,0,1),mar=c(3.3,3.6,1.1,1.1),oma=c(1,1,0,3),mgp=c(2.2,0.45,0),tcl=-0.4)
plot(a0$x,(a0$y-a$y),type='l',col=marker[[1]],xlab="",ylab="Probability of silent circulation differential",xlim=c(-0.2,5),lwd=2,ylim=c(-.1,.5),xaxs="r",yaxs="r")
lines(a0$x,(a0$y-a1$y),col=marker[[2]],lwd=2)
lines(a0$x,(a0$y-a4$y),col=marker[[3]],lwd=2,lty=3)
lines(a0$x,(a0$y-a2$y) ,col=marker[[4]],lwd=2)
lines(a0$x,(a0$y-a5$y),col=marker[[5]],lwd=2,lty=6)
lines(a0$x,(a0$y-a3$y) ,col=marker[[6]],lwd=2)
lines(a0$x,(a0$y-a6$y),col=marker[[7]],lwd=2)
abline(h=0,col='black',lty=3)
legend('topleft',legend=c('32x2k','16x4k','1x32k,8x4k','8x8k','1x32k,4x8k','4x16k','2x32k'),
       col=rev(marker[1:7]),lty=c(1,1,3,1,6,1,1),lwd=2,bty='n')
#mtext("Probability of silent circulation differential",side=2,outer=TRUE)
mtext("Years since last detected paralytic case",side=1,outer=TRUE)
mtext("1x64k\n has\n higher\n probability",side=4,at=.6,las=2,line=-1,outer=TRUE)
mtext("1x64k\n has\n lower\n probability",side=4,at=.2,las=2,line=-1,outer=TRUE)
par(new=TRUE)
par(mar=c(0,0,0,0))
par(fig=c(0.58,0.92,0.5,0.85))
plot(a0$x,(a0$y-a$y),type='l',col=marker[[1]],xlab="",ylab="",xlim=c(3,5),lwd=2,ylim=c(-.04,.005))
lines(a0$x,(a0$y-a1$y),col=marker[[2]],lwd=2)
lines(a0$x,(a0$y-a4$y),col=marker[[3]],lwd=2,lty=3)
lines(a0$x,(a0$y-a2$y) ,col=marker[[4]],lwd=2)
lines(a0$x,(a0$y-a5$y),col=marker[[5]],lwd=2,lty=6)
lines(a0$x,(a0$y-a3$y) ,col=marker[[6]],lwd=2)
lines(a0$x,(a0$y-a6$y),col=marker[[7]],lwd=2)
abline(h=0,col='black',lty=3)
dev.off()

# comparing movement effects plots

sixteenVillage_4000_seasonality1yrperiod = read.table(paste0(dir,'1pcase_16x4000_migRate_0_paper.csv'), col.names=c('time', 'E&D'))
sixteenVillage_4000_seasonality1yrperiod_mig01 = read.table(paste0(dir,'1pcase_16x4000_migRate_0.1_paper.csv'), col.names=c('time', 'E&D'))
sixteenVillage_4000_seasonality1yrperiod_mig005 = read.table(paste0(dir,'1pcase_16x4000_migRate_0.05_paper.csv'), col.names=c('time', 'E&D'))
sixteenVillage_4000_seasonality1yrperiod_mig02 = read.table(paste0(dir,'1pcase_16x4000_migRate_0.2_paper.csv'), col.names=c('time', 'E&D'))
sixteenVillage_4000_seasonality1yrperiod_mig1 = read.table(paste0(dir,'1pcase_16x4000_migRate_1_paper.csv'), col.names=c('time', 'E&D'))


fourVillage_16000_seasonality1yrperiod = read.table(paste0(dir,'1pcase_4x16000_migRate_0_paper.csv'), col.names=c('time', 'E&D'))
fourVillage_16000_seasonality1yrperiod_mig01 = read.table(paste0(dir,'1pcase_4x16000_migRate_0.1_paper.csv'), col.names=c('time', 'E&D'))
fourVillage_16000_seasonality1yrperiod_mig005 = read.table(paste0(dir,'1pcase_4x16000_migRate_0.05_paper.csv'), col.names=c('time', 'E&D'))
fourVillage_16000_seasonality1yrperiod_mig02 = read.table(paste0(dir,'1pcase_4x16000_migRate_0.2_paper.csv'), col.names=c('time', 'E&D'))
fourVillage_16000_seasonality1yrperiod_mig1 = read.table(paste0(dir,'1pcase_4x16000_migRate_1_paper.csv'), col.names=c('time', 'E&D'))

a0 = approx(oneVillage_64000_seasonality1yrperiod$time,oneVillage_64000_seasonality1yrperiod$E.D,xout=seq(0,10,0.01))
a = approx(sixteenVillage_4000_seasonality1yrperiod$time,sixteenVillage_4000_seasonality1yrperiod$E.D,xout=seq(0,10,0.01))
a1 = approx(sixteenVillage_4000_seasonality1yrperiod_mig005$time,sixteenVillage_4000_seasonality1yrperiod_mig005$E.D,xout=seq(0,10,0.01))
a2 = approx(sixteenVillage_4000_seasonality1yrperiod_mig01$time,sixteenVillage_4000_seasonality1yrperiod_mig01$E.D,xout=seq(0,10,0.01))
a3 = approx(sixteenVillage_4000_seasonality1yrperiod_mig02$time,sixteenVillage_4000_seasonality1yrperiod_mig02$E.D,xout=seq(0,10,0.01))
a4 = approx(sixteenVillage_4000_seasonality1yrperiod_mig1$time,sixteenVillage_4000_seasonality1yrperiod_mig1$E.D,xout=seq(0,10,0.01))

a0$y[is.na(a0$y)] = 0
a1$y[is.na(a1$y)]= 0
a$y[is.na(a$y)]=0
a2$y[is.na(a2$y)]=0
a3$y[is.na(a3$y)]=0
a4$y[is.na(a4$y)]=0

marker1=brewer.pal(8,'Greens')



png(paste0(dir2,"SC_16x4000_mig_comp_v2.png"), width=1400, height=800, res=150) 
par(fig=c(0,1,0,1),mar=c(3.3,3.6,1.1,1.1),oma=c(1,0,0,3),mgp=c(2.2,0.45,0),tcl=-0.4,cex=1.2)
plot(sixteenVillage_4000_seasonality1yrperiod,type='l',col=marker1[4],xlab="",ylab="Probability of silent circulation",xlim=c(-.34,5),lwd=2)
lines(sixteenVillage_4000_seasonality1yrperiod_mig005,col=marker1[5],lwd=2)
lines(sixteenVillage_4000_seasonality1yrperiod_mig01,col=marker1[6],lwd=2)
lines(sixteenVillage_4000_seasonality1yrperiod_mig02_test,col=marker1[7],lwd=2)
lines(sixteenVillage_4000_seasonality1yrperiod_mig1,col=marker1[8],lwd=2)
lines(oneVillage_64000_seasonality1yrperiod,col='grey80',lwd=2)
legend('bottomleft',legend=c('Move rate = 0','Move rate = 0.05','Move rate = 0.1','Move rate = 0.2','Move rate = 1','1x64k'),col=c(marker1[4:8],'grey80'),lwd=2,bty='n')
mtext("Years since last detected paralytic case",side=1,outer=TRUE)
par(new=TRUE)
par(mar=c(0,0,0,0))
par(fig=c(0.6,0.95,0.5,0.85))
plot(sixteenVillage_4000_seasonality1yrperiod,type='l',col=marker1[4],xlab="",ylab="Probability of silent circulation",xlim=c(3,5),ylim=c(0,0.035),lwd=2)
lines(sixteenVillage_4000_seasonality1yrperiod_mig005,col=marker1[5],lwd=2)
lines(sixteenVillage_4000_seasonality1yrperiod_mig01,col=marker1[6],lwd=2)
lines(sixteenVillage_4000_seasonality1yrperiod_mig02_test,col=marker1[7],lwd=2)
lines(sixteenVillage_4000_seasonality1yrperiod_mig1,col=marker1[8],lwd=2)
lines(oneVillage_64000_seasonality1yrperiod,col='grey80',lwd=2)
dev.off()


a0 = approx(oneVillage_64000_seasonality1yrperiod$time,oneVillage_64000_seasonality1yrperiod$E.D,xout=seq(0,10,0.01))
a = approx(fourVillage_16000_seasonality1yrperiod$time,fourVillage_16000_seasonality1yrperiod$E.D,xout=seq(0,10,0.01))
a1 = approx(fourVillage_16000_seasonality1yrperiod_mig005$time,fourVillage_16000_seasonality1yrperiod_mig005$E.D,xout=seq(0,10,0.01))
a2 = approx(fourVillage_16000_seasonality1yrperiod_mig01$time,fourVillage_16000_seasonality1yrperiod_mig01$E.D,xout=seq(0,10,0.01))
a3 = approx(fourVillage_16000_seasonality1yrperiod_mig02$time,fourVillage_16000_seasonality1yrperiod_mig02$E.D,xout=seq(0,10,0.01))
a4 = approx(fourVillage_16000_seasonality1yrperiod_mig1$time,fourVillage_16000_seasonality1yrperiod_mig1$E.D,xout=seq(0,10,0.01))



png(paste0(dir2,"SC_4x16000_mig_comp.png"), width=1400, height=800, res=150) 
par(mfrow=c(1,2),oma=c(1,1,1,1),cex=1.2,mar=c(2.5,2.5,1,1))
plot(a0$x,(a0$y-a$y),type='l',col=marker1[4],xlab="",ylab="Probability of silent circulation differential",xlim=c(0,5),lwd=2)
lines(a0$x,(a0$y-a1$y),col=marker1[5],lwd=2)
lines(a0$x,(a0$y-a2$y) ,col=marker1[6],lwd=2)
lines(a0$x,(a0$y-a3$y) ,col=marker1[7],lwd=2)
lines(a0$x,(a0$y-a4$y),col=marker1[8],lwd=2)
abline(h=0,col='black',lty=3)
mtext("Probability of silent circulation differential",side=2,outer=TRUE)
plot(a0$x,(a0$y-a$y),type='l',col=marker1[4],xlab="",ylab="Probability of silent circulation differential",xlim=c(3,5),lwd=2)
lines(a0$x,(a0$y-a1$y),col=marker1[5],lwd=2)
lines(a0$x,(a0$y-a2$y) ,col=marker1[6],lwd=2)
lines(a0$x,(a0$y-a3$y) ,col=marker1[7],lwd=2)
lines(a0$x,(a0$y-a4$y),col=marker1[8],lwd=2)
legend('topright',legend=c('Movement rate = 0','Movement rate = 0.05','Movement rate = 0.1','Movement rate = 0.2','Movement rate = 1'),col=marker1[4:8],lwd=2,bty='n')
mtext("Years since last detected paralytic case",side=1,outer=TRUE)
dev.off()

png(paste0(dir2,"SC_4x16000_mig_comp_v4.png"), width=1400, height=800, res=150) 
par(fig=c(0,1,0,1),mar=c(3.3,3.6,1.1,1.1),oma=c(1,0,0,3),mgp=c(2.2,0.45,0),tcl=-0.4,cex=1.2)
plot(fourVillage_16000_seasonality1yrperiod,type='l',col=marker1[4],xlab="",ylab="Probability of silent circulation",xlim=c(-.34,5),lwd=2)
lines(fourVillage_16000_seasonality1yrperiod_mig005,col=marker1[5],lwd=2)
lines(fourVillage_16000_seasonality1yrperiod_mig01,col=marker1[6],lwd=2)
lines(fourVillage_16000_seasonality1yrperiod_mig02,col=marker1[7],lwd=2)
lines(fourVillage_16000_seasonality1yrperiod_mig1,col=marker1[8],lwd=2)
lines(oneVillage_64000_seasonality1yrperiod,col='grey80',lwd=2)
legend('bottomleft',legend=c('Move rate = 0','Move rate = 0.05','Move rate = 0.1','Move rate = 0.2','Move rate = 1','1x64k'),col=c(marker1[4:8],'grey80'),lwd=2,bty='n')
mtext("Years since last detected paralytic case",side=1,outer=TRUE)
par(new=TRUE)
par(mar=c(0,0,0,0))
par(fig=c(0.6,0.95,0.5,0.85))
plot(fourVillage_16000_seasonality1yrperiod,type='l',col=marker1[4],xlab="",xlim=c(3,5),lwd=2,ylim=c(0,.06))
lines(fourVillage_16000_seasonality1yrperiod_mig005,col=marker1[5],lwd=2)
lines(fourVillage_16000_seasonality1yrperiod_mig01,col=marker1[6],lwd=2)
lines(fourVillage_16000_seasonality1yrperiod_mig02,col=marker1[7],lwd=2)
lines(fourVillage_16000_seasonality1yrperiod_mig1,col=marker1[8],lwd=2)
lines(oneVillage_64000_seasonality1yrperiod,col='grey80',lwd=2)
dev.off()


