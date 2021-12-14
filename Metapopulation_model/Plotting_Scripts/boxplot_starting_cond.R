#boxplot for population distribution at beginning of sim after 50 year burn in

dir2 = 'Metapopulation_model/Figures/'
dir4 = '/Metapopulation_model/SC_statistic_output/'

#1 village of 64,000
oneVillage_64000_seasonality1yrperiod_S = read.csv(paste0(dir4,'S_at_burnin_64000reintRate_0.001000migRate_0.000000_paper.csv'),header=FALSE)
oneVillage_64000_seasonality1yrperiod_I1 = read.csv(paste0(dir4,'I1_at_burnin_64000reintRate_0.001000migRate_0.000000_paper.csv'),header=FALSE)
oneVillage_64000_seasonality1yrperiod_R = read.csv(paste0(dir4,'R_at_burnin_64000reintRate_0.001000migRate_0.000000_paper.csv'),header=FALSE)
oneVillage_64000_seasonality1yrperiod_P = read.csv(paste0(dir4,'P_at_burnin_64000reintRate_0.001000migRate_0.000000_paper.csv'),header=FALSE)
oneVillage_64000_seasonality1yrperiod_Ir = read.csv(paste0(dir4,'IR_at_burnin_64000reintRate_0.001000migRate_0.000000_paper.csv'),header=FALSE)

#filter out sims that didn't have pcase
extInt_64000 <-read.csv(paste0(dir4,"extinction_interval_64000reintRate_0.001000migRate_0.000000_paper.csv"),header=FALSE,fill=TRUE)
ni_64000<-which(is.na(extInt_64000[,1]))
#remove unused rows (i.e. sims that didn't have pcase)
oneVillage64000_S_new <- oneVillage_64000_seasonality1yrperiod_S[-ni_64000,]
oneVillage64000_I1_new <- oneVillage_64000_seasonality1yrperiod_I1[-ni_64000,]
oneVillage64000_IR_new <- oneVillage_64000_seasonality1yrperiod_Ir[-ni_64000,]
oneVillage64000_R_new <- oneVillage_64000_seasonality1yrperiod_R[-ni_64000,]
oneVillage64000_P_new <- oneVillage_64000_seasonality1yrperiod_P[-ni_64000,]


#1 village of 32,000
oneVillage_32000_seasonality1yrperiod_S = read.csv(paste0(dir4,'S_at_burnin_32000reintRate_0.001000migRate_0.000000_paper.csv'),header=FALSE)
oneVillage_32000_seasonality1yrperiod_I1 = read.csv(paste0(dir4,'I1_at_burnin_32000reintRate_0.001000migRate_0.000000_paper.csv'),header=FALSE)
oneVillage_32000_seasonality1yrperiod_R = read.csv(paste0(dir4,'R_at_burnin_32000reintRate_0.001000migRate_0.000000_paper.csv'),header=FALSE)
oneVillage_32000_seasonality1yrperiod_P = read.csv(paste0(dir4,'P_at_burnin_32000reintRate_0.001000migRate_0.000000_paper.csv'),header=FALSE)
oneVillage_32000_seasonality1yrperiod_Ir = read.csv(paste0(dir4,'IR_at_burnin_32000reintRate_0.001000migRate_0.000000_paper.csv'),header=FALSE)

extInt_32000 <-read.csv(paste0(dir4,"extinction_interval_32000reintRate_0.001000migRate_0.000000_paper.csv"),header=FALSE,fill=TRUE)
ni_32000<-which(is.na(extInt_32000[,1]))
#remove unused rows (i.e. sims that didn't have pcase)
oneVillage32000_S_new <- oneVillage_32000_seasonality1yrperiod_S[-ni_32000,]
oneVillage32000_I1_new <- oneVillage_32000_seasonality1yrperiod_I1[-ni_32000,]
oneVillage32000_IR_new <- oneVillage_32000_seasonality1yrperiod_Ir[-ni_32000,]
oneVillage32000_R_new <- oneVillage_32000_seasonality1yrperiod_R[-ni_32000,]
oneVillage32000_P_new <- oneVillage_32000_seasonality1yrperiod_P[-ni_32000,]


oneVillage_16000_seasonality1yrperiod_S = read.csv(paste0(dir4,'S_at_burnin_16000reintRate_0.001000migRate_0.000000_paper.csv'),header=FALSE)
oneVillage_16000_seasonality1yrperiod_I1 = read.csv(paste0(dir4,'I1_at_burnin_16000reintRate_0.001000migRate_0.000000_paper.csv'),header=FALSE)
oneVillage_16000_seasonality1yrperiod_R = read.csv(paste0(dir4,'R_at_burnin_16000reintRate_0.001000migRate_0.000000_paper.csv'),header=FALSE)
oneVillage_16000_seasonality1yrperiod_P = read.csv(paste0(dir4,'P_at_burnin_16000reintRate_0.001000migRate_0.000000_paper.csv'),header=FALSE)
oneVillage_16000_seasonality1yrperiod_Ir = read.csv(paste0(dir4,'IR_at_burnin_16000reintRate_0.001000migRate_0.000000_paper.csv'),header=FALSE)

extInt_16000 <-read.csv(paste0(dir4,"extinction_interval_16000reintRate_0.001000migRate_0.000000_paper.csv"),header=FALSE,fill=TRUE)
ni_16000<-which(is.na(extInt_16000[,1]))
#remove unused rows (i.e. sims that didn't have pcase)
oneVillage16000_S_new <- oneVillage_16000_seasonality1yrperiod_S[-ni_16000,]
oneVillage16000_I1_new <- oneVillage_16000_seasonality1yrperiod_I1[-ni_16000,]
oneVillage16000_IR_new <- oneVillage_16000_seasonality1yrperiod_Ir[-ni_16000,]
oneVillage16000_R_new <- oneVillage_16000_seasonality1yrperiod_R[-ni_16000,]
oneVillage16000_P_new <- oneVillage_16000_seasonality1yrperiod_P[-ni_16000,]


oneVillage_8000_seasonality1yrperiod_S = read.csv(paste0(dir4,'S_at_burnin_8000reintRate_0.001000migRate_0.000000_paper.csv'),header=FALSE)
oneVillage_8000_seasonality1yrperiod_I1 = read.csv(paste0(dir4,'I1_at_burnin_8000reintRate_0.001000migRate_0.000000_paper.csv'),header=FALSE)
oneVillage_8000_seasonality1yrperiod_R = read.csv(paste0(dir4,'R_at_burnin_8000reintRate_0.001000migRate_0.000000_paper.csv'),header=FALSE)
oneVillage_8000_seasonality1yrperiod_P = read.csv(paste0(dir4,'P_at_burnin_8000reintRate_0.001000migRate_0.000000_paper.csv'),header=FALSE)
oneVillage_8000_seasonality1yrperiod_Ir = read.csv(paste0(dir4,'IR_at_burnin_8000reintRate_0.001000migRate_0.000000_paper.csv'),header=FALSE)

extInt_8000 <-read.csv(paste0(dir4,"extinction_interval_8000reintRate_0.001000migRate_0.000000_paper.csv"),header=FALSE,fill=TRUE)
ni_8000<-which(is.na(extInt_8000[,1]))
#remove unused rows (i.e. sims that didn't have pcase)
oneVillage8000_S_new <- oneVillage_8000_seasonality1yrperiod_S[-ni_8000,]
oneVillage8000_I1_new <- oneVillage_8000_seasonality1yrperiod_I1[-ni_8000,]
oneVillage8000_IR_new <- oneVillage_8000_seasonality1yrperiod_Ir[-ni_8000,]
oneVillage8000_R_new <- oneVillage_8000_seasonality1yrperiod_R[-ni_8000,]
oneVillage8000_P_new <- oneVillage_8000_seasonality1yrperiod_P[-ni_8000,]


oneVillage_4000_seasonality1yrperiod_S = read.csv(paste0(dir4,'S_at_burnin_4000reintRate_0.001000migRate_0.000000_paper.csv'),header=FALSE)
oneVillage_4000_seasonality1yrperiod_I1 = read.csv(paste0(dir4,'I1_at_burnin_4000reintRate_0.001000migRate_0.000000_paper.csv'),header=FALSE)
oneVillage_4000_seasonality1yrperiod_R = read.csv(paste0(dir4,'R_at_burnin_4000reintRate_0.001000migRate_0.000000_paper.csv'),header=FALSE)
oneVillage_4000_seasonality1yrperiod_P = read.csv(paste0(dir4,'P_at_burnin_4000reintRate_0.001000migRate_0.000000_paper.csv'),header=FALSE)
oneVillage_4000_seasonality1yrperiod_Ir = read.csv(paste0(dir4,'IR_at_burnin_4000reintRate_0.001000migRate_0.000000_paper.csv'),header=FALSE)

extInt_4000 <-read.csv(paste0(dir4,"extinction_interval_4000reintRate_0.001000migRate_0.000000_paper.csv"),header=FALSE,fill=TRUE)
ni_4000<-which(is.na(extInt_4000[,1]))
#remove unused rows (i.e. sims that didn't have pcase)
oneVillage4000_S_new <- oneVillage_4000_seasonality1yrperiod_S[-ni_4000,]
oneVillage4000_I1_new <- oneVillage_4000_seasonality1yrperiod_I1[-ni_4000,]
oneVillage4000_IR_new <- oneVillage_4000_seasonality1yrperiod_Ir[-ni_4000,]
oneVillage4000_R_new <- oneVillage_4000_seasonality1yrperiod_R[-ni_4000,]
oneVillage4000_P_new <- oneVillage_4000_seasonality1yrperiod_P[-ni_4000,]


oneVillage_2000_seasonality1yrperiod_S = read.csv(paste0(dir4,'S_at_burnin_2000reintRate_0.001000migRate_0.000000_paper.csv'),header=FALSE)
oneVillage_2000_seasonality1yrperiod_I1 = read.csv(paste0(dir4,'I1_at_burnin_2000reintRate_0.001000migRate_0.000000_paper.csv'),header=FALSE)
oneVillage_2000_seasonality1yrperiod_R = read.csv(paste0(dir4,'R_at_burnin_2000reintRate_0.001000migRate_0.000000_paper.csv'),header=FALSE)
oneVillage_2000_seasonality1yrperiod_P = read.csv(paste0(dir4,'P_at_burnin_2000reintRate_0.001000migRate_0.000000_paper.csv'),header=FALSE)
oneVillage_2000_seasonality1yrperiod_Ir = read.csv(paste0(dir4,'IR_at_burnin_2000reintRate_0.001000migRate_0.000000_paper.csv'),header=FALSE)

extInt_2000 <-read.csv(paste0(dir4,"extinction_interval_2000reintRate_0.001000migRate_0.000000_paper.csv"),header=FALSE,fill=TRUE)
ni_2000<-which(is.na(extInt_2000[,1]))
#remove unused rows (i.e. sims that didn't have pcase)
oneVillage2000_S_new <- oneVillage_2000_seasonality1yrperiod_S[-ni_2000,]
oneVillage2000_I1_new <- oneVillage_2000_seasonality1yrperiod_I1[-ni_2000,]
oneVillage2000_IR_new <- oneVillage_2000_seasonality1yrperiod_Ir[-ni_2000,]
oneVillage2000_R_new <- oneVillage_2000_seasonality1yrperiod_R[-ni_2000,]
oneVillage2000_P_new <- oneVillage_2000_seasonality1yrperiod_P[-ni_2000,]


boxplot(oneVillage64000_S_new, oneVillage64000_I1_new,oneVillage64000_IR_new,oneVillage64000_R_new,oneVillage64000_P_new)
points(1,3707,col='red',pch=16)
points(2,93,col='red',pch=16)
points(3,147,col='red',pch=16)
points(4,26308,col='red',pch=16)
points(5,33745,col='red',pch=16)


boxplot(oneVillage32000_S_new, oneVillage32000_I1_new,oneVillage32000_IR_new,oneVillage32000_R_new,oneVillage32000_P_new)
points(1,1854,col='red',pch=16)
points(2,46,col='red',pch=16)
points(3,74,col='red',pch=16)
points(4,13154,col='red',pch=16)
points(5,16872,col='red',pch=16)

boxplot(oneVillage16000_S_new, oneVillage16000_I1_new,oneVillage16000_IR_new,oneVillage16000_R_new,oneVillage16000_P_new)
points(1,927,col='red',pch=16)
points(2,23,col='red',pch=16)
points(3,37,col='red',pch=16)
points(4,6577,col='red',pch=16)
points(5,8436,col='red',pch=16)

boxplot(oneVillage8000_S_new, oneVillage8000_I1_new,oneVillage8000_IR_new,oneVillage8000_R_new,oneVillage8000_P_new)
points(1,463,col='red',pch=16)
points(2,12,col='red',pch=16)
points(3,18,col='red',pch=16)
points(4,3289,col='red',pch=16)
points(5,4218,col='red',pch=16)

boxplot(oneVillage4000_S_new, oneVillage4000_I1_new,oneVillage4000_IR_new,oneVillage4000_R_new,oneVillage4000_P_new)
points(1,232,col='red',pch=16)
points(2,6,col='red',pch=16)
points(3,9,col='red',pch=16)
points(4,1644,col='red',pch=16)
points(5,2109,col='red',pch=16)

png(paste0(dir2,"boxplot_starting_S_v1.png"), width=1600, height=1600, res=180) 
par(oma=c(1,4,1,1),mgp=c(3.2,1,0),cex=1.6,mar=c(3,2.2,1,1))
boxplot(oneVillage64000_S_new,oneVillage32000_S_new,oneVillage16000_S_new,oneVillage8000_S_new,oneVillage4000_S_new,oneVillage2000_S_new,
        names=c("64000","32000","16000","8000","4000","2000"),ylab="Initial number in S compartment")
points(1,3707,pch=15,lwd=2,col='red')
points(2,1854,pch=15,lwd=2,col='red')
points(3,927,pch=15,lwd=2,col='red')
points(4,463,pch=15,lwd=2,col='red')
points(5,232,pch=15,lwd=2,col='red')
points(6,116,pch=15,lwd=2,col='red')
mtext(side=1,"Population size",line=3,cex=1.8)
mtext(side=2,"Initial number in S compartment",line=3,cex=1.8)
dev.off()

png(paste0(dir2,"boxplot_starting_I1_v1.png"), width=1600, height=1600, res=180) 
par(oma=c(1,4,1,1),mgp=c(3.2,1,0),cex=1.6,mar=c(3,2.2,1,1))
boxplot(oneVillage64000_I1_new,oneVillage32000_I1_new,oneVillage16000_I1_new,oneVillage8000_I1_new,oneVillage4000_I1_new,oneVillage2000_I1_new,
        names=c("64000","32000","16000","8000","4000","2000"),ylab="Initial number in I1 compartment")
points(1,93,pch=15,lwd=2,col='red')
points(2,46,pch=15,lwd=2,col='red')
points(3,23,pch=15,lwd=2,col='red')
points(4,6,pch=15,lwd=2,col='red')
points(5,6,pch=15,lwd=2,col='red')
points(6,3,pch=15,lwd=2,col='red')
mtext(side=1,"Population size",line=3,cex=1.8)
mtext(side=2,"Initial number in I1 compartment",line=3,cex=1.8)
dev.off()

png(paste0(dir2,"boxplot_starting_R_v1.png"), width=1600, height=1600, res=180) 
par(oma=c(1,4,1,1),mgp=c(3.2,1,0),cex=1.6,mar=c(3,2.2,1,1))
boxplot(oneVillage64000_R_new,oneVillage32000_R_new,oneVillage16000_R_new,oneVillage8000_R_new,oneVillage4000_R_new,oneVillage2000_R_new,
        names=c("64000","32000","16000","8000","4000","2000"),ylab="Initial number in R compartment")
points(1,26308,pch=15,lwd=2,col='red')
points(2,13154,pch=15,lwd=2,col='red')
points(3,6577,pch=15,lwd=2,col='red')
points(4,3289,pch=15,lwd=2,col='red')
points(5,1644,pch=15,lwd=2,col='red')
points(6,822,pch=15,lwd=2,col='red')
mtext(side=1,"Population size",line=3,cex=1.8)
mtext(side=2,"Initial number in R compartment",line=3,cex=1.8)
dev.off()

png(paste0(dir2,"boxplot_starting_P_v1.png"), width=1600, height=1600, res=180) 
par(oma=c(1,4,1,1),mgp=c(3.2,1,0),cex=1.6,mar=c(3,2.2,1,1))
boxplot(oneVillage64000_P_new,oneVillage32000_P_new,oneVillage16000_P_new,oneVillage8000_P_new,oneVillage4000_P_new,oneVillage2000_P_new,
        names=c("64000","32000","16000","8000","4000","2000"),ylab="Initial number in P compartment")
points(1,33745,pch=15,lwd=2,col='red')
points(2,16872,pch=15,lwd=2,col='red')
points(3,8436,pch=15,lwd=2,col='red')
points(4,4218,pch=15,lwd=2,col='red')
points(5,2109,pch=15,lwd=2,col='red')
points(6,1054,pch=15,lwd=2,col='red')
mtext(side=1,"Population size",line=3,cex=1.8)
mtext(side=2,"Initial number in P compartment",line=3,cex=1.8)
dev.off()

png(paste0(dir2,"boxplot_starting_IR_v1.png"), width=1600, height=1600, res=180) 
par(oma=c(1,4,1,1),mgp=c(3.2,1,0),cex=1.6,mar=c(3,2.2,1,1))
boxplot(oneVillage64000_IR_new,oneVillage32000_IR_new,oneVillage16000_IR_new,oneVillage8000_IR_new,oneVillage4000_IR_new,oneVillage2000_IR_new,
        names=c("64000","32000","16000","8000","4000","2000"),ylab="Initial number in Ir compartment")
points(1,147,pch=15,lwd=2,col='red')
points(2,74,pch=15,lwd=2,col='red')
points(3,37,pch=15,lwd=2,col='red')
points(4,18,pch=15,lwd=2,col='red')
points(5,9,pch=15,lwd=2,col='red')
points(6,5,pch=15,lwd=2,col='red')
mtext(side=1,"Population size",line=3,cex=1.8)
mtext(side=2,"Initial number in Ir compartment",line=3,cex=1.8)
dev.off()


