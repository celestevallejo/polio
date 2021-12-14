#compare number of infections between case to geometric dist
#compare number of infections between last case and extinction to geometric dist

dir5 = "/Users/Celeste/Desktop/multipatch_model/sim_results/"
dir3= "/Users/Celeste/Desktop/multipatch_model/"
dir2 = "/Users/Celeste/Desktop/multipatch_model/Figures/"
dir4 = '/Users/Celeste/Desktop/multipatch_model/metapopulation_polio/'

oneVillage64000_pcase_I1 = read.csv(paste0(dir4,"num_I1_inf_bet_case_64000reintRate_0.001000migRate_0.000000_paper.csv"),sep=",",header = FALSE,col.names = paste0("V",seq_len(9000)), fill=TRUE)
oneVillage64000_pcase_IR = read.csv(paste0(dir4,"num_IR_inf_bet_case_64000reintRate_0.001000migRate_0.000000_paper.csv"),sep=",",header = FALSE,col.names = paste0("V",seq_len(9000)), fill=TRUE)
oneVillage64000_pcase_ext_I1 = read.csv(paste0(dir4,"num_I1_inf_bet_case_ext_64000reintRate_0.001000migRate_0.000000_paper.csv"),sep=",",header = FALSE,col.names = paste0("V",seq_len(2)), fill=TRUE)
oneVillage64000_pcase_ext_IR = read.csv(paste0(dir4,"num_IR_inf_bet_case_ext_64000reintRate_0.001000migRate_0.000000_paper.csv"),sep=",",header = FALSE,col.names = paste0("V",seq_len(2)), fill=TRUE)
extInt_64000 <-read.csv(paste0(dir4,"extinction_interval_64000reintRate_0.001000migRate_0.000000_paper.csv"),header=FALSE,fill=TRUE)
ni_64000<-which(is.na(extInt_64000[,1]))

oneVillage64000_pcase_I1_new = oneVillage64000_pcase_I1[-ni_64000,]
oneVillage64000_pcase_IR_new = oneVillage64000_pcase_IR[-ni_64000,]
oneVillage64000_pcase_ext_I1_new = oneVillage64000_pcase_ext_I1[-ni_64000,]
oneVillage64000_pcase_ext_IR_new = oneVillage64000_pcase_ext_IR[-ni_64000,]

oneVillage32000_pcase_I1 = read.csv(paste0(dir4,"num_I1_inf_bet_case_32000reintRate_0.001000migRate_0.000000_paper.csv"),sep=",",header = FALSE,col.names = paste0("V",seq_len(1000)), fill=TRUE)
oneVillage32000_pcase_IR = read.csv(paste0(dir4,"num_IR_inf_bet_case_32000reintRate_0.001000migRate_0.000000_paper.csv"),sep=",",header = FALSE,col.names = paste0("V",seq_len(1000)), fill=TRUE)
oneVillage32000_pcase_ext_I1 = read.csv(paste0(dir4,"num_I1_inf_bet_case_ext_32000reintRate_0.001000migRate_0.000000_paper.csv"),sep=",",header = FALSE,col.names = paste0("V",seq_len(2)), fill=TRUE)
oneVillage32000_pcase_ext_IR = read.csv(paste0(dir4,"num_IR_inf_bet_case_ext_32000reintRate_0.001000migRate_0.000000_paper.csv"),sep=",",header = FALSE,col.names = paste0("V",seq_len(2)), fill=TRUE)
extInt_32000 <-read.csv(paste0(dir4,"extinction_interval_32000reintRate_0.001000migRate_0.000000_paper.csv"),header=FALSE,fill=TRUE)
ni_32000<-which(is.na(extInt_32000[,1]))

oneVillage32000_pcase_I1_new = oneVillage32000_pcase_I1[-ni_32000,]
oneVillage32000_pcase_IR_new = oneVillage32000_pcase_IR[-ni_32000,]
oneVillage32000_pcase_ext_I1_new = oneVillage32000_pcase_ext_I1[-ni_32000,]
oneVillage32000_pcase_ext_IR_new = oneVillage32000_pcase_ext_IR[-ni_32000,]

oneVillage16000_pcase_I1 = read.csv(paste0(dir4,"num_I1_inf_bet_case_16000reintRate_0.001000migRate_0.000000_paper.csv"),sep=",",header = FALSE,col.names = paste0("V",seq_len(9000)), fill=TRUE)
oneVillage16000_pcase_IR = read.csv(paste0(dir4,"num_IR_inf_bet_case_16000reintRate_0.001000migRate_0.000000_paper.csv"),sep=",",header = FALSE,col.names = paste0("V",seq_len(9000)), fill=TRUE)
oneVillage16000_pcase_ext_I1 = read.csv(paste0(dir4,"num_I1_inf_bet_case_ext_16000reintRate_0.001000migRate_0.000000_paper.csv"),sep=",",header = FALSE,col.names = paste0("V",seq_len(2)), fill=TRUE)
oneVillage16000_pcase_ext_IR = read.csv(paste0(dir4,"num_IR_inf_bet_case_ext_16000reintRate_0.001000migRate_0.000000_paper.csv"),sep=",",header = FALSE,col.names = paste0("V",seq_len(2)), fill=TRUE)
extInt_16000 <-read.csv(paste0(dir4,"extinction_interval_16000reintRate_0.001000migRate_0.000000_paper.csv"),header=FALSE,fill=TRUE)
ni_16000<-which(is.na(extInt_16000[,1]))

oneVillage16000_pcase_I1_mig01 = read.csv(paste0(dir4,"num_I1_inf_bet_case_16000160001600016000reintRate_0.001000migRate_0.100000_paper.csv"),sep=",",header = FALSE,col.names = paste0("V",seq_len(9000)), fill=TRUE)
oneVillage16000_pcase_IR_mig01 = read.csv(paste0(dir4,"num_IR_inf_bet_case_16000160001600016000reintRate_0.001000migRate_0.100000_paper.csv"),sep=",",header = FALSE,col.names = paste0("V",seq_len(9000)), fill=TRUE)
oneVillage16000_pcase_ext_I1_mig01 = read.csv(paste0(dir4,"num_I1_inf_bet_case_ext_16000160001600016000reintRate_0.001000migRate_0.100000_paper.csv"),sep=",",header = FALSE,col.names = paste0("V",seq_len(2)), fill=TRUE)
oneVillage16000_pcase_ext_IR_mig01 = read.csv(paste0(dir4,"num_IR_inf_bet_case_ext_16000160001600016000reintRate_0.001000migRate_0.100000_paper.csv"),sep=",",header = FALSE,col.names = paste0("V",seq_len(2)), fill=TRUE)
extInt_16000_mig01 <-read.csv(paste0(dir4,"extinction_interval_16000160001600016000reintRate_0.001000migRate_0.100000_paper.csv"),header=FALSE,fill=TRUE)
ni_16000_mig01<-which(is.na(extInt_16000_mig01[,1]))

oneVillage16000_pcase_I1_mig0 = read.csv(paste0(dir4,"num_I1_inf_bet_case_16000160001600016000reintRate_0.001000migRate_0.000000_paper.csv"),sep=",",header = FALSE,col.names = paste0("V",seq_len(9000)), fill=TRUE)
oneVillage16000_pcase_IR_mig0 = read.csv(paste0(dir4,"num_IR_inf_bet_case_16000160001600016000reintRate_0.001000migRate_0.000000_paper.csv"),sep=",",header = FALSE,col.names = paste0("V",seq_len(9000)), fill=TRUE)
oneVillage16000_pcase_ext_I1_mig0 = read.csv(paste0(dir4,"num_I1_inf_bet_case_ext_16000160001600016000reintRate_0.001000migRate_0.000000_paper.csv"),sep=",",header = FALSE,col.names = paste0("V",seq_len(2)), fill=TRUE)
oneVillage16000_pcase_ext_IR_mig0 = read.csv(paste0(dir4,"num_IR_inf_bet_case_ext_16000160001600016000reintRate_0.001000migRate_0.000000_paper.csv"),sep=",",header = FALSE,col.names = paste0("V",seq_len(2)), fill=TRUE)
extInt_16000_mig0 <-read.csv(paste0(dir4,"extinction_interval_16000160001600016000reintRate_0.001000migRate_0.000000_paper.csv"),header=FALSE,fill=TRUE)
ni_16000_mig0<-which(is.na(extInt_16000_mig0[,1]))


oneVillage16000_pcase_I1_new = oneVillage16000_pcase_I1[-ni_16000,]
oneVillage16000_pcase_IR_new = oneVillage16000_pcase_IR[-ni_16000,]
oneVillage16000_pcase_ext_I1_new = oneVillage16000_pcase_ext_I1[-ni_16000,]
oneVillage16000_pcase_ext_IR_new = oneVillage16000_pcase_ext_IR[-ni_16000,]

oneVillage16000_pcase_I1_new_mig01 = oneVillage16000_pcase_I1_mig01[-ni_16000_mig01,]
oneVillage16000_pcase_IR_new_mig01 = oneVillage16000_pcase_IR_mig01[-ni_16000_mig01,]
oneVillage16000_pcase_ext_I1_new_mig01 = oneVillage16000_pcase_ext_I1_mig01[-ni_16000_mig01,]
oneVillage16000_pcase_ext_IR_new_mig01 = oneVillage16000_pcase_ext_IR_mig01[-ni_16000_mig01,]

oneVillage16000_pcase_I1_new_mig0 = oneVillage16000_pcase_I1_mig0[-ni_16000_mig0,]
oneVillage16000_pcase_IR_new_mig0 = oneVillage16000_pcase_IR_mig0[-ni_16000_mig0,]
oneVillage16000_pcase_ext_I1_new_mig0 = oneVillage16000_pcase_ext_I1_mig0[-ni_16000_mig0,]
oneVillage16000_pcase_ext_IR_new_mig0 = oneVillage16000_pcase_ext_IR_mig0[-ni_16000_mig0,]


oneVillage8000_pcase_I1 = read.csv(paste0(dir4,"num_I1_inf_bet_case_8000reintRate_0.001000migRate_0.000000_paper.csv"),sep=",",header = FALSE,col.names = paste0("V",seq_len(9000)), fill=TRUE)
oneVillage8000_pcase_IR = read.csv(paste0(dir4,"num_IR_inf_bet_case_8000reintRate_0.001000migRate_0.000000_paper.csv"),sep=",",header = FALSE,col.names = paste0("V",seq_len(9000)), fill=TRUE)
oneVillage8000_pcase_ext_I1 = read.csv(paste0(dir4,"num_I1_inf_bet_case_ext_8000reintRate_0.001000migRate_0.000000_paper.csv"),sep=",",header = FALSE,col.names = paste0("V",seq_len(2)), fill=TRUE)
oneVillage8000_pcase_ext_IR = read.csv(paste0(dir4,"num_IR_inf_bet_case_ext_8000reintRate_0.001000migRate_0.000000_paper.csv"),sep=",",header = FALSE,col.names = paste0("V",seq_len(2)), fill=TRUE)
extInt_8000 <-read.csv(paste0(dir4,"extinction_interval_8000reintRate_0.001000migRate_0.000000_paper.csv"),header=FALSE,fill=TRUE)
ni_8000<-which(is.na(extInt_8000[,1]))

oneVillage8000_pcase_I1_new = oneVillage8000_pcase_I1[-ni_8000,]
oneVillage8000_pcase_IR_new = oneVillage8000_pcase_IR[-ni_8000,]
oneVillage8000_pcase_ext_I1_new = oneVillage8000_pcase_ext_I1[-ni_8000,]
oneVillage8000_pcase_ext_IR_new = oneVillage8000_pcase_ext_IR[-ni_8000,]

oneVillage4000_pcase_I1 = read.csv(paste0(dir4,"num_I1_inf_bet_case_4000reintRate_0.001000migRate_0.000000_paper.csv"),sep=",",header = FALSE,col.names = paste0("V",seq_len(1000)), fill=TRUE)
oneVillage4000_pcase_IR = read.csv(paste0(dir4,"num_IR_inf_bet_case_4000reintRate_0.001000migRate_0.000000_paper.csv"),sep=",",header = FALSE,col.names = paste0("V",seq_len(1000)), fill=TRUE)
oneVillage4000_pcase_ext_I1 = read.csv(paste0(dir4,"num_I1_inf_bet_case_ext_4000reintRate_0.001000migRate_0.000000_paper.csv"),sep=",",header = FALSE,col.names = paste0("V",seq_len(2)), fill=TRUE)
oneVillage4000_pcase_ext_IR = read.csv(paste0(dir4,"num_IR_inf_bet_case_ext_4000reintRate_0.001000migRate_0.000000_paper.csv"),sep=",",header = FALSE,col.names = paste0("V",seq_len(2)), fill=TRUE)
extInt_4000 <-read.csv(paste0(dir4,"extinction_interval_4000reintRate_0.001000migRate_0.000000_paper.csv"),header=FALSE,fill=TRUE)
ni_4000<-which(is.na(extInt_4000[,1]))

oneVillage4000_pcase_I1_new = oneVillage4000_pcase_I1[-ni_4000,]
oneVillage4000_pcase_IR_new = oneVillage4000_pcase_IR[-ni_4000,]
oneVillage4000_pcase_ext_I1_new = oneVillage4000_pcase_ext_I1[-ni_4000,]
oneVillage4000_pcase_ext_IR_new = oneVillage4000_pcase_ext_IR[-ni_4000,]

oneVillage4000_pcase_I1_mig0 = read.csv(paste0(dir4,"num_I1_inf_bet_case_ext_4000400040004000400040004000400040004000400040004000400040004000reintRate_0.001000migRate_0.000000_paper.csv"),sep=",",header = FALSE,col.names = paste0("V",seq_len(4000)), fill=TRUE)
oneVillage4000_pcase_I1_mig01 = read.csv(paste0(dir4,"num_I1_inf_bet_case_ext_4000400040004000400040004000400040004000400040004000400040004000reintRate_0.001000migRate_0.100000_paper.csv"),sep=",",header = FALSE,col.names = paste0("V",seq_len(4000)), fill=TRUE)
extInt_4000_mig0 <-read.csv(paste0(dir4,"extinction_interval_4000400040004000400040004000400040004000400040004000400040004000reintRate_0.001000migRate_0.000000_paper.csv"),header=FALSE,fill=TRUE)
ni_4000_mig0<-which(is.na(extInt_4000_mig0[,1]))

extInt_4000_mig01 <-read.csv(paste0(dir4,"extinction_interval_4000400040004000400040004000400040004000400040004000400040004000reintRate_0.001000migRate_0.100000_paper.csv"),header=FALSE,fill=TRUE)
ni_4000_mig01<-which(is.na(extInt_4000_mig01[,1]))

oneVillage4000_pcase_I1_new_mig0 = oneVillage4000_pcase_I1_mig0[-ni_4000_mig0,]
oneVillage4000_pcase_I1_new_mig01 = oneVillage4000_pcase_I1_mig01[-ni_4000_mig01,]

twoVillage32000_pcase_I1 = read.csv(paste0(dir4,"num_I1_inf_bet_case_3200032000reintRate_0.001000migRate_0.000000_paper.csv"),sep=",",header = FALSE,col.names = paste0("V",seq_len(5000)), fill=TRUE)
twoVillage32000_pcase_IR = read.csv(paste0(dir4,"num_IR_inf_bet_case_3200032000reintRate_0.001000migRate_0.000000_paper.csv"),sep=",",header = FALSE,col.names = paste0("V",seq_len(5000)), fill=TRUE)
twoVillage32000_pcase_ext_I1 = read.csv(paste0(dir4,"num_I1_inf_bet_case_ext_3200032000reintRate_0.001000migRate_0.000000_paper.csv"),sep=",",header = FALSE,col.names = paste0("V",seq_len(5000)), fill=TRUE)
twoVillage32000_pcase_ext_IR = read.csv(paste0(dir4,"num_IR_inf_bet_case_ext_3200032000reintRate_0.001000migRate_0.000000_paper.csv"),sep=",",header = FALSE,col.names = paste0("V",seq_len(5000)), fill=TRUE)

extInt_32000 <-read.csv(paste0(dir4,"extinction_interval_3200032000reintRate_0.001000migRate_0.000000_paper.csv"),header=FALSE,fill=TRUE)
ni_32000<-which(is.na(extInt_32000[,1]))

twoVillage32000_pcase_ext_I1_new = twoVillage32000_pcase_ext_I1[-ni_32000,]
twoVillage32000_pcase_ext_IR_new = twoVillage32000_pcase_ext_IR[-ni_32000,]
twoVillage32000_pcase_I1_new = twoVillage32000_pcase_I1[-ni_32000,]
twoVillage32000_pcase_IR_new = twoVillage32000_pcase_IR[-ni_32000,]

fourVillage16000_pcase_I1 = read.csv(paste0(dir4,"num_I1_inf_bet_case_16000160001600016000reintRate_0.001000migRate_0.000000_paper.csv"),sep=",",header = FALSE,col.names = paste0("V",seq_len(5000)), fill=TRUE)
fourVillage16000_pcase_IR = read.csv(paste0(dir4,"num_IR_inf_bet_case_16000160001600016000reintRate_0.001000migRate_0.000000_paper.csv"),sep=",",header = FALSE,col.names = paste0("V",seq_len(5000)), fill=TRUE)
fourVillage16000_pcase_ext_I1 = read.csv(paste0(dir4,"num_I1_inf_bet_case_ext_16000160001600016000reintRate_0.001000migRate_0.000000_paper.csv"),sep=",",header = FALSE,col.names = paste0("V",seq_len(5000)), fill=TRUE)
fourVillage16000_pcase_ext_IR = read.csv(paste0(dir4,"num_IR_inf_bet_case_ext_16000160001600016000reintRate_0.001000migRate_0.000000_paper.csv"),sep=",",header = FALSE,col.names = paste0("V",seq_len(5000)), fill=TRUE)
extInt_16000 <-read.csv(paste0(dir4,"extinction_interval_16000160001600016000reintRate_0.001000migRate_0.000000_paper.csv"),header=FALSE,fill=TRUE)
ni<-which(is.na(extInt_16000[,1]))

fourVillage16000_pcase_I1_new <-fourVillage16000_pcase_I1[-ni,]
fourVillage16000_pcase_IR_new <-fourVillage16000_pcase_IR[-ni,]
fourVillage16000_pcase_ext_I1_new <-fourVillage16000_pcase_ext_I1[-ni,]
fourVillage16000_pcase_ext_IR_new <-fourVillage16000_pcase_ext_IR[-ni,]

eightVillage8000_pcase_I1 = read.csv(paste0(dir4,"num_I1_inf_bet_case_80008000800080008000800080008000reintRate_0.001000migRate_0.000000_paper.csv"),sep=",",header = FALSE,col.names = paste0("V",seq_len(5000)), fill=TRUE)
eightVillage8000_pcase_IR = read.csv(paste0(dir4,"num_IR_inf_bet_case_80008000800080008000800080008000reintRate_0.001000migRate_0.000000_paper.csv"),sep=",",header = FALSE,col.names = paste0("V",seq_len(5000)), fill=TRUE)
eightVillage8000_pcase_ext_I1 = read.csv(paste0(dir4,"num_I1_inf_bet_case_ext_80008000800080008000800080008000reintRate_0.001000migRate_0.000000_paper.csv"),sep=",",header = FALSE,col.names = paste0("V",seq_len(5000)), fill=TRUE)
eightVillage8000_pcase_ext_IR = read.csv(paste0(dir4,"num_IR_inf_bet_case_ext_80008000800080008000800080008000reintRate_0.001000migRate_0.000000_paper.csv"),sep=",",header = FALSE,col.names = paste0("V",seq_len(5000)), fill=TRUE)
extInt_8000 <-read.csv(paste0(dir4,"extinction_interval_80008000800080008000800080008000reintRate_0.001000migRate_0.000000_paper.csv"),header=FALSE,fill=TRUE)
ni<-which(is.na(extInt_8000[,1]))

eightVillage8000_pcase_I1_new <-eightVillage8000_pcase_I1[-ni,]
eightVillage8000_pcase_IR_new <-eightVillage8000_pcase_IR[-ni,]
eightVillage8000_pcase_ext_I1_new <-eightVillage8000_pcase_ext_I1[-ni,]
eightVillage8000_pcase_ext_IR_new <-eightVillage8000_pcase_ext_IR[-ni,]

sixteenVillage4000_pcase_I1 = read.csv(paste0(dir4,"num_I1_inf_bet_case_4000400040004000400040004000400040004000400040004000400040004000reintRate_0.001000migRate_0.000000_paper.csv"),sep=",",header = FALSE,col.names = paste0("V",seq_len(5000)), fill=TRUE)
sixteenVillage4000_pcase_IR = read.csv(paste0(dir4,"num_IR_inf_bet_case_4000400040004000400040004000400040004000400040004000400040004000reintRate_0.001000migRate_0.000000_paper.csv"),sep=",",header = FALSE,col.names = paste0("V",seq_len(5000)), fill=TRUE)
sixteenVillage4000_pcase_ext_I1 = read.csv(paste0(dir4,"num_I1_inf_bet_case_ext_4000400040004000400040004000400040004000400040004000400040004000reintRate_0.001000migRate_0.000000_paper.csv"),sep=",",header = FALSE,col.names = paste0("V",seq_len(5000)), fill=TRUE)
sixteenVillage4000_pcase_ext_IR = read.csv(paste0(dir4,"num_IR_inf_bet_case_ext_4000400040004000400040004000400040004000400040004000400040004000reintRate_0.001000migRate_0.000000_paper.csv"),sep=",",header = FALSE,col.names = paste0("V",seq_len(5000)), fill=TRUE)
extInt_4000 <-read.csv(paste0(dir4,"extinction_interval_4000400040004000400040004000400040004000400040004000400040004000reintRate_0.001000migRate_0.000000_paper.csv"),header=FALSE,fill=TRUE)
ni<-which(is.na(extInt_4000[,1]))

sixteenVillage4000_pcase_I1_new <-sixteenVillage4000_pcase_I1[-ni,]
sixteenVillage4000_pcase_IR_new <-sixteenVillage4000_pcase_IR[-ni,]
sixteenVillage4000_pcase_ext_I1_new <-sixteenVillage4000_pcase_ext_I1[-ni,]
sixteenVillage4000_pcase_ext_IR_new <-sixteenVillage4000_pcase_ext_IR[-ni,]

thirtytwoVillage2000_pcase_I1 = read.csv(paste0(dir4,"num_I1_inf_bet_case_20002000200020002000200020002000200020002000200020002000200020002000200020002000200020002000200020002000200020002000200020002000reintRate_0.001000migRate_0.000000_paper.csv"),sep=",",header = FALSE,col.names = paste0("V",seq_len(5000)), fill=TRUE)
thirtytwoVillage2000_pcase_IR = read.csv(paste0(dir4,"num_IR_inf_bet_case_20002000200020002000200020002000200020002000200020002000200020002000200020002000200020002000200020002000200020002000200020002000reintRate_0.001000migRate_0.000000_paper.csv"),sep=",",header = FALSE,col.names = paste0("V",seq_len(5000)), fill=TRUE)
thirtytwoVillage2000_pcase_ext_I1 = read.csv(paste0(dir4,"num_I1_inf_bet_case_ext_20002000200020002000200020002000200020002000200020002000200020002000200020002000200020002000200020002000200020002000200020002000reintRate_0.001000migRate_0.000000_paper.csv"),sep=",",header = FALSE,col.names = paste0("V",seq_len(5000)), fill=TRUE)
thirtytwoVillage2000_pcase_ext_IR = read.csv(paste0(dir4,"num_IR_inf_bet_case_ext_20002000200020002000200020002000200020002000200020002000200020002000200020002000200020002000200020002000200020002000200020002000reintRate_0.001000migRate_0.000000_paper.csv"),sep=",",header = FALSE,col.names = paste0("V",seq_len(5000)), fill=TRUE)
extInt_2000 <-read.csv(paste0(dir4,"extinction_interval_32x2000reintRate_0.001000migRate_0.000000_paper.csv"),header=FALSE,fill=TRUE)
ni<-which(is.na(extInt_2000[,1]))

thirtytwoVillage2000_pcase_I1_new <-thirtytwoVillage2000_pcase_I1[-ni,]
thirtytwoVillage2000_pcase_IR_new <-thirtytwoVillage2000_pcase_IR[-ni,]
thirtytwoVillage2000_pcase_ext_I1_new <-thirtytwoVillage2000_pcase_ext_I1[-ni,]
thirtytwoVillage2000_pcase_ext_IR_new <-thirtytwoVillage2000_pcase_ext_IR[-ni,]



fourVillage16000_pcase_I1 = read.csv(paste0(dir5,"num_I1_inf_bet_case_16000160001600016000reintRate_0.001000migRate_0.000000_seasonality_period_1yr_1pcase_inf_bet_case.csv"),sep=",",header = FALSE,col.names = paste0("V",seq_len(2000)), fill=TRUE)
fourVillage16000_pcase_IR = read.csv(paste0(dir5,"num_IR_inf_bet_case_16000160001600016000reintRate_0.001000migRate_0.000000_seasonality_period_1yr_1pcase_inf_bet_case.csv"),sep=",",header = FALSE,col.names = paste0("V",seq_len(2000)), fill=TRUE)
fourVillage16000_pcase_ext_I1 = read.csv(paste0(dir5,"num_I1_inf_bet_case_ext_16000160001600016000reintRate_0.001000migRate_0.000000_seasonality_period_1yr_1pcase_inf_bet_case.csv"),sep=",",header = FALSE,col.names = paste0("V",seq_len(2000)), fill=TRUE)
fourVillage16000_pcase_ext_IR = read.csv(paste0(dir5,"num_IR_inf_bet_case_ext_16000160001600016000reintRate_0.001000migRate_0.000000_seasonality_period_1yr_1pcase_inf_bet_case.csv"),sep=",",header = FALSE,col.names = paste0("V",seq_len(2000)), fill=TRUE)

fourVillage16000_pcase_I1_mig_0.1 = read.csv(paste0(dir5,"num_I1_inf_bet_case_16000160001600016000reintRate_0.001000migRate_0.100000_seasonality_period_1yr_1pcase_inf_bet_case.csv"),sep=",",header = FALSE,col.names = paste0("V",seq_len(2000)), fill=TRUE)
fourVillage16000_pcase_IR_mig_0.1 = read.csv(paste0(dir5,"num_IR_inf_bet_case_16000160001600016000reintRate_0.001000migRate_0.100000_seasonality_period_1yr_1pcase_inf_bet_case.csv"),sep=",",header = FALSE,col.names = paste0("V",seq_len(2000)), fill=TRUE)
fourVillage16000_pcase_ext_I1_mig_0.1 = read.csv(paste0(dir5,"num_I1_inf_bet_case_ext_16000160001600016000reintRate_0.001000migRate_0.100000_seasonality_period_1yr_1pcase_inf_bet_case.csv"),sep=",",header = FALSE,col.names = paste0("V",seq_len(2000)), fill=TRUE)
fourVillage16000_pcase_ext_IR_mig_0.1 = read.csv(paste0(dir5,"num_IR_inf_bet_case_ext_16000160001600016000reintRate_0.001000migRate_0.100000_seasonality_period_1yr_1pcase_inf_bet_case.csv"),sep=",",header = FALSE,col.names = paste0("V",seq_len(2000)), fill=TRUE)


extInt_16000 <-read.csv(paste0(dir5,"extinction_interval_16000reintRate_0.001000migRate_0.000000_seasonality_period_1yr_1pcase_inf_bet_case.csv"),header=FALSE,fill=TRUE)
ni<-which(is.na(extInt_16000[,1]))
oneVillage16000_pcase_I1_new <-oneVillage16000_pcase_I1[-ni,]
oneVillage16000_pcase_IR_new <-oneVillage16000_pcase_IR[-ni,]
oneVillage16000_pcase_ext_I1_new <-oneVillage16000_pcase_ext_I1[-ni,]
oneVillage16000_pcase_ext_IR_new <-oneVillage16000_pcase_ext_IR[-ni,]

extInt_16000 <-read.csv(paste0(dir5,"extinction_interval_16000160001600016000reintRate_0.001000migRate_0.000000_seasonality_period_1yr_1pcase_inf_bet_case.csv"),header=FALSE,fill=TRUE)
ni<-which(is.na(extInt_16000[,1]))
fourVillage16000_pcase_I1_new <-fourVillage16000_pcase_I1[-ni,]
fourVillage16000_pcase_IR_new <-fourVillage16000_pcase_IR[-ni,]
fourVillage16000_pcase_ext_I1_new <-fourVillage16000_pcase_ext_I1[-ni,]
fourVillage16000_pcase_ext_IR_new <-fourVillage16000_pcase_ext_IR[-ni,]

extInt_16000 <-read.csv(paste0(dir5,"extinction_interval_16000160001600016000reintRate_0.001000migRate_0.100000_seasonality_period_1yr_1pcase_inf_bet_case.csv"),header=FALSE,fill=TRUE)
ni<-which(is.na(extInt_16000[,1]))
fourVillage16000_pcase_I1_new_mig_0.1 <-fourVillage16000_pcase_I1_mig_0.1[-ni,]
fourVillage16000_pcase_IR_new_mig_0.1 <-fourVillage16000_pcase_IR_mig_0.1[-ni,]
fourVillage16000_pcase_ext_I1_new_mig_0.1 <-fourVillage16000_pcase_ext_I1_mig_0.1[-ni,]
fourVillage16000_pcase_ext_IR_new_mig_0.1 <-fourVillage16000_pcase_ext_IR_mig_0.1[-ni,]



g<-rgeom(n=10000,p=(1/200))
d<-qqplot(t(twoVillage32000_pcase_I1_new),g)
d1<-qqplot(t(fourVillage16000_pcase_I1_new),g)
d2<-qqplot(t(eightVillage8000_pcase_I1_new),g)
d3<-qqplot(t(sixteenVillage4000_pcase_I1_new),g)
d4<-qqplot(t(thirtytwoVillage2000_pcase_I1_new),g)

marker=colorRampPalette(c('orange','red','purple','royalblue'))(7)

png(paste0(dir2,"qqplot_64000_metapop_i1_pcase_paper_v2.png"), width=1600, height=1600, res=180) 
par(fig=c(0,1,0,1),mar=c(3.3,3.6,1.1,1.1),oma=c(1,1,0,3),mgp=c(2.2,0.45,0),tcl=-0.4,cex=1.5)
qqplot((t(oneVillage64000_pcase_I1_new)),(g),xlab="Empirical quantiles from simulation",ylab="Quantiles from geometric distribution",pch=20,col='black',xlim=c(0,3000),ylim=c(0,3000),cex.axis=1.2,cex.lab=1.2)
abline(a=0,b=1,col='grey80',lwd=2)
points((d$x),(d$y),col=marker[1],pch=20)
points((d1$x),(d1$y),col=marker[2],pch=20)
points((d2$x),(d2$y),col=marker[4],pch=20)
points((d3$x),(d3$y),col=marker[6],pch=20)
points((d4$x),(d4$y),col=marker[7],pch=20)
legend('topleft',legend=rev(c('1x64k','2x32k','4x16k','8x8k','16x4k','32x2k')),col=rev(c('black',marker[1:2],marker[4],marker[6:7])),lwd=3,bty='n',cex=1.5,pt.cex=1)
dev.off()

g<-rgeom(n=10000,p=(1/200))
d<-qqplot(t(twoVillage32000_pcase_IR_new),g)
d1<-qqplot(t(fourVillage16000_pcase_IR_new),g)
d2<-qqplot(t(eightVillage8000_pcase_IR_new),g)
d3<-qqplot(t(sixteenVillage4000_pcase_IR_new),g)
d4<-qqplot(t(thirtytwoVillage2000_pcase_IR_new),g)

png(paste0(dir2,"qqplot_64000_metapop_ir_pcase_paper_v2.png"), width=1600, height=1600, res=180) 
par(fig=c(0,1,0,1),mar=c(3.3,3.6,1.1,1.1),oma=c(1,1,0,3),mgp=c(2.2,0.45,0),tcl=-0.4,cex=1.5)
qqplot((t(oneVillage64000_pcase_IR_new)),(g),xlab="Empirical quantiles from simulation",ylab="Quantiles from geometric distribution",pch=20,col='black',xlim=c(0,3000),ylim=c(0,3000),cex.axis=1.2,cex.lab=1.2)
abline(a=0,b=1,col='grey80',lwd=2)
points((d$x),(d$y),col=marker[1],pch=20)
points((d1$x),(d1$y),col=marker[2],pch=20)
points((d2$x),(d2$y),col=marker[4],pch=20)
points((d3$x),(d3$y),col=marker[6],pch=20)
points((d4$x),(d4$y),col=marker[7],pch=20)
legend('topleft',legend=rev(c('1x64k','2x32k','4x16k','8x8k','16x4k','32x2k')),col=rev(c('black',marker[1:2],marker[4],marker[6:7])),lwd=3,bty='n',cex=1.5,pt.cex=1)
dev.off()

g<-rgeom(n=10000,p=(1/200))
d<-qqplot(t(twoVillage32000_pcase_ext_I1_new),g)
d1<-qqplot(t(fourVillage16000_pcase_ext_I1_new),g)
d2<-qqplot(t(eightVillage8000_pcase_ext_I1_new),g)
d3<-qqplot(t(sixteenVillage4000_pcase_ext_I1_new),g)
d4<-qqplot(t(thirtytwoVillage2000_pcase_ext_I1_new),g)

png(paste0(dir2,"qqplot_64000_metapop_i1_pcase_ext_paper_v2.png"), width=1600, height=1600, res=180) 
par(fig=c(0,1,0,1),mar=c(3.3,3.6,1.1,1.1),oma=c(1,1,0,3),mgp=c(2.2,0.45,0),tcl=-0.4,cex=1.5)
qqplot((t(oneVillage64000_pcase_ext_I1_new)),(g),xlab="Empirical quantiles from simulation",ylab="Quantiles from geometric distribution",pch=20,col='black',xlim=c(0,3000),ylim=c(0,3000),cex.axis=1.2,cex.lab=1.2)
abline(a=0,b=1,col='grey80',lwd=2)
points((d$x),(d$y),col=marker[1],pch=20)
points((d1$x),(d1$y),col=marker[2],pch=20)
points((d2$x),(d2$y),col=marker[4],pch=20)
points((d3$x),(d3$y),col=marker[6],pch=20)
points((d4$x),(d4$y),col=marker[7],pch=20)
legend('topleft',legend=rev(c('1x64k','2x32k','4x16k','8x8k','16x4k','32x2k')),col=rev(c('black',marker[1:2],marker[4],marker[6:7])),lwd=3,bty='n',cex=1.5,pt.cex=1)
dev.off()

g<-rgeom(n=10000,p=(1/200))
d<-qqplot(t(twoVillage32000_pcase_ext_IR_new),g)
d1<-qqplot(t(fourVillage16000_pcase_ext_IR_new),g)
d2<-qqplot(t(eightVillage8000_pcase_ext_IR_new),g)
d3<-qqplot(t(sixteenVillage4000_pcase_ext_IR_new),g)
d4<-qqplot(t(thirtytwoVillage2000_pcase_ext_IR_new),g)

png(paste0(dir2,"qqplot_64000_metapop_ir_pcase_ext_paper_v2.png"), width=1600, height=1600, res=180) 
par(fig=c(0,1,0,1),mar=c(3.3,3.6,1.1,1.1),oma=c(1,1,0,3),mgp=c(2.2,0.45,0),tcl=-0.4,cex=1.5)
qqplot((t(oneVillage64000_pcase_ext_IR_new)),(g),xlab="Empirical quantiles from simulation",ylab="Quantiles from geometric distribution",pch=20,col='black',xlim=c(0,3000),ylim=c(0,3000),cex.axis=1.2,cex.lab=1.2)
abline(a=0,b=1,col='grey80',lwd=2)
points((d$x),(d$y),col=marker[1],pch=20)
points((d1$x),(d1$y),col=marker[2],pch=20)
points((d2$x),(d2$y),col=marker[4],pch=20)
points((d3$x),(d3$y),col=marker[6],pch=20)
points((d4$x),(d4$y),col=marker[7],pch=20)
legend('topleft',legend=rev(c('1x64k','2x32k','4x16k','8x8k','16x4k','32x2k')),col=rev(c('black',marker[1:2],marker[4],marker[6:7])),lwd=3,bty='n',cex=1.5,pt.cex=1)
dev.off()




png(paste0(dir2,"qqplot_64000_ir_pcase_paper.png"), width=1600, height=1600, res=180) 
qqplot(t(oneVillage64000_pcase_IR_new),g,xlab="Simulated",ylab="Theoretical")
abline(a=0,b=1,col='red',lwd=2)
dev.off()

g<-rgeom(n=nrow(oneVillage64000_pcase_ext_I1_new),p=(1/200))
png(paste0(dir2,"qqplot_64000_i1_pcase_ext_paper.png"), width=1600, height=1600, res=180) 
qqplot(t(oneVillage64000_pcase_ext_I1_new),g,xlab="Simulated",ylab="Theoretical")
abline(a=0,b=1,col='red',lwd=2)
dev.off()

png(paste0(dir2,"qqplot_64000_ir_pcase_ext_paper.png"), width=1600, height=1600, res=180) 
qqplot(t(oneVillage64000_pcase_ext_IR_new),g,xlab="Simulated",ylab="Theoretical")
abline(a=0,b=1,col='red',lwd=2)
dev.off()

g<-rgeom(n=nrow(oneVillage32000_pcase_I1_new),p=(1/200))
png(paste0(dir2,"qqplot_32000_i1_pcase_paper.png"), width=1600, height=1600, res=180) 
par(mfrow=c(1,1))
qqplot(t(oneVillage32000_pcase_I1_new),g,xlab="Simulated",ylab="Theoretical")
abline(a=0,b=1,col='red',lwd=2)
dev.off()

png(paste0(dir2,"qqplot_32000_ir_pcase_paper.png"), width=1600, height=1600, res=180) 
qqplot(t(oneVillage32000_pcase_IR_new),g,xlab="Simulated",ylab="Theoretical")
abline(a=0,b=1,col='red',lwd=2)
dev.off()

g<-rgeom(n=nrow(oneVillage32000_pcase_ext_I1_new),p=(1/200))
png(paste0(dir2,"qqplot_32000_i1_pcase_ext_paper.png"), width=1600, height=1600, res=180) 
qqplot(t(oneVillage32000_pcase_ext_I1_new),g,xlab="Simulated",ylab="Theoretical")
abline(a=0,b=1,col='red',lwd=2)
dev.off()

png(paste0(dir2,"qqplot_32000_ir_pcase_ext_paper.png"), width=1600, height=1600, res=180) 
qqplot(t(oneVillage32000_pcase_ext_IR_new),g,xlab="Simulated",ylab="Theoretical")
abline(a=0,b=1,col='red',lwd=2)
dev.off()

g<-rgeom(n=nrow(oneVillage16000_pcase_I1_new),p=(1/200))
png(paste0(dir2,"qqplot_16000_i1_pcase_paper.png"), width=1600, height=1600, res=180) 
par(mfrow=c(1,1))
qqplot(t(oneVillage16000_pcase_I1_new),g,xlab="Simulated",ylab="Theoretical")
abline(a=0,b=1,col='red',lwd=2)
dev.off()

g<-rgeom(n=nrow(oneVillage16000_pcase_I1_new_mig01),p=(1/200))
png(paste0(dir2,"qqplot_16000_mig01_i1_pcase_paper.png"), width=1600, height=1600, res=180) 
par(mfrow=c(1,1))
qqplot(t(oneVillage16000_pcase_I1_new_mig01),g,xlab="Simulated",ylab="Theoretical")
abline(a=0,b=1,col='red',lwd=2)
dev.off()

g<-rgeom(n=nrow(oneVillage16000_pcase_I1_new_mig0),p=(1/200))
png(paste0(dir2,"qqplot_16000_mig0_i1_pcase_paper.png"), width=1600, height=1600, res=180) 
par(mfrow=c(1,1))
qqplot(t(oneVillage16000_pcase_I1_new_mig0),g,xlab="Simulated",ylab="Theoretical")
abline(a=0,b=1,col='red',lwd=2)
dev.off()



png(paste0(dir2,"qqplot_16000_ir_pcase_paper.png"), width=1600, height=1600, res=180) 
qqplot(t(oneVillage16000_pcase_IR_new),g,xlab="Simulated",ylab="Theoretical")
abline(a=0,b=1,col='red',lwd=2)
dev.off()

g<-rgeom(n=nrow(oneVillage16000_pcase_ext_I1_new),p=(1/200))
png(paste0(dir2,"qqplot_16000_i1_pcase_ext_paper.png"), width=1600, height=1600, res=180) 
qqplot(t(oneVillage16000_pcase_ext_I1_new),g,xlab="Simulated",ylab="Theoretical")
abline(a=0,b=1,col='red',lwd=2)
dev.off()

png(paste0(dir2,"qqplot_16000_ir_pcase_ext_paper.png"), width=1600, height=1600, res=180) 
qqplot(t(oneVillage16000_pcase_ext_IR_new),g,xlab="Simulated",ylab="Theoretical")
abline(a=0,b=1,col='red',lwd=2)
dev.off()

g<-rgeom(n=nrow(oneVillage8000_pcase_I1_new),p=(1/200))
png(paste0(dir2,"qqplot_8000_i1_pcase_paper.png"), width=1600, height=1600, res=180) 
par(mfrow=c(1,1))
qqplot(t(oneVillage8000_pcase_I1_new),g,xlab="Simulated",ylab="Theoretical")
abline(a=0,b=1,col='red',lwd=2)
dev.off()

png(paste0(dir2,"qqplot_8000_ir_pcase_paper.png"), width=1600, height=1600, res=180) 
qqplot(t(oneVillage8000_pcase_IR_new),g,xlab="Simulated",ylab="Theoretical")
abline(a=0,b=1,col='red',lwd=2)
dev.off()

g<-rgeom(n=nrow(oneVillage8000_pcase_ext_I1_new),p=(1/200))
png(paste0(dir2,"qqplot_8000_i1_pcase_ext_paper.png"), width=1600, height=1600, res=180) 
qqplot(t(oneVillage8000_pcase_ext_I1_new),g,xlab="Simulated",ylab="Theoretical")
abline(a=0,b=1,col='red',lwd=2)
dev.off()

png(paste0(dir2,"qqplot_8000_ir_pcase_ext_paper.png"), width=1600, height=1600, res=180) 
qqplot(t(oneVillage8000_pcase_ext_IR_new),g,xlab="Simulated",ylab="Theoretical")
abline(a=0,b=1,col='red',lwd=2)
dev.off()

g<-rgeom(n=nrow(oneVillage4000_pcase_I1_new),p=(1/200))
png(paste0(dir2,"qqplot_4000_i1_pcase_paper.png"), width=1600, height=1600, res=180) 
par(mfrow=c(1,1))
qqplot(t(oneVillage4000_pcase_I1_new),g,xlab="Simulated",ylab="Theoretical")
abline(a=0,b=1,col='red',lwd=2)
dev.off()

g<-rgeom(n=nrow(oneVillage4000_pcase_I1_new_mig0),p=(1/200))
png(paste0(dir2,"qqplot_4000_mig0_i1_pcase_paper.png"), width=1600, height=1600, res=180) 
par(mfrow=c(1,1))
qqplot(t(oneVillage4000_pcase_I1_new_mig0),g,xlab="Simulated",ylab="Theoretical")
abline(a=0,b=1,col='red',lwd=2)
dev.off()

g<-rgeom(n=nrow(oneVillage4000_pcase_I1_new_mig01),p=(1/200))
png(paste0(dir2,"qqplot_4000_mig01_i1_pcase_paper.png"), width=1600, height=1600, res=180) 
par(mfrow=c(1,1))
qqplot(t(oneVillage4000_pcase_I1_new_mig01),g,xlab="Simulated",ylab="Theoretical")
abline(a=0,b=1,col='red',lwd=2)
dev.off()


png(paste0(dir2,"qqplot_4000_ir_pcase_paper.png"), width=1600, height=1600, res=180) 
qqplot(t(oneVillage4000_pcase_IR_new),g,xlab="Simulated",ylab="Theoretical")
abline(a=0,b=1,col='red',lwd=2)
dev.off()

g<-rgeom(n=nrow(oneVillage4000_pcase_ext_I1_new),p=(1/200))
png(paste0(dir2,"qqplot_4000_i1_pcase_ext_paper.png"), width=1600, height=1600, res=180) 
qqplot(t(oneVillage4000_pcase_ext_I1_new),g,xlab="Simulated",ylab="Theoretical")
abline(a=0,b=1,col='red',lwd=2)
dev.off()

png(paste0(dir2,"qqplot_4000_ir_pcase_ext_paper.png"), width=1600, height=1600, res=180) 
qqplot(t(oneVillage4000_pcase_ext_IR_new),g,xlab="Simulated",ylab="Theoretical")
abline(a=0,b=1,col='red',lwd=2)
dev.off()



png(paste0(dir2,"qqplot_4x16000_i1_pcase.png"), width=1600, height=1600, res=180) 
qqplot(t(fourVillage16000_pcase_I1_new)[1:1000,],g,main="4x16000 - Number of I1 infections between pcases",xlab="Simulated",ylab="rgeom(10000,1/200)")
abline(a=0,b=1,col='red',lwd=2)
dev.off()

png(paste0(dir2,"qqplot_4x16000_ir_pcase.png"), width=1600, height=1600, res=180) 
qqplot(t(fourVillage16000_pcase_IR_new)[1:1000,],g,main="4x16000 - Number of IR infections between pcases",xlab="Simulated",ylab="rgeom(10000,1/200)")
abline(a=0,b=1,col='red',lwd=2)
dev.off()

png(paste0(dir2,"qqplot_4x16000_i1_pcase_ext.png"), width=1600, height=1600, res=180) 
qqplot(t(fourVillage16000_pcase_ext_I1_new)[1:1000,],g,main="4x16000 - Number of I1 infections between pcase and extinction",xlab="Simulated",ylab="rgeom(10000,1/200)")
abline(a=0,b=1,col='red',lwd=2)
dev.off()

png(paste0(dir2,"qqplot_4x16000_ir_pcase_ext.png"), width=1600, height=1600, res=180) 
qqplot(t(fourVillage16000_pcase_ext_IR_new)[1:1000,],g,main="4x16000 - Number of IR infections between pcase and extinction",xlab="Simulated",ylab="rgeom(10000,1/200)")
abline(a=0,b=1,col='red',lwd=2)
dev.off()

png(paste0(dir2,"qqplot_4x16000_mig_0.1_i1_pcase.png"), width=1600, height=1600, res=180) 
qqplot(t(fourVillage16000_pcase_I1_new_mig_0.1)[1:1000,],g,main="4x16000 mig 0.1 - Number of I1 infections between pcases",xlab="Simulated",ylab="rgeom(10000,1/200)")
abline(a=0,b=1,col='red',lwd=2)
dev.off()

png(paste0(dir2,"qqplot_4x16000_mig_0.1_ir_pcase.png"), width=1600, height=1600, res=180) 
qqplot(t(fourVillage16000_pcase_IR_new_mig_0.1)[1:1000,],g,main="4x16000 mig 0.1 - Number of IR infections between pcases",xlab="Simulated",ylab="rgeom(10000,1/200)")
abline(a=0,b=1,col='red',lwd=2)
dev.off()

png(paste0(dir2,"qqplot_4x16000_mig_0.1_i1_pcase_ext.png"), width=1600, height=1600, res=180) 
qqplot(t(fourVillage16000_pcase_ext_I1_new_mig_0.1)[1:1000,],g,main="4x16000 mig 0.1 - Number of I1 infections between pcase and extinction",xlab="Simulated",ylab="rgeom(10000,1/200)")
abline(a=0,b=1,col='red',lwd=2)
dev.off()

png(paste0(dir2,"qqplot_4x16000_mig_0.1_ir_pcase_ext.png"), width=1600, height=1600, res=180) 
qqplot(t(fourVillage16000_pcase_ext_IR_new_mig_0.1)[1:1000,],g,main="4x16000 mig 0.1 - Number of IR infections between pcase and extinction",xlab="Simulated",ylab="rgeom(10000,1/200)")
abline(a=0,b=1,col='red',lwd=2)
dev.off()




#calculate average time between infections (no seasonality)
beta = 135
n = 16000
gamma = 13
kappa = 0.4179
pir = 0.005
seasonalAmp = 0
s = 927 #16000
i1 = 23 #16000
ir = 37 #16000
#s = 1#3707 #64000
#ir = 63999#93 #64000
#i1 = 0#147 #64000
#seasonalBeta = beta*(1+seasonalAmp*sin((2*pi)*time))
foi = beta*(i1+kappa*ir)/n
firstInf = s*foi
timeBet = 1/firstInf
timeTo200 = 250*timeBet

p = as.vector(apply(oneVillage16000_pcase_I1_new,1, function(x) length(which(!is.na(x)))))
p1 = as.vector(apply(oneVillage64000_pcase_I1,1, function(x) length(which(!is.na(x)))))



e<-rcomppois(10000,1/mean(p))

qqplot(p,e)
abline(a=0,b=1,col='red',lwd=2)




