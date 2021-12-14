#I1 = read.table('N=10000,beta=135,fast_I1vec_time=15_birth=.02_valid.csv', header=F, fill=T, col.names = paste0("V", seq_len(1000)))
#Ir = read.table('N=10000,beta=135,fast_Irvec_time=15_birth=.02_valid.csv', header=F, fill=T, col.names = paste0("V", seq_len(1000)))
#S = read.table('N=10000,beta=135,fast_Svec_time=15_birth=.02_valid.csv', header=F, fill=T, col.names = paste0("V", seq_len(1000)))
#S = read.table('N=10000,beta=135,fast_Svec_test_time=15_birth=.02_valid.csv', header=F, fill=T, col.names = paste0("V", seq_len(1000)))
#R = read.table('N=10000,beta=135,fast_Rvec_test_time=15_birth=.02_valid.csv', header=F, fill=T, col.names = paste0("V", seq_len(1000)))
#P = read.table('N=10000,beta=135,fast_Pvec_test_time=15_birth=.02_valid.csv', header=F, fill=T, col.names = paste0("V", seq_len(1000)))

plot_class = function(d, label) {
    d       = d[,1:800]
    means   = colMeans(d,na.rm=T)
    ses     = apply(d,MARGIN=2,FUN=sd,na.rm=T)/sqrt(dim(d)[2])
    top     = means + 1.96*ses
    bottom  = means - 1.96*ses
    #print(c(min(bottom,na.rm=T),max(top,na.rm=F)))
    x       = (0:(length(means)-1))*3+1
    plot(x, means, ylab=label, xlab='', type='n', ylim=c(min(bottom,na.rm=T),max(top,na.rm=T)))
    polygon(c(x,rev(x)), c(bottom,rev(top)), col='#88888888',border='NA')
    lines(x, means)
}


png('gillespie_test2.png',width=1000, height=2000, res=250)
par(mfrow=c(5,1), mar=c(2.1,4.1,1,1), oma=c(2,0,0,0), las=1)
plot_class(I1, 'I1')
plot_class(Ir, 'Ir')
plot_class(S,  'S')
plot_class(R,  'R')
plot_class(P,  'P')
mtext('Step', side=1, outer=T, cex=0.8, line=0.5)
dev.off()
