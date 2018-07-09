# real data
d<-read.csv("./k5.csv", sep=",")
summary(d)

#simulated data
q<-read.csv("./k5Model.csv", sep=",")
summary(q)

#simulated data
u<-read.csv("./k5Unif.csv", sep=",")
summary(u)


#if not too many NAs
d<-d[complete.cases(d),]

par(mfrow=c(1,2))
#plot(table(as.numeric(as.character(rnbinom(100,mu=44,size=400)))))
plot(table(as.numeric(as.character(d[d$Steps==0,]$Length))),col='red')
plot(table(as.numeric(as.character(d[d$Steps==7,]$Length))),col='red')



points(table(as.numeric(as.character(52-rpois(1000,lambda=4)))))

                                        #comparison

# length plot, need to reset graphical parameters depending on number of experiments
rc<-54
iters<-dim(q)[1]/(max(q$Steps)+1)
plot.new()
par(mfrow=c(4,4))
for (i in 0:7) {
    model<-as.numeric(as.character(q[q$Steps==i,]$Length))
    modelAug<-table(c(model,0:rc))
    unif<-as.numeric(as.character(u[u$Steps==i,]$Length))
    unifAug<-table(c(unif,0:rc))
    real<-as.numeric(as.character(d[d$Steps==i,]$Length))
    realAug<-table(c(real,0:rc))
    plot(table(model),xlab="Length",ylab="Frequency (100 samples)",main=paste(toString(i)," steps, TV = ", toString(sum(abs(modelAug-unifAug))/(2*iters)), sep=""),xlim=c(min(q[q$Steps==i,]$Length),rc),col='red');
    points(table(unif))
    plot(table(real),xlab="Length",ylab="Frequency (100 samples)",main=paste(toString(i)," steps, TV = ", toString(sum(abs(realAug-unifAug))/(2*iters)), sep=""),xlim=c(min(q[q$Steps==i,]$Length),rc),col='blue');
    points(table(unif))
}

dev.copy(png,"5buscompletelength.png")
dev.off()

par(mfrow=c(1,2))
hist(u[u$Steps==0,]$AvgLen)
hist(d[d$Steps==0,]$AvgLen)

rc<-54
iters<-dim(q)[1]/(max(q$Steps)+1)
plot.new()
par(mfrow=c(4,4))
for (i in 0:7) {
    model<-as.numeric(as.character(q[q$Steps==i,]$AvgLen))
    modelAug<-table(c(model,0:rc))
    unif<-as.numeric(as.character(u[u$Steps==i,]$AvgLen))
    unifAug<-table(c(unif,0:rc))
    real<-as.numeric(as.character(d[d$Steps==i,]$AvgLen))
    realAug<-table(c(real,0:rc))
    plot(table(model),xlab="AvgLen",ylab="Frequency (100 samples)",main=paste("Model length w #steps = ",toString(i)," , sep=""),xlim=c(min(q[q$Steps==i,]$AvgLen),rc),col='red');
    points(table(unif))
    plot(table(real),xlab="AvgLen",ylab="Frequency (100 samples)",main=paste("Two-edge loop length w #steps = ",toString(i) ,  sep=""),xlim=c(min(q[q$Steps==i,]$AvgLen),rc),col='blue');
    points(table(unif))
}
