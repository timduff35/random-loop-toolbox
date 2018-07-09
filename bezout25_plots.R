# real data
d<-read.csv("./bezout2vars5deg1000iters8steps.csv", sep=",")
summary(d)

#simulated data
q<-read.csv("./bezoutModel.csv", sep=",")

summary(q)

#simulated data
u<-read.csv("./bezoutUnif.csv", sep=",")
summary(u)


#if not too many NAs
d<-d[complete.cases(d),]

plot(table(as.numeric(as.character(d[d$Steps==0,]$Length))),col='red')
points(table(as.numeric(as.character(25-rpois(1000,lambda=8.31)))))


#h<-hist(rnorm(1000,mean = 25-log(25),sd=sqrt(log(25))),plot=FALSE)


plot(table(as.numeric(as.character(54-rpois(1000,lambda=11)))))

plot.new()
par(mfrow=c(4,4))
for (i in 0:7) {
    model<-as.numeric(as.character(q[q$Steps==i,]$Length))
    modelAug<-table(c(model,0:25))
    unif<-as.numeric(as.character(u[u$Steps==i,]$Length))
    unifAug<-table(c(unif,0:25))
    real<-as.numeric(as.character(d[d$Steps==i,]$Length))
    realAug<-table(c(real,0:25))
    plot(table(model),xlab="Length",ylab="Frequency (1000 samples)",main=paste(toString(i)," steps, TV = ", toString(sum(abs(modelAug-unifAug))/2000), sep=""),xlim=c(min(q[q$Steps==i,]$Length),28));
    points(table(unif),col='red')
    plot(table(real),xlab="Length",ylab="Frequency (1000 samples)",main=paste(toString(i)," steps, TV = ", toString(sum(abs(realAug-unifAug))/2000), sep=""),xlim=c(min(q[q$Steps==i,]$Length),28));
    points(table(unif),col='blue')
}

dev.copy(png,"bezout25.png")
dev.off()
    
for (i in 0:7) {                                        
