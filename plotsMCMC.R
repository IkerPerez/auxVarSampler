
Output1 <- read.csv("OutputReal.csv", stringsAsFactors=FALSE)
Output1 <- Output1[,1:12]
summary(Output1)
reals <- colMeans(Output1)

par(mfrow=c(2,2),oma=c(0,0,2,0),mar=c(5,5,1,3))
for(i in 1:12){
  ts.plot(Output1[501:1000,i],xlab="iteration",ylab=" ",main=" ")
  plot(density(Output1[501:1000,i]),main=" ")
  acf(Output1[501:1000,i],main=" ")
  pacf(Output1[501:1000,i],main=" ")
  title(colnames(Output1)[i], outer=TRUE)
}

Output1 <- read.csv("Output2cp07m05.csv", stringsAsFactors=FALSE)
Output2 <- read.csv("Output15cp05m05.csv", stringsAsFactors=FALSE)
Output3 <- read.csv("Output12cp02m05.csv", stringsAsFactors=FALSE)
Output1 <- Output1[,1:12]; Output2 <- Output2[,1:12]; Output3 <- Output3[,1:12]
summary(Output1); summary(Output2); summary(Output3)


round(cor(Output1[1:100000,]),2); # Massive correlation as expected!
round(cor(Output2[1:100000,]),2); # Massive correlation as expected!
round(cor(Output3[1:100000,]),3); # Massive correlation as expected!

# Thinning and burn-in...
burnin = 10000
thinning = 10
Output1New <- Output1[seq(burnin+1,nrow(Output1),thinning),]; summary(Output1New); 
Output2New <- Output2[seq(burnin+1,nrow(Output2),thinning),]; summary(Output2New); 
Output3New <- Output3[seq(burnin+1,nrow(Output3),thinning),]; summary(Output3New); 
par(mfrow=c(2,2),oma=c(1,0,0,0),mar=c(2.5,5,5,3))
for(i in 4:12){
  acf(Output1New[,i],ylim=c(0,1),main=colnames(Output1)[i],col="blue");
  lines(seq(0,39)+0.3,acf(Output2New[,i],ylim=c(0,1),plot = F)[[1]],type="h",col="red")
  lines(seq(0,39)+0.6,acf(Output3New[,i],ylim=c(0,1),plot = F)[[1]],type="h",col="orange")
  pacf(Output1New[,i],ylim=c(0,1),main=colnames(Output1)[i],col="blue"); 
  lines(seq(1,39)+0.3,pacf(Output2New[,i],ylim=c(0,1),plot = F)[[1]],type="h",col="red")
  lines(seq(1,39)+0.6,pacf(Output3New[,i],ylim=c(0,1),plot = F)[[1]],type="h",col="orange")
}
#Save EPSat 1000x1414

par(mfrow=c(2,2),oma=c(1,0,0,0),mar=c(2.5,5,5,3))
for(i in 4:12){
  x <- Output1New[,i]
  y <- Output2New[,i]
  z <- Output3New[,i]
  ts.plot(Output1New[,i],col="blue",xlab="iteration",ylab=" ",main=colnames(Output1)[i],ylim=c(min(c(x,y,z)),max(c(x,y,z))))
  lines(1:nrow(Output2New),Output2New[,i],col="red")
  lines(1:nrow(Output3New),Output3New[,i],col="orange")
  abline(reals[i],0,lwd=2,lty=2)
  x <- density(Output1New[,i])
  y <- density(Output2New[,i])
  z <- density(Output3New[,i])
  plot(density(Output1New[,i]),main=colnames(Output1)[i],col="blue",xlim=c(min(c(x$x,y$x,z$x)),max(c(x$x,y$x,z$x))))
  lines(c(reals[i],reals[i]),c(0,100),lwd=1,lty=2)
  lines(density(Output2New[,i]),col="red")
  lines(density(Output3New[,i]),col="orange")
  #title(colnames(Output1)[i], outer=TRUE)
}
#Save EPSat 1000x1414


library(coda)  
combinedchains = mcmc.list(as.mcmc(Output1New), as.mcmc(Output2New), as.mcmc(Output3New))  
plot(combinedchains)  
gelman.diag(combinedchains)
gelman.plot(combinedchains)  

summary(combinedchains)
crosscorr(combinedchains)

