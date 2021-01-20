# Set working directory
setwd("~/Documents/UCLA/Winter 2021/biostat234/lab 3/")
getwd()

#LOAD NECESSARY PACKAGES
library(R2jags)
library(lattice)
#Function to help with the burn in 
load("AddBurnin.RData")


# useful function
mysummary = function(invector) {
c(mean(invector), sd(invector), quantile(invector, .025), 
	quantile(invector,.975),
	length(invector[invector>0])/length(invector))
}


# Save the data file TraumaData.txt to your working directory. 

#load in the data 

TraumaData = read.table("bcj97data.txt")

#Give the columns useful names 
colnames(TraumaData) <- c("death", "n",	"intercept",	"iss",	"rts", "age", "ti", "age*ti")

#first 6 observations from proto-typical cases 
head(TraumaData, n=6)

#First we divide up the dataset into the 6 proto-typical cases and the observed data

#For the 6 proto-typical cases define Xp the design matrix, Yp the outcomes (death=1), and np the number of trials
Xp = as.matrix(TraumaData[1:6,3:8])
Yp = TraumaData[1:6,1]
np = TraumaData[1:6,2]

#note a=Yp+1 and b=np-a+2 in the corresponding beta distributions 

#define the inverse of the Xp matrix to be used to convert the prior distributions on pi to distributions on the regression coefficients
invXp = solve(Xp)

#For the observed data define the design matrix, outcomes, and number of trials

Xobs = as.matrix(TraumaData[7:306,3:8])
Yobs = TraumaData[7:306,1]
nobs = TraumaData[7:306,2]

#To see how prior distributions on the betas we define a simple model which maps the  
#distributions on pi to distributions


#Store the model in the file LAB3.Priors.txt

sink("LAB3.Priors.txt")
cat("
model{
	
	betas<-invXp %*% logitp[]

	for(j in 1:6){
		logitp[j]<-logit(pie[j])
	}
	pie[1]~dbeta(1.1,8.5)
	pie[2]~dbeta(3.0,11.0)
	pie[3]~dbeta(5.9,1.7)
	pie[4]~dbeta(1.3,12.9)
	pie[5]~dbeta(1.1,4.9)
	pie[6]~dbeta(1.5,5.5)

}

  ",fill = TRUE)
sink()

#In the code we write p (or pie) for pi since pi already has a meaning in R

ex1.data=list(invXp=invXp)
ex1.inits=rep(list(list(pie=c(0.5,0.5,0.5,0.5,0.5,0.5))),5)
ex1.parameters = c("betas", "pie[1:6]")

#Run the JAGS model
set.seed(1234)
ex1.out = jags(ex1.data, ex1.inits, ex1.parameters, "LAB3.Priors.txt", 
	n.chains=5, n.iter=11000, n.burnin=0, n.thin=2, DIC=F)

names(ex1.out) #components of ex1.out

#Treat the first 1000 iterations as a burn in	
Output1 = AddBurnin(ex1.out$BUGSoutput$sims.array,burnin=1000,n.thin=2)


print(Output1$Burnin.Summary)

#Now we incorporate the data into the model

sink("Lab3.Posteriors.txt")
cat("
model{
	
	betas<-invXp %*% logitp[]

	for(j in 1:6){
		logitp[j]<-logit(pie[j])
	}
	pie[1]~dbeta(1.1,8.5)
	pie[2]~dbeta(3.0,11.0)
	pie[3]~dbeta(5.9,1.7)
	pie[4]~dbeta(1.3,12.9)
	pie[5]~dbeta(1.1,4.9)
	pie[6]~dbeta(1.5,5.5)

	
		for(i in 1:T){
		y[i] ~ dbern(p[i])
		p[i]<-ilogit(inprod(x[i,],betas[]))
		}
			
	
}
  ",fill = TRUE)
sink()

#fit the model

ex2.data = list(x=Xobs, y=Yobs, T=300, invXp=invXp)
ex2.inits = rep(list(list(pie=c(0.5,0.5,0.5,0.5,0.5,0.5))),5)
ex2.parameters = c("betas", "pie[1:6]")

set.seed(1234)
ex2.out = jags(ex2.data, ex2.inits, ex2.parameters, "Lab3.Posteriors.txt", 
	n.chains=5, n.iter=11000, n.burnin=0, n.thin=2, DIC=F)

names(ex2.out)

#Treat the first 1000 iterations as a burn in	
Output2 = AddBurnin(ex2.out$BUGSoutput$sims.array,burnin=1000,n.thin=2)

names(Output2)

print(Output2$Burnin.Summary)

#Graphs!

temp3=Output1$Burnin.sims.matrix
temp4=Output2$Burnin.sims.matrix
par(mfrow=c(2,3))
plot(density(temp4[,7]),xlim=c(0,1),main="",xlab="pi1")  
# pi1.  
lines(density(temp3[,7],bw=.055),lty=2,lwd=2, col="blue")  
# prior
legend(.3, 8, legend=c("post","prior"), col=c("black", "blue"), lty=1:2 , lwd=c(1,2))
plot(density(temp4[,8]),xlim=c(0,1),main="",xlab="pi2")  
# pi2.  
lines(density(temp3[,8],bw=.05),lty=2,lwd=2, col="blue")  
# prior
plot(density(temp4[,9]),xlim=c(0,1),main="",xlab="pi3")  
# pi3.  
lines(density(temp3[,9],bw=.075),lty=2,lwd=2, col="blue")  
# prior
plot(density(temp4[,10]),xlim=c(0,1),main="",xlab="pi4")  
# pi4.  
lines(density(temp3[,10],bw=.045),lty=2,lwd=2, col="blue")  
# prior
plot(density(temp4[,11]),xlim=c(0,1),main="",xlab="pi5")  
# pi5.  
lines(density(temp3[,11],bw=.085),lty=2,lwd=2, col="blue")  
# prior
plot(density(temp4[,12]),xlim=c(0,1),main="",xlab="pi6")  
# pi6.  
lines(density(temp3[,12],bw=.085),lty=2,lwd=2, col="blue")  
# prior


par(mfrow=c(2,3))
acf(ex2.out$BUGSoutput$sims.array[1:5000,1,1], main="beta1")
acf(ex2.out$BUGSoutput$sims.array[1:5000,1,2], main="beta2")
acf(ex2.out$BUGSoutput$sims.array[1:5000,1,3], main="beta3")
acf(ex2.out$BUGSoutput$sims.array[1:5000,1,4], main="beta4")
acf(ex2.out$BUGSoutput$sims.array[1:5000,1,5], main="beta5")
acf(ex2.out$BUGSoutput$sims.array[1:5000,1,6], main="beta6")

par(mfrow=c(2,3))
acf(ex2.out$BUGSoutput$sims.array[1:5000,1,7], main="pi1")
acf(ex2.out$BUGSoutput$sims.array[1:5000,1,8], main="pi2")
acf(ex2.out$BUGSoutput$sims.array[1:5000,1,9], main="pi3")
acf(ex2.out$BUGSoutput$sims.array[1:5000,1,10], main="pi4")
acf(ex2.out$BUGSoutput$sims.array[1:5000,1,11], main="pi5")
acf(ex2.out$BUGSoutput$sims.array[1:5000,1,12], main="pi6")

par(mfrow=c(2,3))
plot(ex2.out$BUGSoutput$sims.array[4500:5000,1,1], ylab="beta1", type="l", main="chain 1")
plot(ex2.out$BUGSoutput$sims.array[4500:5000,2,1], ylab="beta1", type="l", col="blue", main="chain 2")
plot(ex2.out$BUGSoutput$sims.array[4500:5000,3,1], ylab="beta1", type="l", col="steelblue", main="chain 3")
plot(ex2.out$BUGSoutput$sims.array[4500:5000,4,1], ylab="beta1", type="l", col="aquamarine", main="chain 4")
plot(ex2.out$BUGSoutput$sims.array[4500:5000,5,1], ylab="beta1", type="l", col="cyan", main="chain 5")
plot(ex2.out$BUGSoutput$sims.array[4500:5000,1,1], ylab="beta1", type="l", main="All chains")
lines(ex2.out$BUGSoutput$sims.array[4500:5000,2,1], type="l", col="blue")
lines(ex2.out$BUGSoutput$sims.array[4500:5000,3,1], type="l", col="steelblue")
lines(ex2.out$BUGSoutput$sims.array[4500:5000,4,1], type="l", col="aquamarine")
lines(ex2.out$BUGSoutput$sims.array[4500:5000,5,1], type="l", col="cyan")


par(mfrow=c(2,3))
plot(ex2.out$BUGSoutput$sims.array[4500:5000,1,1], ylab="beta1", type="l")
plot(ex2.out$BUGSoutput$sims.array[4500:5000,1,2], ylab="beta2", type="l", col="blue", main="chain 1")
plot(ex2.out$BUGSoutput$sims.array[4500:5000,1,3], ylab="beta3", type="l", col="steelblue")
plot(ex2.out$BUGSoutput$sims.array[4500:5000,1,4], ylab="beta4", type="l", col="aquamarine")
plot(ex2.out$BUGSoutput$sims.array[4500:5000,1,5], ylab="beta5", type="l", col="cyan")
plot(ex2.out$BUGSoutput$sims.array[4500:5000,1,6], ylab="beta6", type="l", col="skyblue")

par(mfrow=c(2,3))
plot(ex2.out$BUGSoutput$sims.array[,1,1], ylab="beta1", type="l")
plot(ex2.out$BUGSoutput$sims.array[,1,2], ylab="beta2", type="l", col="blue", main="chain 1")
plot(ex2.out$BUGSoutput$sims.array[,1,3], ylab="beta3", type="l", col="steelblue")
plot(ex2.out$BUGSoutput$sims.array[,1,4], ylab="beta4", type="l", col="aquamarine")
plot(ex2.out$BUGSoutput$sims.array[,1,5], ylab="beta5", type="l", col="cyan")
plot(ex2.out$BUGSoutput$sims.array[,1,6], ylab="beta6", type="l", col="skyblue")

