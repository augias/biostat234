
run2 = jags(priordata, inits, parameters, "lab4model.txt", 
	n.chains=5, n.iter=5100, n.burnin=0, n.thin=1, DIC=F)


names(run2)

Output2=AddBurnin(run2$BUGSoutput$sims.array,burnin=100,n.thin=1)

print(Output2$Burnin.Summary)

temp = Output2$Burnin.sims.matrix

dim(temp)
round(t(apply(temp,2,mysummary)),3)  # the t() function is the transpose
                              # and round rounds to 3 decimal places.  

###
### posterior densities
###

#you may need to set the limits of the x-axes with xlim=c(,) to get clearer graphs

par(mfrow=c(3,3))
par(mar=c(3.1,4.1,2.1,2.1))

plot(density(temp[,1]), main="alpha[1]", xlim=c(2.6,3.6))
plot(density(temp[,2]), main="alpha[2]", xlim=c(-0.3,1.0))
abline(v=0)
plot(density(temp[,3]), main="alpha[3]", xlim=c(-0.6,0.8))
abline(v=0)
plot(density(temp[,4]), main="alpha[4]",  xlim=c(-0.6,0.6))
abline(v=0)
plot(density(temp[,5]), main="alpha[5]",  xlim=c(-0.8,0.6))
abline(v=0)
plot(density(temp[,6]), main="alpha[6]",  xlim=c(-0.8,0.5))
abline(v=0)
plot(density(temp[,7]), main="alpha[7]",  xlim=c(-0.3,1))
abline(v=0)
plot(density(temp[,8]), main="alpha[8]",  xlim=c(-0.9,0.4))
abline(v=0)
hist(temp[,9], main="alphasign[1]")


#dev.off()   # closes the last graphics window.  Don't
		# execute if you want to keep the last plot.
		
par(mfrow=c(3,3))
par(mar=c(3.1,4.1,2.1,2.1))		
hist(temp[,10], main="alphasign[2]")
hist(temp[,11], main="alphasign[3]")
hist(temp[,12], main="alphasign[4]")
hist(temp[,13], main="alphasign[5]")
hist(temp[,14], main="alphasign[6]")
hist(temp[,15], main="alphasign[7]")
hist(temp[,16], main="alphasign[8]")
plot(density(temp[,26]), main="tau.e",   xlim=c(4,9))
plot(density(temp[,25]), main="tau.b",   xlim=c(1,5))

#dev.off() 
par(mfrow=c(3,3))
par(mar=c(3.1,4.1,2.1,2.1))
plot(density(temp[,26]), main="tau.e",  xlim=c(3.5,9))
plot(density(temp[,25]), main="tau.b",  xlim=c(1,5))
plot(density(temp[,23]), main="sigma",  xlim=c(.32,0.5))
plot(density(temp[,24]), main="sqrtD",  xlim=c(0.4,0.9))
plot(density(temp[,22]), main="rho",  xlim=c(0.52,0.85))
plot(density(temp[,17]), main="beta[1]",  xlim=c(-1.2,0.6))
plot(density(temp[,18]), main="beta[2]",  xlim=c(-1.3,0.5))
plot(density(temp[,19]), main="beta[3]",  xlim=c(-1.5,0.3))
plot(density(temp[,20]), main="beta[4]",  xlim=c(-1.5,1))

#dev.off()
par(mfrow=c(2,2))
par(mar=c(4.1,4.1,2.1,2.1))
plot(density(temp[,27]), main="y[4,3]", xlim=c(1,5), xlab="log secs")
plot(density(temp[,28]), main="y[4,4]", xlim=c(1,5), xlab="log secs")
plot(density(exp(temp[,27])), main="exp y[4,3]", xlab="seconds", xlim=c(0,80))
plot(density(exp(temp[,28])), main="exp y[4,4]", xlab="seconds", xlim=c(0,80))

#  Those are somewhat automatic, not of great interest -- making sure
#  samples are ok, posteriors look sensible, and looking at the posterior 
#  predictions for two missing observations. 

#look at the different chains
temp2= Output2$Burnin.sims.array
dim(temp2)

par(mfrow=c(3,3))      # sets plots in a 3x3 grid
par(mar=c(3.1,4.1,2.1,2.1))   # removes space from between plots and increases
					# the size of the plots.  
acf(temp2[,1,1], main="", lag.max=35)     #autocorrelation plots
mtext("alpha[1]",side=3, line=1, cex=.8)
acf(temp2[,1,2], main="", lag.max=35) 
mtext("alpha[2]",side=3, line=1, cex=.8)
acf(temp2[,1,3], main="") 
mtext("alpha[3]",side=3, line=1, cex=.8)
acf(temp2[,1,4], main="") 
mtext("alpha[4]",side=3, line=1, cex=.8)
acf(temp2[,1,5], main="") 
mtext("alpha[5]",side=3, line=1, cex=.8)
acf(temp2[,1,6], main="") 
mtext("alpha[6]",side=3, line=1, cex=.8)
acf(temp2[,1,7], main="") 
mtext("alpha[7]",side=3, line=1, cex=.8)
acf(temp2[,1,8], main="") 
mtext("alpha[8]",side=3, line=1, cex=.8)

#dev.off()

par(mfrow=c(3,3))      # sets plots in a 3x3 grid
par(mar=c(3.1,4.1,2.1,2.1))   # removes space from between plots and increases
					# the size of the plots.  
acf(temp2[,1,1], main="", lag.max=200)     #autocorrelation plots
mtext("alpha[1]",side=3, line=1, cex=.8)
acf(temp2[,1,2], main="", lag.max=200) 
mtext("alpha[2]",side=3, line=1, cex=.8)
acf(temp2[,1,3], main="") 
mtext("alpha[3]",side=3, line=1, cex=.8)
acf(temp2[,1,4], main="") 
mtext("alpha[4]",side=3, line=1, cex=.8)
acf(temp2[,1,5], main="") 
mtext("alpha[5]",side=3, line=1, cex=.8)
acf(temp2[,1,6], main="") 
mtext("alpha[6]",side=3, line=1, cex=.8)
acf(temp2[,1,7], main="") 
mtext("alpha[7]",side=3, line=1, cex=.8)
acf(temp2[,1,8], main="") 
mtext("alpha[8]",side=3, line=1, cex=.8)
acf(temp2[1:1000*5,1,1], main="", lag.max=200)     #autocorrelation plots
mtext("alpha[1] lag 5",side=3, line=1, cex=.8)


#dev.off()

par(mfrow=c(3,3))      # re-sets plots in a 3x3 grid
#acf(temp2[,1,9], main="")     #you cannot plot this one, why?
# mtext("alphastep[1]",side=3, line=1, cex=.8)
acf(temp2[,1,10], main="") 
mtext("alphastep[2]",side=3, line=1, cex=.8)
acf(temp2[,1,11], main="") 
mtext("alphastep[3]",side=3, line=1, cex=.8)
acf(temp2[,1,12], main="") 
mtext("alphastep[4]",side=3, line=1, cex=.8)
acf(temp2[,1,13], main="") 
mtext("alphastep[5]",side=3, line=1, cex=.8)
acf(temp2[,1,14], main="") 
mtext("alphastep[6]",side=3, line=1, cex=.8)
acf(temp2[,1,15], main="") 
mtext("alphastep[7]",side=3, line=1, cex=.8)
acf(temp2[,1,16], main="") 
mtext("alphastep[8]",side=3, line=1, cex=.8)

#dev.off()

par(mfrow=c(3,3))      # re-sets plots in a 3x3 grid
acf(temp2[,1,26], main="")     #autocorrelation plots
mtext("tau.e",side=3, line=1, cex=.8)
acf(temp2[,1,25], main="") 
mtext("tau.b",side=3, line=1, cex=.8)
acf(temp2[,1,24], main="") 
mtext("sqrt(D)",side=3, line=1, cex=.8)
acf(temp2[,1,22], main="") 
mtext("rho",side=3, line=1, cex=.8)
acf(temp2[,1,27], main="") 
mtext("beta[1]",side=3, line=1, cex=.8)
acf(temp2[,1,18], main="") 
mtext("beta[2]",side=3, line=1, cex=.8)
acf(temp2[,1,27], main="") 
mtext("y[4,3]",side=3, line=1, cex=.8)
acf(temp2[,1,28], main="") 
mtext("y[4,4]",side=3, line=1, cex=.8)


###
### scatterplot matrices.  
### 

splom(temp[1:1000,1:8],pch=".")
splom(temp[1:1000,c(1,2,26,25,22,17)], pch=".")





###
### time series plots.  
###

#dev.off()

par(mfrow=c(4,1))
par(mar=c(4.1,4.1,2.1,2.1))

plot(1:1000,temp2[1:1000,1,1], type="l", xlab="iteration", ylab="")  
mtext("alpha[1]", side=3, line=1, cex=.8)
plot(1:1000,temp2[1:1000,1,2], type="l", xlab="iteration", ylab="")  
mtext("alpha[2]", side=3, line=1, cex=.8)
plot(1:1:1000,temp2[1:1000,1,3], type="l", xlab="iteration", ylab="")  
mtext("alpha[3]", side=3, line=1, cex=.8)
plot(1:length(temp2[,1,4]),temp2[,1,4], type="l", xlab="iteration", ylab="")  
mtext("alpha[4]", side=3, line=1, cex=.8)

#dev.off()

par(mfrow=c(4,1))
par(mar=c(4.1,4.1,2.1,2.1))

plot(1:1000,temp2[1:1000,1,5], type="l", xlab="iteration", ylab="")  
mtext("alpha[5]", side=3, line=1, cex=.8)
plot(1:1000,temp2[1:1000,1,6], type="l", xlab="iteration", ylab="")  
mtext("alpha[6]", side=3, line=1, cex=.8)
plot(1:1000,temp2[1:1000,1,9], type="l", xlab="iteration", ylab="")  
mtext("alphastep[1]", side=3, line=1, cex=.8)
plot(1:1000,temp2[1:1000,1,10], type="l", xlab="iteration", ylab="")  
mtext("alphastep[2]", side=3, line=1, cex=.8)


#dev.off()

par(mfrow=c(4,1))
par(mar=c(4.1,4.1,2.1,2.1))
plot(1:1000,temp2[1:1000,1,26], type="l", xlab="iteration", ylab="")  
mtext("tau.e", side=3, line=1, cex=.8)
plot(1:1000,temp2[1:1000,1,25], type="l", xlab="iteration", ylab="")  
mtext("tau.b", side=3, line=1, cex=.8)
plot(1:1000,temp2[1:1000,1,22], type="l", xlab="iteration", ylab="")  
mtext("rho", side=3, line=1, cex=.8)
plot(1:1000,temp2[1:1000,1,17], type="l", xlab="iteration", ylab="")  
mtext("beta[1]", side=3, line=1, cex=.8)




### Predictions

#dev.off()

### These next sets of plots are quite useful. 
### For EVERY one of these next plots, you need
### to be able to verbalize exactly what the 
### plot shows, and what the conclusion is from 
### the plot. 

par(mfrow=c(1,1))
plot(density(exp(temp[,28])), main="predictions", 
		xlab="seconds", ylab="Density", col="green", xlim=c(0,120))
lines(density(exp(temp[,27])), lty=2, col="blue")
legend(60, .05, c("y[4,3]", "y[4,4]"), lty=c(2,1), col=c("blue", "green"))

### Baseline medians

#dev.off()

plot(density(exp(temp[,1] + temp[,2]),bw=1), main="Baseline", 
		xlab="seconds", ylab="Density", col="green",lty=2,ylim=c(0,.15))
lines(density(exp(temp[,1]), bw=1), lty=1, col="blue")
legend(40, .15, c("Attender", "Distracter"), lty=c(1,2), 
	col=c( "blue","green"))

### Treatment effects for Attenders and Distracters

#Treatment effects among attenders 
plot(density(temp[,3]), main="Attenders")
lines(density(temp[,4]), lty=2, col="blue" )
lines(density(temp[,5]), lty=3, col="green")
abline(v=0)
legend(.3, 2.6, legend= c("Attend", "Distract", "Null"), 
		col=c("black", "blue", "green"), 
		lty=1:3)


#dev.off()
par(mfrow=c(1,2))

plot(density(temp[,3], bw=.1), main="", 
		xlab="Attenders", ylab="Density", col="blue", 
		lty=1, xlim=c(-1.3,1) )
lines(density(temp[,4], bw=.1), lty=2, col="green")
lines(density(temp[,5], bw=.1), lty=3, col="red")
legend(-1.3, 2.3, c("Attend", "Distract", "Null"), lty=c(1,2,3), 
	col=c( "blue","green", "red"), cex=.7)

plot(density(temp[,6], bw=.1), main="", 
		xlab="Distracters", ylab="Density", col="blue", 
		lty=1, xlim=c(-1.3,1) )
lines(density(temp[,7], bw=.1), lty=2, col="green")
lines(density(temp[,8], bw=.1), lty=3, col="red")
legend(-1.3, 2.3, c("Attend", "Distract", "Null"), lty=c(1,2,3), 
	col=c( "blue","green", "red"), cex=.7)

###
### exponential scale
### multiplicative effects.  
###

#dev.off()
par(mfrow=c(1,2))
plot(density(exp(temp[,3]), bw=.1), main="", 
		xlab="Attenders", ylab="Density", col="blue", 
		lty=1, xlim=c(0,2.5),ylim=c(0,2.7) )
lines(density(exp(temp[,4]), bw=.1), lty=2, col="green")
lines(density(exp(temp[,5]), bw=.1), lty=3, col="red")
legend(1.5, 2.5, c("Attend", "Distract", "Null"), lty=c(1,2,3), 
	col=c( "blue","green", "red"), cex=.7)

plot(density(exp(temp[,6]), bw=.1), main="", 
		xlab="Distracters", ylab="Density", col="blue", 
		lty=1, xlim=c(0,2.5),ylim=c(0,2.7) )
lines(density(exp(temp[,7]), bw=.1), lty=2, col="green")
lines(density(exp(temp[,8]), bw=.1), lty=3, col="red")
legend(1.5, 2.5, c("Attend", "Distract", "Null"), lty=c(1,2,3), 
	col=c( "blue","green", "red"), cex=.7)

###
### Treatment differences for attenders and distracters
### Calculating new functions of the parameters
###

#dev.off()

par(mfrow=c(1,2))

plot(density(temp[,4] - temp[,3], bw=.1), main="", 
		xlab="Attenders", ylab="Density", col="blue", 
		lty=1, xlim=c(-1.,1.5) )
lines(density(temp[,4] - temp[,5], bw=.1), lty=2, col="green")
lines(density(temp[,3] - temp[,5], bw=.1), lty=3, col="red")
legend(.7, 1.8, c("D - A", "D - N", "A - N"), lty=c(1,2,3), 
	col=c( "blue", "green", "red"), cex=.7)
abline(v=0, col="lightgrey")

plot(density(temp[,7] - temp[,6], bw=.1), main="", 
		xlab="Distracters", ylab="Density", col="blue", 
		lty=1, xlim=c(-1,1.5) )
lines(density(temp[,7] - temp[,8], bw=.1), lty=2, col="green")
lines(density(temp[,6] - temp[,8], bw=.1), lty=3, col="red")
legend(-1, 1.8, c("D - A", "D - N", "A - N"), lty=c(1,2,3), 
	col=c( "blue", "green", "red"), cex=.7)
abline(v=0, col="lightgrey")
title("Treatment Effect Differences", outer=T, line=-2)





