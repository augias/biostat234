#LOAD NECESSARY PACKAGES
library(R2jags)
library(lattice) # Needed for scatterplot matrix
# Set working directory
setwd("~/Documents/UCLA/Winter 2021/biostat234/lab 2/")
getwd()

# Save the data file to your working directory. 

#READ IN DATA
housing=read.table("housingdata2.txt")
housing   ## Always look at your data
## Anything funny about any of the columns?

#Give the columns useful names 
colnames(housing) <- c("cost", "eaves", "windows", "yard", "roof")

#SEPARATE X & Y
y<-housing[,1]
x<-as.matrix(housing[,2:5])

## Remember: look at your data.  

## Did you look at the recently defined x and y?  

#Create the model

reg = lm(y~x)
summary(reg) # classical regression
# There's a problem, what is it?
# Why is there a problem?

sink("housingmodel.txt")
cat("
model
{ 
   for(i in 1:N) {
	     y[i] ~ dnorm( mu[i] , tau )
	     mu[i] <- beta0 + inprod(x[i,] , beta[] )
		}
	      
	 beta0 ~ dnorm( mbeta0 , precbeta0)
	 
for (j in 1:K) {
	 beta[j] ~ dnorm( m[j] , prec[j] )
		}
	   tau ~ dgamma( tau.a , tau.b )
	   sigma <- 1 / sqrt( tau )
	}
  ",fill = TRUE)
sink()

# This copies the jags program to a file. Please take 
# the time to understand the model.

#DIFFERENT PRIORS TO TRY
dataA<-list(N=21, K=4, m=c(1.6053, 1.2556, 2.3413, 3.6771), 
prec=c(.2164, .1105, .2061, .1337), tau.a=17,
tau.b=1128, mbeta0= -5.682, precbeta0=.05464, x=x, y=y)

dataB<-list(N=21, K=4, m=c(1.6053, 1.2556, 2.3413, 3.6771), 
prec=c(.02774, .014160, .02642, .01714), tau.a=2.1795,
tau.b=144.6, mbeta0= -5.682, precbeta0=.007005, x=x, y=y)

dataC<-list(N=21, K=4, m=c(1.6053, 1.2556, 2.3413, 3.6771), 
prec=c(.005549, .002832, .005284, .003428), tau.a=.4359,
tau.b=28.92, mbeta0= -5.682, precbeta0=.00140, x=x, y=y)

#SET UP INITAL VALUES
# This creates a list with 5 copies of the initial values.
# Change the number 5 to match the n.chains variable in 
# the jags call.  

inits<-rep(list(list(beta0=0, beta=c(1,1,1,1),tau=1)),5)

#DEFINE PARAMETERS TO MONITOR
parameters <- c("beta0", "beta" , "tau")

#RUN THE JAGS PROGRAM, SAVING DATA TO LAB2.SIM
set.seed(1234)
lab2.sim <- jags (dataA, inits, parameters, "housingmodel.txt", n.chains=5, 
	n.iter=5100, n.burnin=100, n.thin=1, DIC=FALSE)

#
print(lab2.sim)
plot(lab2.sim)
temp2=lab2.sim$BUGSoutput$sims.matrix
head(temp2)
interceptbeta = temp2[,5]
eavesbeta = temp2[,1]
windowsbeta = temp2[,2]
yardbeta = temp2[,3]
roofbeta = temp2[,4]
tau = temp2[,6]
plot(acf(interceptbeta))
plot(acf(tau))

length(interceptbeta)
plot(interceptbeta, type="l")	#Trace Plot		
plot(acf(eavesbeta))		#Autocorrelation plot
plot(density(yardbeta))	#	Density plot
plot(interceptbeta,yardbeta)		#Single correlation plot of beta0 v 1
splom(temp2[1:5000,1:5],pch=".")	#	Scatterplot matrix of correlation plots
          # 25000 takes a long time. "." is a better character to plot with
mean(eavesbeta)
sd(eavesbeta)
length(eavesbeta[eavesbeta>0])/length(eavesbeta)
apply(temp2, 2, mean)			#apply command is very handy and worth knowing.
apply(temp2, 2, sd)
apply(temp2, 2, quantile, c(.025,.975))

# Extra.  
# We can define our own summary measures

mysummary = function(invector) {
c(mean(invector), sd(invector), quantile(invector, .025), quantile(invector,.975),
	length(invector[invector>0])/length(invector))
}
apply(temp2, 2, mysummary)

# I use this function often.
t(round(apply(temp2, 2, mysummary),2))
# t() is the transpose function
# rounding to two digits accuracy is fine for everything except tau
# and the probability >0 might wish to have an extra digit. 

# Why (or why not?) does it make sense to plot all four regression coefficients 
# on a single plot in this particular example?

plot(density(yardbeta), main="Housing Data Posteriors", 
      xlab="Regression Coefficient", ylab="Density")  #	Density plot
lines(density(eavesbeta), col="blue")
lines(density(windowsbeta), col="red")
lines(density(roofbeta), col="green")

# What is the interpretation of this plot?

###############################################
# Extra stuff
#

######
# 
# Generating the prior predictor distribution for the example of lab 2.
# 
# 

#1.  Move y from the data to the inits.   

dataPP<-list(N=21, K=4, m=c(1.6053, 1.2556, 2.3413, 3.6771), 
             prec=c(.2164, .1105, .2061, .1337), tau.a=17,
             tau.b=1128, mbeta0= -5.682, precbeta0=.05464, x=x)  # remove y

initsPP<-rep(list(list(beta0=0, beta=c(1,1,1,1),tau=1, y=y)),5)

#2.  Add y to the parameters that we track
parametersPP <- c("beta0", "beta" , "tau", "sigma", "y")

#3.  Include the keyword DIC=FALSE


lab2.simPP <- jags (dataPP, initsPP, parametersPP, "housingmodel.txt", n.chains=5, 
                    n.iter=5100, n.burnin=100, n.thin=1, DIC=FALSE)


print(lab2.simPP)

  
  # The means and sd for the parameters: regression coefficients, 
  # sigma/tau should exactly reproduce the prior.  
  #  
  # The distributions for the data are predictions of what you would expect
  # to see prior to seeing the data.  
  # The predictions do use the x values but do not use the y values (except to initialize
  # the computations.)  

#####################
#####################

# In a related note, when we create the futurefit, futureobs and futuretail
# variables in the regular regression, we need to (i) add an initial value for 
# futureobs and (ii) include the three variables in our output parameters.  

# Added jags code:  
futurefit <- beta0 + beta[1] + beta[2] + beta[3]*2 + beta[4]*2
futureobs ~ dnorm(futurefit, tau)
futuretail <- beta0 + beta[1] + beta[2] + beta[3]*2 + beta[4]*2 + 1.645*sigma

# New example R commands:
inits<-rep(list(list(beta0=0, beta=c(1,1,1,1),tau=1, futureobs=10)),5)
parameters <- c("beta0", "beta" , "tau", "sigma", "futurefit", "futureobs", "futuretail")



