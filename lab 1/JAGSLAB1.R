
# Next command, uncomment if R2jags not already installed. 
# install.packages("R2jags")
library(R2jags)

setwd("C:\\Users\\Rob\\Courses\\Bayes\\labsjags\\lab1")
getwd()

## set your working directory with setwd()
## Note double backslashes. 
## Your working directory here. 
## Suggested: pen/flash drive or desktop while in computer lab.  
## To identify a directory where you have a file you need, 
## drag and drop the file into the R Console -- it 
## (hopefully) won't run but R will give you the proper directory and format
## of the file name, which you can copy and paste.

## Create jags model file. This uses the bugs language 
## to specify the prior times likelihood.
## Our model is simple linear regression
## This is virtually the same model in the jags manual section 2.1. 
## The prior for beta is different here. 
## Our prior for tau is "superior". 

sink("model1.txt")
cat("
    model
	{
	for(i in 1:N){
		y[i] ~ dnorm(mu[i],tau)
		mu[i] <- alpha + beta * (x[i]-x.bar)
		}
		alpha ~ dnorm(0, 0.0001)
		beta ~ dnorm(1,1)
		tau ~ dgamma(.25,.25)
    sigma <- 1/sqrt(tau)
	}
    ",fill = TRUE)
sink()

## input data
## We specify a simple linear regression data set with 5 observations 
## Here is the data:

x     = c(1,2,3,4,5)
y     = c(1,3,3,3,5)
N     = 5 
x.bar = 3
jags.data = list("x","y","N","x.bar")

##
## Input parameter names
## These are the parameters whose posterior distribution(s) we are interested in.
## If we aren't interested in a particular parameter and don't list it, JAGS will 
## not keep track of information about that parameter, though it is still part of 
## the posterior. Here we list all parameters, including both tau and sigma, though
## they are functions of each other. 

jags.params = c("alpha", "beta", "tau", "sigma")

## TODO: identify where each parameter is used in the model. 

##
## Input initital values
## These are starting values for the fitting algorithm. 
## They are not part of the model or data set. 
## They are NOT part of our assumptions. 

jags.inits = function(){
    list("alpha" = 0, "beta" = 1, "tau" = 1)
}

##
## Fit model
## This assumes that the model1.txt file is in the working directory. 
## The variable "lab1.sim" will contain all the output from the jags run. 

lab1.sim = jags(jags.data, jags.inits, jags.params, 
              model.file = "model1.txt", 
              n.chains = 3, n.iter = 11000, n.burnin = 1000)

##
## Tools are provided to automatically summarize the results. 
## These tools aren't perfect, and we will need to edit/adapt the output
## depending on our needs. 

print(lab1.sim)

plot(lab1.sim)

## For easier access, we grab the posterior samples for each parameter. 
temp1 = lab1.sim$BUGSoutput$sims.matrix
head(temp1[,1])

alpha = temp1[,1]
beta  = temp1[,2]
tau   = temp1[,3]
sigma = temp1[,4]

## Example summaries for various parameters. 
## quantiles, density plots, and two time series plots. 

quantile(alpha, .1)

plot(density(beta))

plot(tau)
plot(tau , type = "l")

plot(sigma[1:1000] , type = "l")

## Other useful code to try out.
## What does each command do?

length(alpha)   ## How many samples? 

mean(alpha)
sd(alpha)
quantile(alpha,c(.025,.25, .5, .75, .975))

plot(density(alpha))
plot(sigma, type = "l")
autocorr.plot(sigma)

## repeat these last 3 for each parameter in turn. 

## TODO: Try out the alternative jags.inits definition
jags.inits = function(){
  list("alpha" = rnorm(1), "beta" = rnorm(1),  "tau" = runif(1))
}

## Rerun, saving the jags output in a different file. 
## Can you find _any_ important differences in your
## conclusions? 
## (Note 1: How might you evaluate if there are any differences
## important or not between the two runs? )
## (Note 2: Starting values are part of the algorithm, not
## the model or data, so they aren't supposed to affect the conclusions.
## What would it mean if the starting values had a major effect
## on conclusions? )
