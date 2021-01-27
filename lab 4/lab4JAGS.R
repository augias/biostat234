#LOAD NECESSARY PACKAGES
library(R2jags)
load("AddBurnin.RData")
library(lattice)

#CHANGE WORKING DIRECTORY
setwd("C:\\Users\\Rob\\Courses\\Bayes\\labsjags\\lab4")
getwd()

# useful function
mysummary = function(invector) {
c(mean(invector), sd(invector), quantile(invector, .025), 
	quantile(invector,.975),
	length(invector[invector>0])/length(invector))
}

# load the data.  

load("lab4_data.RData")

# inspect the data please.  
# variables are
#  y
#  x
#  priordata   # note:  long!
#  inits       # even longer.  look at just the first one!  

#Create the model

sink("lab4model.txt")
cat("
model
        {               
                for( i in 1 : 64 ) {
                        for( j in 1 : 4 ) {
                            s[i, j]<-4*(i-1)+j
                            y[i, j] ~ dnorm(mu[i , j],tau.e)
                                mu[i , j] <- inprod(x[s[i,j],],alpha[])+beta[i]
                        }
                        beta[i]~dnorm(0, tau.b)
                }

for( k in 1:8) {
                alpha[k]~dnorm(m[k],varinv[k])
                alphasign[k] <- step(alpha[k])
}

                tau.e ~ dgamma(ea,eb)
                tau.b~dgamma(ba,bb)

                sigma <- 1 /sqrt( tau.e)
                sqrtD <- 1 /sqrt( tau.b)
                rho <- sqrtD*sqrtD/(sigma*sigma + sqrtD *sqrtD)

        }

    ",fill = TRUE)
sink()

proc.time()
run1 = jags(priordata, inits, parameters, "lab4model.txt", 
	n.chains=5, n.iter=1100, n.burnin=0, n.thin=1)
proc.time()
# 1100 iterations takes about 3 seconds on my computer.  
# 11000 iterations takes a little under 15 seconds on my computer.  

names(run1)
Output1=AddBurnin(run1$BUGSoutput$sims.array,burnin=100,n.thin=1)

print(Output1$Burnin.Summary)





