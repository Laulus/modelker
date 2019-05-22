## Not run: 
## remove (almost) everything in the working environment.
## You will get no warning, so don't do this unless you are really sure.
rm(list = ls())
## WORKING DIRECTORY

# work.dir<-paste('D:/THESE/ANALYSES/modele-complet/modelker')

# work.dir<-paste('C:/Users/laulus/Documents/Mirror/THESE/MODELE-COMPLET/modelker')

#setwd(work.dir)
library(rjags)

#--------------- DATA --------------#
library(readr)

#data<-read.jagsdata("data/data.txt")
#data<-read.bugsdata("data/data.txt")
# data<-read.bugsdata("data/data500ssX.txt") # data recent sample
# data<-read.bugsdata("data/data500-1.txt") # data all sample
#data<-read.bugsdata("data/data500-2.txt") # data Eddy sample
data<-read.bugsdata("data/data-all.txt") # data all (no sampling) # 622 recent fish + 2990 fish

# Compute vector with occasions of first capture
get.last <- function(x) max(which(!is.na(x)))
last <- apply(data$radius, 1, get.last) #

# # Max age in freshwater
# agefr=NULL
# for (i in 1:dim(data$radius)[1]) {
#   agefr[i] <- ifelse(last[i]>data$smoltAge[i],data$smoltAge[i], last[i] )
# }

# add values for age 0 (emergence) =  size and radius at emergence
bodySize <- cbind(rep(28,data$nScale),data$bodySize)
radius <- cbind(rep(0,data$nScale),data$radius)

# on retire les individus capturés entre aout et novembre / car dernier annuli non marqué / sosu-estimation âge
toremove <- which(data$age < data$age.tot)
dataToJags <- list(nScale=data$nScale - length(toremove)
                  ,age=data$age[-toremove] +1
                 ,smoltAge=data$smoltAge[-toremove] +1
                 ,bodySize=bodySize[-toremove,]
                 ,radius=radius[-toremove,]
                 ,fishMaxSize=data$fishMaxSize[-toremove]
                 ,scaleMaxSize=data$scaleMaxSize[-toremove]
                 ,smolt=data$smolt[-toremove,]
                 )

# TO RUN FOR eddy data ,because toremove is integer(0)
# dataToJags <- list(nScale=data$nScale
#                    ,age=data$age
#                    ,smoltAge=data$smoltAge
#                    ,bodySize=data$bodySize
#                    ,radius=data$radius
#                    ,fishMaxSize=data$fishMaxSize
#                    ,scaleMaxSize=data$scaleMaxSize
#                    ,smolt=data$smolt
# )

# dataToJags <- list(nScale=data$nScale
#                    ,age=last#agefr #data$age
#                    ,smoltAge=data$smoltAge
#                    ,bodySize=data$bodySize
#                    ,radius=data$radius
#                    ,fishMaxSize=data$fishMaxSize
#                    ,scaleMaxSize=data$scaleMaxSize
# )


plot(NULL,xlim=c(0,max(data$radius,na.rm=T)),ylim=c(0,max(data$bodySize, na.rm=T)),xlab="radius",ylab="size")
points(data$radius,data$bodySize)


##-----------------------------MODEL ----------------------------------##

write("
      model {
      
     for (i in 1:nScale) { # 20 scales belonging to 10 individuals, unbalanced. 
     
           Y[i]<- ifelse(age[i]>smoltAge[i],2,1) # 2: at sea / 1: in freshwater
           
		for (j in 1:age[i]) { # maxAge of each scale
      
      ### likelihood
      # 1. retrocalcul
      bodySize[i,j] ~ dnorm(allo[i,j],tau[1]) #T(,fishMaxSize[i]+1)
      allo[i,j] <- a[Y[i]] + (X[i,j]) 
       
      X[i,j] <- (fishMaxSize[i]-a[Y[i]])*(radius[i,j]/scaleMaxSize[i]) # Fraser-Lee model
      #X[i,j] <- ((fishMaxSize[i]-a[Y[i]])*(radius[i,j]/scaleMaxSize[i]))^b[Y[i]] #eddy / to check
      
      } # END OF LOOP j
      } # END OF LOOP i
      
      ### priors
      a[1]~dunif(0, 100) # length of the fish at the time of scale formation
      a[2] <- a[1] + delta[1]
      delta[1]~dnorm(0, 0.01) # difference between a[1] and a[2]

      b[1]~dunif(0.9, 1.1) 
      b[2] <- b[1] + delta[2]
      delta[2]~dnorm(0, 0.01)  # difference between b[1] and b[2]
      
      tau[1]<-1/(sigma[1]^2)
      sigma[1]~dunif(0,100)
      
      # Assess model fit using a sum-of-squares-type discrepancy
      for (i in 1:nScale){
      #sq[i] <- pow(bodySize[i,age[i]]-allo[i,age[i]], 2)      # Squared residuals
      #res[i] <- (bodySize[i,age[i]] - allo[i,age[i]])/sigma[1] # standardised residual
      # p.res[i] <- phi(res[i]) # p-value
      
      ## Generate replicate data and compute fit statistics for them
      y.pred[i]~dnorm(allo[i,age[i]], tau[1])        # One new data set at each MCMC iteration
      #sq.pred[i] <- pow(y.pred[i] - allo[i,age[i]], 2)  # Squared residuals for new data
      #p.pred[i] <- step(bodySize[i] - y.pred[i]) 
      }

      #fit <- sum(sq[])              # Sum of squared residuals for actual data set
      #fit.new <- sum(sq.pred[])      # Sum of squared residuals for new data set
      # test[j] <- step(fit.new[j]-fit[j])     # Test whether new data set more extreme
      # bpvalue[j] <- mean(test[j])         # Bayesian p-value
      
      
      } # END OF MODEL
      
      ", "code/mymodel3-jags.R")

#____________________________________________________________________#

## -------------------------- ANALYSIS ------------------------------#

require(rjags)

nChains = 2 # Number of chains to run.
adaptSteps = 1000 # Number of steps to "tune" the samplers.
burnInSteps = 5000 # Number of steps to "burn-in" the samplers.
niter=10000 # Total number of steps in chains to save.


## PARAMETERS TO SAVE
parameters=c("a","b","sigma","delta")

#______________________________________________________________________#

## ------------------------------INITS---------------------------------#

inits<-function(){
  list(sigma=c(10)
       , a=c(30,NA)
       , b=c(1.01,NA)
       ,delta=c(0,0)
       )}


#________________________________________________________________________#
## Compile & adapt
### Start of the run ###
## Compile & adapt
# Create, initialize, and adapt the model:

jagsfit=jags.model(
  'code/mymodel3-jags.R',
  dataToJags,inits,
  n.chains=nChains,
  n.adapt = adaptSteps)

# Burn-in:
cat( "Burning in the MCMC chain.\n" )
#jagsfit$recompile()
update(jagsfit, n.iter=burnInSteps,progress.bar="text")

# The saved MCMC chain:
cat( "Sampling final MCMC chain.\n" )
fit.mcmc<-coda.samples(jagsfit,variable.names=parameters, n.iter=niter, thin=2)
# To save individuals len:
# len.mcmc <- coda.samples(jagsfit,
# variable.names = c("mu.retro"),
# n.iter = 1,thin = 1)
pred.mcmc<-coda.samples(jagsfit,variable.names=c("y.pred"), n.iter=niter, thin=1)
size.mcmc<-coda.samples(jagsfit,variable.names=c("bodySize[1:50,1:12]", "allo[1:50,1:12]"), n.iter=niter, thin=1)


## BACKUP
save(fit.mcmc,file=paste0("results/Results-all_","model3-jags",".RData"))


#______________________________________________________________________________#
# EXAMINE THE RESULTS

#summary(fit.mcmc)
# # summary(len.mcmc)
#gelman.diag(fit.mcmc,multivariate=F) # test stat de m?lange chaine mcmc, convergence? multivariate psrf <1.1
# 
# 
library(mcmcplots)
# traplot(fit.mcmc,"a")
# traplot(fit.mcmc,"b")
# traplot(fit.mcmc,"sigmae")

traplot(fit.mcmc,parameters)
denplot(fit.mcmc,parameters)
caterplot(fit.mcmc,"a")
caterplot(fit.mcmc,"b")

mcmc<-as.matrix(fit.mcmc)
plot(mcmc[,"a[1]"],mcmc[,"a[2]"])
plot(mcmc[,"a[1]"],mcmc[,"sigma"])
plot(mcmc[,"b[1]"],mcmc[,"sigma"])

mean(mcmc[,"delta[1]"]>0) # proportion de valeurs de alpha > 0
mean(mcmc[,"delta[2]"]>0) # proportion de valeurs de alpha > 0


## fit
mcmc<-as.matrix(pred.mcmc)
pred <- apply(mcmc,2,quantile, probs=c(0.025, .5, 0.975))

plot(NULL,xlim=c(0,800),ylim=c(0,800), xlab="Observed",ylab="Predicted")
for (i in 1:dataToJags$nScale){
  points(dataToJags$bodySize[i,last[i]+1], pred["50%",i], col="lightgrey")
  segments(dataToJags$bodySize[i,last[i]+1], pred["2.5%",i], dataToJags$bodySize[i,last[i]+1], pred["97.5%",i], col="lightgrey")
}
abline(0, 1, lwd = 1, col = "black")


# Plot size
mcmc<-as.matrix(size.mcmc)
res <- summary(size.mcmc)[[2]]
#size <- apply(mcmc, 2,quantile, probs=0.5)

plot(NULL,xlim=c(0,dim(dataToJags$bodySize)[2]),ylim=c(0,1000), xlab="Age",ylab="Body size")

last <- apply(dataToJags$bodySize, 1, get.last) #
for (i in 1:dataToJags$nScale){
points(rnorm(1,0,.1)+last[i]-1, dataToJags$bodySize[i,last[i]], col="lightgrey")
}
points(0:(max(last)-1),apply(dataToJags$bodySize,2,mean, na.rm=TRUE),pch=16,col="darkgrey",cex=2)


for (i in 1:50){
  
  points(dataToJags$smoltAge[i]-1,dataToJags$bodySize[i,dataToJags$smoltAge[i]],pch=17,col=i)
  points(dataToJags$age[i]-1,dataToJags$bodySize[i,dataToJags$age[i]],pch=16,col=i)
  text(dataToJags$age[i]-1,dataToJags$bodySize[i,dataToJags$age[i]]+1,i)
  points(1:dataToJags$age[i]-1,res[paste0("bodySize[",i,",",1:dataToJags$age[i],"]"),"50%"],col=i, type="b")
  
  #segments(1:dataToJags$age[i],res[paste0("bodySize[",i,",",1:dataToJags$age[i],"]"),"2.5%"],1:dataToJags$age[i],res[paste0("bodySize[",i,",",1:dataToJags$age[i],"]"),"97.5%"],col=i)
  }
#points(0:(max(last)-1),summary(size.mcmc)[[2]][paste0("bodySize[",i,",",1:dataToJags$age[i],"]"),"50%"],col=i, type="b")


# ### ALLO
# mcmc<-as.matrix(size.mcmc)
# #mcmc<-as.matrix(allo.mcmc)
# smp <- sample(1:niter,1)
# col <- c("steelblue1","royalblue3")
# plot(NULL,xlim=c(0,dim(data$bodySize)[2]),ylim=c(0,1000), xlab="Age",ylab="Body size")
# for (i in 1:20){
# 
#   points(dataToJags$smoltAge[i],dataToJags$bodySize[i,dataToJags$smoltAge[i]],pch=17,col=i)
#   points(dataToJags$age[i],dataToJags$bodySize[i,dataToJags$age[i]],pch=16,col=i)
#   #points(1:dataToJags$age[i],mcmc[smp, paste0("bodySize[",i,",",1:dataToJags$age[i],"]")],col=i, type="b")
#   points(1:dataToJags$age[i],exp(mcmc[smp, paste0("allo[",i,",",1:dataToJags$age[i],"]")]),col=i, type="b")
# 
# }



# # RECONSTRUCTION DES mu_retro # solution temporaire, pas propre
# 
# a.retro<-summary(fit.mcmc)$quantiles[paste0("a.retro"),"50%"]
# b.retro<-summary(fit.mcmc)$quantiles[paste0("b.retro"),"50%"]
# 
# # mu.retro<-(summary(len.mcmc)$empirical[paste0("mu.retro"),"Mean"])
# 
# 
# mu.retro<-rep(0,length(dataframe$lf))
# for ( i in 1:length(dataframe$lf)){
#   mu.retro[i]<-a.retro + ((dataframe$lf[i]*dataframe$RM[i]/dataframe$RT[i])-a.retro)*(((dataframe$radi[i])/dataframe$RM[i])^b.retro)
# }
# dataframe<-cbind(dataframe,mu.retro)
# save(dataframe,file="dataframeeddy.RData")
