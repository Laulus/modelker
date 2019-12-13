## Not run: 
## remove (almost) everything in the working environment.
## You will get no warning, so don't do this unless you are really sure.
rm(list = ls())


## WORKING DIRECTORY
#work.dir<-paste('C:/Users/laulus/Documents/Mirror/THESE/MODELE-COMPLET/modelker')
#setwd(work.dir)

####------------------------------- PACKAGES -----------------------------------####
library(readr)
library(jagsUI)
library(rjags)
load.module('glm')
library(mcmcplots)

invlogit<-function(x) {1/(1+exp(-(x)))}
#invlogit <- function(x) exp(x) / (1+exp(x))

####------------------------------- DATA -----------------------------------####


## Load data
#data<-read.jagsdata("data/data.txt")
data<-read.bugsdata("data/data-all.txt")
# data<-read.bugsdata("data/data500ssX.txt") # data recent
#data<-read.bugsdata("data/data500-1.txt") # data all 
#data<-read.bugsdata("data/data500-2.txt") # data Eddy
#data<-read.bugsdata("data/data-all-norvegienne.txt") # 1297 poissons


### Recode rivers name
# 19 Château -- 1962 -> 1962
# 50 Norvégienne -- 1968 -> 1964
# 3 Albatros -- 1968 -> 1967
# 1 Acoena -- 1983 -> 1993
# 42 Manchots -- 1990 -> 1983
# 59 Rohan -- 2000 -> 2005
# 47 Nord -- 1986 -> 1989
#river.names <- c("CHA", "NORV", "ALB", "ACOE", "MAN","ROH", "NOR")
tmp<-as.factor(data$river)
levels(tmp)<-c("ACOENA","ALBATROS","CHATEAU","MANCHOT","NORD","NORVEGIENNE","ROHAN")
data$river <- tmp


# redefinition of the state of each individual
for (i in 1:length(data$age)){
  for (j in 1:max(data$age)){ #maxAge
    if (j <=data$age[i]){data$smolt[i,j]<-data$smolt[i,j]}
    else {data$smolt[i,j]<-NA} # we don't know the future state of the individual once it has been captured
  }}

# Compute vector with occasions of first capture
get.last <- function(x) max(which(!is.na(x)))
last <- apply(data$radius, 1, get.last) #

# # Max age in freshwater
# agefr=NULL
# for (i in 1:dim(data$radius)[1]) {
#   agefr[i] <- ifelse(last[i]>data$smoltAge[i],data$smoltAge[i], last[i] )
# }

# add values for age 0 (emergence) =  size and radius at emergence
bodySize <- cbind(rep(28,data$nScale),data$bodySize) ### size at emergence = 28 mm (see Eddy Beal Thesis)
radius <- cbind(rep(0,data$nScale),data$radius)
smolt <- cbind(rep(0,data$nScale),data$smolt)

### CLEANING
toremove <- which(
  data$age < data$age.tot # on retire les individus capturés entre aout et novembre / car dernier annuli non marqué / sosu-estimation âge
  | data$post < 0 # on retire les individus dont l'origine ou la date de colonisation est incertaine (i.e. age indv. > age colonization) # remove indi with age post-colonization <0
)

# Change indices (from 1 to ...):
a = data$ID[-toremove] # for individuals
ID = match(a, unique(sort(a)))
# ID<-factor(ID)
# levels(ID)<-1:data$nID
# ID <- as.numeric(ID)

tmp = data$river[-toremove] # for rivers
river = match(tmp, unique(sort(tmp)))



### Storing data for jags
dataToJags <- list(nScale=data$nScale - length(toremove)
                   ,ID = ID
                   ,nID = max(ID)
                   ,AGE=data$age[-toremove] +1
                   #,smoltAge=data$smoltAge[-toremove] +1
                   ,bodySize=bodySize[-toremove,]
                   ,radius=radius[-toremove,]
                   ,fishMaxSize=data$fishMaxSize[-toremove]
                   ,scaleMaxSize=data$scaleMaxSize[-toremove]
                   ,smolt=smolt[-toremove,]
                   ,invSmolt = 1-smolt[-toremove,]
                   ,m.bodysize = colMeans(bodySize[-toremove,],na.rm=TRUE)
                   ,river = river
                   ,nPOP = max(river)
                   #,colonizationAge = data$post[-toremove]
)
#save(dataToJags, file="data/dataToJags-all-norvegienne.Rdata")



####-----------------------------MODEL ----------------------------------####

write("
      model {
      
      for (i in 1:nScale) { # 20 scales belonging to 10 individuals, unbalanced. 
      
      for (age in 1:AGE[i]) { # maxAge of each scale
      
      #Y[i,age]<- ifelse(age>smoltAge[i],2,1) # 2: at sea / 1: in freshwater
      #Y[i,age]<- ifelse(age+1>smoltAge[i],2,1) # 2: at sea / 1: in freshwater
      #Y[i,age] <- 1
      
      
      ### LIKELIHOOD
      
      ## 1. Back-Calculation model
      bodySize[i,age] ~ dnorm(allo[i,age],tau[1]) #T(,fishMaxSize[i]+1)
      
      allo[i,age] <- a[1] + (fishMaxSize[i]-a[1])*(radius[i,age]/scaleMaxSize[i]) # Fraser-Lee model
      #allo[i,age] <- a[Y[i,age]] + ((fishMaxSize[i]-a[Y[i,age]])*pow((radius[i,age]/scaleMaxSize[i]),b[Y[i,age]])) #eddy / to check
      #allo[i,age] <-  fishMaxSize[i]+((a[Y[i,age]]+b[Y[i,age]]*radius[i,age])/(a[Y[i,age]]+b[Y[i,age]]*scaleMaxSize[i])) # body proportional hypothesis (BPH)
      #allo[i,age] <- (- a[Y[i,age]]/b[Y[i,age]]) + (fishMaxSize[i]+(a[Y[i,age]]/b[Y[i,age]]))*(radius[i,age]/scaleMaxSize[i]) # scale proportional hypothesis (SPH)
      
      } # END OF LOOP age
      
      ## 2. Reaction norm
      # Emergence (age = 1)
      #smolt[i,1]~dbern(p[i,1])	# edit jlabonne
      #p[i,1] <- m[i,1]
      #logit(m[i,1]) <- alpha[ID[i],river[i]]*1 + beta[ID[i],river[i]]*(bodySize[i,1] - m.bodysize[1]) +  gamma[ID[i],river[i]]*1*(bodySize[i,1] - m.bodysize[1])

      for (age in 2:AGE[i]) { # maxAge of each scale

      smolt[i,age]~dbern(p[i,age])	# edit jlabonne
      ##invSmolt[i,age-1]~ dbern(q)  # for now, we assume it is stochastic # edit jlabonne
      ## if state is already smolt, then p= 1. # edit jlabonne
      p[i,age] <- invSmolt[i,age-1]*(m[i,age]) + (1-invSmolt[i,age-1]) # edit jlabonne
      # m is what you really want to look at. # edit jlabonne

      logit(m[i,age]) <- alpha[ID[i],river[i]]*age + beta[ID[i],river[i]]*(bodySize[i,age] - m.bodysize[age]) +  gamma[ID[i],river[i]]*age*(bodySize[i,age] - m.bodysize[age])
      #logit(m[i,age]) <- alpha[ID[i],age,river[i]] + beta[ID[i],age,river[i]]*(bodySize[i,age] - m.bodysize[age]) # hierarchical model
      } # END OF LOOP age
      } # END OF LOOP i
      
      
      ### PRIORS
      a[1]~dunif(0, 100) # length of the fish at the time of scale formation
      # a[2] <- a[1] + delta[1]
      # a[2]<-0 # at that step isometry between scale and size of fish should be reached
      # delta[1]~dnorm(0, 0.01) # difference between a[1] and a[2]
      # delta[1]<-0
      # b[1]~dunif(0.1, 10) 
      # b[2] <- b[1] + delta[2]
      # delta[2]~dnorm(0, 0.01)  # difference between b[1] and b[2]
      
      tau[1]<-pow(sigma[1],-2)
      sigma[1]~dunif(0,100)

      tau[2]<-pow(sigma[2],-2)
      sigma[2]~dunif(0,100)
    
      
      # these are the individual random effects.
      # the mean of these things could be space and time dependent, assuming individuals 
      # deviate equally from the mean through space and time. 
      # these are poorly estimated with my fake data.
      # note that we don't really care about these effects anyway. 
      # means through space and time: that is what we care all about !
      # could have a A[m] where m is a population ID
      # A[m] <- alpha[m] +beta * colonizationAge[m] + gamma.
      # then you only have to check if alpha!=0 and/or beta!=0. 
      
    for (pop in 1:nPOP){  
      for (id in 1:nID){
      #alpha[id]~dnorm(A,tau[3])
      #beta[id]~dnorm(B,tau[4])
      #gamma[id]~dnorm(C,tau[5])

      alpha[id,pop]~dnorm(A[pop],tau[3])
      beta[id,pop]~dnorm(B[pop],tau[4])  #T(0,)
      gamma[id,pop]~dnorm(C[pop],tau[5])
      } # end of loop id
      
      A[pop]~dnorm(muA,tau[6])
      B[pop]~dnorm(muB,tau[7]) #T(0,)
      C[pop]~dnorm(muC,tau[8])
    } # end loop pop
      
      #A~dnorm(0,0.01)
      #B~dnorm(0,0.01) T(0,)
      #C~dnorm(0,0.01)
      
      # with Cauchy distribution: dt(0, pow(2.5,-2), 1)
      muA~dt(0, pow(2.5,-2), 1) #dnorm(0,0.01)
      muB~dt(0, pow(2.5,-2), 1) #dnorm(0,0.01) #T(0,)
      muC~dt(0, pow(2.5,-2), 1) #dnorm(0,0.01)
      
      for (k in 3:8){
      tau[k]~dgamma(0.1, 0.1)
      sigma[k] <- sqrt(1/tau[k])
      #tau[k] <- pow(sigma[k], -2)
      #sigma[k]~dunif(0,10)
      } # end of loop k
      

      
      # Assess model fit using a sum-of-squares-type discrepancy
      #for (i in 1:nScale){
      #sq[i] <- pow(bodySize[i,age[i]]-allo[i,age[i]], 2)      # Squared residuals
      #res[i] <- (bodySize[i,age[i]] - allo[i,age[i]])/sigma[1] # standardised residual
      #p.res[i] <- phi(res[i]) # p-value
      
      ## Generate replicate data and compute fit statistics for them
      #y.pred[i]~dnorm(allo[i,age[i]], tau[1])        # One new data set at each MCMC iteration
      #sq.pred[i] <- pow(y.pred[i] - allo[i,age[i]], 2)  # Squared residuals for new data
      #p.pred[i] <- step(bodySize[i] - y.pred[i]) 
      #} # end of loop nscale
      
      #fit <- sum(sq[])              # Sum of squared residuals for actual data set
      #fit.new <- sum(sq.pred[])      # Sum of squared residuals for new data set
      # test[j] <- step(fit.new[j]-fit[j])     # Test whether new data set more extreme
      # bpvalue[j] <- mean(test[j])         # Bayesian p-value
      
      } # END OF MODEL
      
      ", "code/mymodel5-jags.R")


## -------------------------- ANALYSIS ------------------------------#
nchains = 2 # Number of chains to run.
adaptSteps = 1000 # Number of steps to "tune" the samplers.
burnInSteps = 10000 # Number of steps to "burn-in" the samplers.
niter=20000 # Total number of steps in chains to save.
nthin=20

## PARAMETERS TO SAVE
parameters=c("a"
             #,"b"
             ,"sigma"
             #,"delta"
             , "A", "B", "C"
             , "muA", "muB", "muC"
)


#### ------------------------------INITS---------------------------------#####
inits<-function(){
  list(
    #sigma=c(10,10,.1,.1,.1)
    sigma=c(10,10,rep(NA,6))
    ,tau=c(NA,NA,rep(100,6))
    , a=30#c(30,NA)
    #, b=c(1.01,NA)
    # ,delta=0#c(0,0)
    #, A = -4
    #, B = -0.15
    #, C = 0.1
  )}


####---------------------ANALYSIS-----------------------####
## Compile & adapt
### Start of the run ###
## Compile & adapt
# Create, initialize, and adapt the model:

#jagsfit=jags.model(
#  'code/mymodel5-jags.R',
#  dataToJags,inits,
#  n.chains = nchains,
#  n.adapt = adaptSteps)

# # Burn-in:
# cat( "Burning in the MCMC chain.\n" )
# #jagsfit$recompile()
# update(jagsfit, n.iter=burnInSteps,progress.bar="text")
# 
# # The saved MCMC chain:
# cat( "Sampling final MCMC chain.\n" )
# fit.mcmc<-coda.samples(jagsfit,variable.names=parameters, n.iter=niter, thin=nthin)
# 
# ## BACKUP
# save(fit.mcmc,file=paste0("results/Results_","model5-jags","data-all",".RData"))
# 
# # To save individuals len:
# # len.mcmc <- coda.samples(jagsfit,
# # variable.names = c("mu.retro"),
# # n.iter = 1,thin = 1)
# pred.mcmc<-coda.samples(jagsfit,variable.names=c("y.pred"), n.iter=5000, thin=1)
# size.mcmc<-coda.samples(jagsfit,variable.names=c("bodySize[1:50,1:12]", "allo[1:50,1:12]"), n.iter=5000, thin=1)



### PARALLEL
#Call jags function; specify number of chains, number of adaptive iterations,
#the length of the burn-in period, total iterations, and the thin rate.
fit <- jags(data = dataToJags,
            inits = inits,
            parameters.to.save = parameters,
            model.file =    'code/mymodel5-jags.R',
            parallel=TRUE,
            n.chains = nchains,
            n.adapt = adaptSteps,
            n.iter = niter,
            n.burnin = burnInSteps,
            n.thin = nthin)
#Arguments will be passed to JAGS; you will see progress bars
#and other information
#Examine output summary
#fit
fit.mcmc <- fit$samples
## BACKUP
save(fit.mcmc,file=paste0("results/Results_","model5-jags","data-all",".RData"))


#####-------------------RESULTS----------------####

#summary(fit.mcmc)
# # summary(len.mcmc)
gelman.diag(fit.mcmc,multivariate=F) # test stat de m?lange chaine mcmc, convergence? multivariate psrf <1.1

#traplot(fit.mcmc,parameters)
traplot(fit.mcmc,"a")
traplot(fit.mcmc,"sigma")
traplot(fit.mcmc,c("A", "B", "C"))
traplot(fit.mcmc,c("muA", "muB", "muC"))

traplot(fit.mcmc,"C[2]")

#denplot(fit.mcmc,parameters)
denplot(fit.mcmc,"sigma")
denplot(fit.mcmc,c("A", "B", "C"))

caterplot(fit.mcmc,"A",reorder = FALSE)
caterplot(fit.mcmc,"B",reorder = FALSE)
caterplot(fit.mcmc,"C",reorder = FALSE)
caterplot(fit.mcmc,c("A", "B", "C"),reorder = FALSE)


mcmc<-as.matrix(fit.mcmc)
plot(mcmc[,"A[1]"],mcmc[,"B[1]"])
plot(mcmc[,"A[2]"],mcmc[,"B[2]"])
plot(mcmc[,"A[3]"],mcmc[,"B[3]"])

plot(mcmc[,"C[1]"],mcmc[,"B[1]"])
plot(mcmc[,"C[2]"],mcmc[,"B[2]"])
plot(mcmc[,"C[3]"],mcmc[,"B[3]"])
plot(mcmc[,"C[4]"],mcmc[,"B[4]"])
plot(mcmc[,"C[5]"],mcmc[,"B[5]"])

#mean(mcmc[,"delta[1]"]>0) # proportion de valeurs de alpha > 0
#mean(mcmc[,"delta[2]"]>0) # proportion de valeurs de alpha > 0





####--------------------FIGURES-------------------####
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



#### REACTION NORMS

res <- summary(fit.mcmc)[[2]]
A <- res[paste0("A[",1:7,"]"),"50%"]
B <- res[paste0("B[",1:7,"]"),"50%"]
C <- res[paste0("C[",1:7,"]"),"50%"]

bodySize <- seq(20,900,10)
age <- seq(1,14,1)
m <- array(,dim=c(length(age),length(bodySize),7))
for ( pop in 1:7){
  for ( j in 1:length(age)){
    for ( i in 1:length(bodySize)){
      m[j,i,pop] <- A[pop]*j + B[pop]*(bodySize[i] - dataToJags$m.bodysize[j]) +  C[pop]*j*(bodySize[i] - dataToJags$m.bodysize[j])
      #m[j,i,pop] <- invlogit(m[j,i,pop])
    }}}
m <- invlogit(m)



# library(plotly)
# plotList=list()
# for ( river in 1:7){
#   df <- data.frame(Age=dataToJags$age[dataToJags$river==river], BodySize=dataToJags$bodySize[which(dataToJags$bodySize > 28 & dataToJags$river==river)], River=river)
# #p <- plot_ly(
# #  x = ~age, 
# #  y = ~bodySize, 
# #  z = m[,,river] ,
# #  type = "contour" 
# #)
# plotList[[river]] =  plotly_build(plot_ly(  x = ~age
#                                             ,y = ~bodySize
#                                             ,z = m[,,river]
#                                             ,type = "contour" 
#                                             ,colors=c("white","steelblue","yellow","tomato")
#                                             , name=river)
#                                   ) %>% add_lines()
#   #add_trace(data = df, x = ~Age, y = ~BodySize, z = ~0, inherit=TRUE,mode = "markers",  
#   #          marker = list(size = 3.5, color = "red", symbol = 10))
# 
# }

#sbp = subplot(plotList, nrows=2)
#print(sbp)


par(mfrow = c(3,3))
for ( river in 1:7){
  last <- apply(dataToJags$bodySize, 1, get.last) #
  last.age <-last[dataToJags$river==river]
  lf <- dataToJags$bodySize[dataToJags$river==river,]
  smolt <- dataToJags$smolt[dataToJags$river==river,]
  z <- m[,,river]
  #z <- t(m[,,river])
  df <- data.frame(Age=last.age, BodySize=lf, River=river, Smolt=smolt)
  
  cols = (colorRampPalette(c("white","steelblue","tomato"))(10))
  col.smolt=c("#00BFFF", "#3A5FCD")
  contour(
    x=age #nrow(z)
    ,y=bodySize #ncol(z)
    ,z=z
    ,ylab="Size (mm)",xlab="Age", main=levels(data$river)[river],col=cols)
  #points(dataToJags$AGE[dataToJags$river==river],dataToJags$bodySize[which(dataToJags$bodySize > 28 & dataToJags$river==river)], pch=16, col=col.smolt[df$Smolt+1]) #<- New
  
  for ( i in 1:dim(lf)[1]){
  for ( j in 1:dim(lf)[2]){
  points(j,lf[i,j], pch=16, col=col.smolt[smolt[i,j]+1]) #<- New
  }
  }
  }

