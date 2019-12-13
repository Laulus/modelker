model {
  
  for (i in 1:nScale) { # 20 scales belonging to 10 individuals, unbalanced. 
    
    #Y[i,j]<- ifelse(age[i]>smoltAge[i],2,1) # 2: at sea / 1: in freshwater
    
    for (j in 1:age[i]) { # maxAge of each scale
      Y[i,j]<- ifelse(j>smoltAge[i],2,1) # 2: at sea / 1: in freshwater
      # Y[i,j]<- ifelse(j+1>smoltAge[i],2,1) # 2: at sea / 1: in freshwater
      
      ### LIKELIHOOD
      ## 1. Back-Calculation models
      bodySize[i,j] ~ dnorm(allo[i,j],tau[Y[i,j]]) #T(,fishMaxSize[i]+1)
      
      allo[i,j] <- a[Y[i,j]] + (fishMaxSize[i]-a[Y[i,j]])*(radius[i,j]/scaleMaxSize[i]) # Fraser-Lee model
      #allo[i,j] <- a[Y[i,j]] + ((fishMaxSize[i]-a[Y[i,j]])*pow((radius[i,j]/scaleMaxSize[i]),b[Y[i,j]])) #eddy / to check
      #allo[i,j] <-  fishMaxSize[i]+((a[Y[i,j]]+b[Y[i,j]]*radius[i,j])/(a[Y[i,j]]+b[Y[i,j]]*scaleMaxSize[i])) # body proportional hypothesis (BPH)
      #allo[i,j] <- (- a[Y[i,j]]/b[Y[i,j]]) + (fishMaxSize[i]+(a[Y[i,j]]/b[Y[i,j]]))*(radius[i,j]/scaleMaxSize[i]) # scale proportional hypothesis (SPH)
      
    } # END OF LOOP j
    
    ## 2, Reaction norm
    
    # Age 1
    smolt[i,1]~dbern(p[i,1])	# edit jlabonne
    p[i,1] <- m[i,1]
    logit(m[i,1]) <- alpha[ID[i]]*1 + beta[ID[i]]*(bodySize[i,1] - m.bodysize[1]) +  gamma[ID[i]]*1*(bodySize[i,1] - m.bodysize[1])
    
    for (j in 2:age[i]) { # maxAge of each scale
      
      smolt[i,j]~dbern(p[i,j])	# edit jlabonne
      ##invSmolt[i,j-1]~ dbern(q)  # for now, we assume it is stochastic # edit jlabonne
      ## if state is already smolt, then p= 1. # edit jlabonne
      p[i,j] <- invSmolt[i,j-1]*(m[i,j]) + (1-invSmolt[i,j-1]) # edit jlabonne
      # m is what you really want to look at. # edit jlabonne
      
      logit(m[i,j]) <- alpha[ID[i]]*j + beta[ID[i]]*(bodySize[i,j] - m.bodysize[j]) +  gamma[ID[i]]*j*(bodySize[i,j] - m.bodysize[j])
      #logit(m[i,j]) <- alpha[ID[i],j] + beta[ID[i],j]*(bodySize[i,j] - m.bodysize[j]) # hierarchical model
    } # END OF LOOP j
  } # END OF LOOP i
  
  ### PRIORS
  a[1]~dunif(0, 100) # length of the fish at the time of scale formation
  # a[2] <- a[1] + delta[1]
  a[2]<-0 # at that step isometry between scale and size of fish should be reached
  # delta[1]~dnorm(0, 0.01) # difference between a[1] and a[2]
  # delta[1]<-0
  
  #b[1]~dunif(0.1, 10) 
  #b[2] <- b[1] + delta[2]
  #delta[2]~dnorm(0, 0.01)  # difference between b[1] and b[2]
  
  tau[1]<-1/(sigma[1]^2)
  sigma[1]~dunif(0,100)
  
  tau[2]<-1/(sigma[2]^2)
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
  
  for (id in 1:nID){
    alpha[id]~dnorm(A,tau[3])
    beta[id]~dnorm(B,tau[4])
    gamma[id]~dnorm(C,tau[5])
    
    #alpha[id,j]~dnorm(A[j],tau[3])
    #beta[id,j]~dnorm(B[j],tau[4])
    #gamma[id,j]~dnorm(C[j],tau[5])
  } # end of loop id
  
  #	taua~dgamma(0.001, 0.001) #a bit informative, else it runs like shit. 
  #	taub~dgamma(0.001, 0.001)
  #	tauc~dgamma(0.001, 0.001)
  
  for (k in 3:5){
    tau[k]<-1/(sigma[k]^2)
    sigma[k]~dunif(0,10)
  } # end of loop k
  
  
  A~dnorm(0,0.01) #T(-10,0)
  B~dnorm(0,0.01) #T(-1,1)
  C~dnorm(0,0.01)	T(-1,1)
  
  
  # Assess model fit using a sum-of-squares-type discrepancy
  for (i in 1:nScale){
    #sq[i] <- pow(bodySize[i,age[i]]-allo[i,age[i]], 2)      # Squared residuals
    #res[i] <- (bodySize[i,age[i]] - allo[i,age[i]])/sigma[1] # standardised residual
    # p.res[i] <- phi(res[i]) # p-value
    
    ## Generate replicate data and compute fit statistics for them
    y.pred[i]~dnorm(allo[i,age[i]], tau[Y[i,age[i]]])        # One new data set at each MCMC iteration
    #sq.pred[i] <- pow(y.pred[i] - allo[i,age[i]], 2)  # Squared residuals for new data
    #p.pred[i] <- step(bodySize[i] - y.pred[i]) 
  } # end of loop nscale
  
  #fit <- sum(sq[])              # Sum of squared residuals for actual data set
  #fit.new <- sum(sq.pred[])      # Sum of squared residuals for new data set
  # test[j] <- step(fit.new[j]-fit[j])     # Test whether new data set more extreme
  # bpvalue[j] <- mean(test[j])         # Bayesian p-value
  
  
} # END OF MODEL
      
