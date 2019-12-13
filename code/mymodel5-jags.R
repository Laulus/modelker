
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
      
      
