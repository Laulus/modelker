
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
      
      
