
      model {
      
     for (i in 1:nScale) { # 20 scales belonging to 10 individuals, unbalanced. 
     
           Y[i]<- ifelse(age[i]>smoltAge[i],2,1) # 2: at sea / 1: in freshwater
           
		for (j in 1:age[i]) { # maxAge of each scale
      
      ### likelihood
      # 1. retrocalcul
      bodySize[i,j] ~ dlnorm(allo[i,j],tau[1]) #T(,fishMaxSize[i]+1)
      allo[i,j] <- a[Y[i]] + X[i,j,1]* b[Y[i]] * (fishMaxSize[i]*X[i,j,2] -a[Y[i]])
      
      X[i,j,1] <- radius[i,j]/radius[i,smoltAge[i]]
      X[i,j,2] <- radius[i,smoltAge[i]]/scaleMaxSize[i]
      #X[i,j,1] <- ifelse(age[i]>smoltAge[i], radius[i,j]/radius[i,smoltAge[i]], radius[i,j]/radius[i,age[i]])
      #X[i,j,2] <- ifelse(age[i]>smoltAge[i], radius[i,smoltAge[i]]/scaleMaxSize[i], 1)

      } # END OF LOOP j
      } # END OF LOOP i
      
      ### priors
      #a[1]~dnorm(0, 0.01) 
      a[1]~dunif(0, 10) 
      #a[2]~dunif(0, 10)
      a[2] <- a[1] #+ delta[1]
      delta[1]~dnorm(0, 0.01) 

      #b[1]~dnorm(0, 0.01) 
      b[1]~dunif(0, 10) 
      b[2]~dunif(0, 10) 
      #b[2] <- b[1] #+ delta[2]
      delta[2]~dnorm(0, 0.01) 
      
      tau[1]<-1/(sigma[1]^2)
      sigma[1]~dunif(0,10)
      
      } # END OF MODEL
      
      
