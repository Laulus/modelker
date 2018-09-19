

# len should only be size in freshwater
# len can either be directly growth increments => growthincrements or backcalculated sizes => len
      

model {
  
  # likelihood
  
  for (i in 1:N)  {
    len[i] ~ dnorm(mu[i],tau)
    mu[i] <- Linf[i]*(1-exp(-k[i]*(X[i]-t0)))
    # Linf[i]<- a[river[i]] + c[ID[i]] # effet Linf par rivi?re et aussi effet random fish sur Linf
    Linf[i]<- c[ID[i]] # on retire l'effet rivi?re sur Linf
    k[i]<-b[river[i]] + d[ID[i]]
    
  }
  
  # some probas
  # a12<-step(a[1]-a[2])
  # a13<-step(a[1]-a[3])
  # a23<-step(a[2]-a[3])
  b12<-step(b[1]-b[2])
  b13<-step(b[1]-b[3])
  b23<-step(b[2]-b[3])
  
  #predictions
  for (i in 1:3) {
    for (j in 1:3) {
      for (h in 1:9) {
        frypred[i,j,h] ~ dnorm(mupred[i,j,h],tau)
        # mupred[i,j,h] <- a[i]*(1-exp(-b[j]*(h-t0)))
        mupred[i,j,h] <- (1-exp(-b[j]*(h-t0)))
      }
    }   
  }
  
  # priors
  for (i in idfish) {
    c[i]~dnorm(0,tauki)
    d[i]~dnorm(0,tauli)       
  }
  for (i in idriver) {   
    # a[i]~ dunif(50,1500)
    b[i]~ dnorm(0,0.001)       
  }
  
  tau~dgamma(0.001,0.001)
  tauki~dgamma(0.001,0.001)
  tauli~dgamma(0.001,0.001)
  
  t0~ dunif(-10,0.5)
  #a~ dunif(80,1200)
}

      
