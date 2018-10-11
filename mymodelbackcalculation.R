

data {  for (i in 1:N) { # for each lines
         len.retro[i] <- a.retro + ((length[i]*RM[i]/rt[i])-a.retro)*((radi[i]/RM[i])^b.retro) # attention a la matrice de depart pour les r[i,j]
}
}

model {

  # likelihood Backcalculation
   for (i in 1:50){
     len[i]~dnorm(mu.tot[i], tau.tot)
     mu.tot[i]<-a.retro + c.retro* RT[i]^b.retro
  } # end of loop i

  # priors
  a.retro~dnorm(100,0.01) #T(0,)
  b.retro~dnorm(1,0.01)
  c.retro~dnorm(0.1,0.01)

  tau.tot~dgamma(0.001,0.001)
# sigma.tot<-1/tau.tot

  #predictions
  # WARNING RM=RT for sedentary individuals !!!!!!!!!!!!!!!!!

## DANS CETTE BOUCLE LA IL FAUT RAJOUTER DE LA VARIANCE SUR LES RADI,
# et ajouter de la variance sur les RT
# admettons qui augmente avec l'age j
# on aura donc une len retro i j moyenne et des intervalles de confiance



  # for (i in 1:N) { # for each lines
  #       len.retro[i] <- a.retro + ((length[i]*RM[i]/rt[i])-a.retro)*(((radi[i]+eps[aget[i]])/RM[i])^b.retro) # attention a la matrice de depart pour les r[i,j]
  # }

# for (j in 1:max(aget[])){
# eps[j]~dnorm(0, 0.001)
# }


# --- GROWTH --- #
  for (i in 1:N)  { #for each lines
    len.retro[i] ~ dnorm(mu[i],tau)
    mu[i] <- Linf*(1-exp(-k*(agei[i]-t0)))
    # Linf[i]<- a[river[i]] + c[ID[i]] # effet Linf par riviere et aussi effet random fish sur Linf
    # Linf[i]<- c[ID[i]] # on retire l'effet rivi?re sur Linf
    # k[i]<-b[river[i]] + d[ID[i]]
    # k[i]<-d[ID[i]]

  }

# # priors
#   for (i in idfish) {
#     # c[i]~dnorm(0,tauki)
#     d[i]~dnorm(0,tauli)
#   }

k~dnorm(0,0.01)
Linf~dnorm(1000, 0.001)
tau~dgamma(0.001,0.001)
tauki~dgamma(0.001,0.001)
tauli~dgamma(0.001,0.001)
t0~ dunif(-10,0.5)

} # END OF MODEL
  
