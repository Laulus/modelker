


model {

  # likelihood Backcalculation
   for (i in 1:length(idfish)){
     len[i]~dnorm(mu.tot[i], tau.tot)
     mu.tot[i]<-a.retro + c.retro* RT[i]^b.retro
  } # end of loop i

  # priors
  a.retro~dnorm(30,0.01) #T(0,)
  b.retro~dnorm(1,0.01)
  c.retro~dnorm(3,0.01)


# tau.tot <- 1/sqrt(sigma.tot)
# sigma.tot~dunif(0,10)

  tau.tot~dgamma(0.001,0.001)

  #predictions
  # WARNING RM=RT for sedentary individuals !!!!!!!!!!!!!!!!!

## DANS CETTE BOUCLE LA IL FAUT RAJOUTER DE LA VARIANCE SUR LES RADI,
# et ajouter de la variance sur les RT
# admettons qui augmente avec l'age j
# on aura donc une len retro i j moyenne et des intervalles de confiance



   for (i in 1:N) { # for each lines
         # len.retro[i] <- a.retro + ((length[i]*RM[i]/rt[i])-a.retro)*(((radi[i]+eps[aget[i]])/RM[i])^b.retro) # attention a la matrice de depart pour les r[i,j]
len.retro[i]~dnorm(mu.retro[i], tau.retro)
mu.retro[i]<-a.retro + ((length[i]*RM[i]/rt[i])-a.retro)*(((radi[i])/RM[i])^b.retro) # attention a la matrice de depart pour les r[i,j]  
   }

# tau.retro <- 1/sqrt(sigma.retro)
# sigma.retro~dunif(0,10)
tau.retro~dgamma(0.001,0.001)

# for (j in 1:max(aget[])){
# eps[j]~dnorm(0, 0.001)
# }


} # END OF MODEL

  
