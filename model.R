rm(list=ls())   # Clear memory



## WORKING DIRECTORY

work.dir<-paste('G:/ANALYSES DONNEES/Norme_reaction',sep="")

setwd(work.dir)
library(rjags)
load.module("glm")



#--------------- DATA --------------#
library(readr)
bdd_norm <- read_delim("G:/ANALYSES DONNEES/Norme_reaction/bdd_norm.csv", 
                       ";", escape_double = FALSE, trim_ws = TRUE)
bdd_norm<-subset(bdd_norm, Pheno_ec!="TE") # retirer les TE de l'analyse
bdd_norm<-bdd_norm[, c("Fish","River","Scale","X","Fry","Milieu")]
#bdd_norm<-subset(bdd_norm,River!="M") # without MANCHOT RIVER
#bdd_norm<-subset(bdd_norm,Len>0)
#bdd_norm<-droplevels(bdd_norm) # don't forget to drop the level 
bdd_norm$R<-bdd_norm$River
# bdd_norm$River<-as.factor(bdd_norm$River)
# bdd_norm$Fish<-as.factor(bdd_norm$Fish)
# bdd_norm$X<-as.numeric(bdd_norm$X)
# bdd_norm$Scale<-as.factor(bdd_norm$Scale)
# bdd_norm$X<-as.numeric(bdd_norm$X)
bdd_norm<-subset(bdd_norm,bdd_norm$X<11) # try to remove the 6 scales of age 11 and the one of age 10


# for (a in 1 : length(bdd_norm$River)){
#   if (bdd_norm$River[a]=="M")
#     bdd_norm$R[a]<-1
#   else
#     if (bdd_norm$River[a]=="N")
#       bdd_norm$R[a]<-2
#     else 
#       bdd_norm$R[a]<-3
# } # end of loop a, to change River into factor



bdd_norm<-na.omit(bdd_norm)
View(bdd_norm)

N=nrow(bdd_norm)




## DATA

data<-list(N=N, Y=bdd_norm$Milieu, len=bdd_norm$Fry, age=bdd_norm$X, fish=bdd_norm$Fish, scale=bdd_norm$Scale,river=as.factor(bdd_norm$R), idriver=unique(bdd_norm$R), idfish=unique(bdd_norm$Fish))


data$Y=ifelse(data$Y>0,1,0) # poissons TE et TM sont considérés comme migrants
# data$Y=ifelse(data$Y<=1,1,0) # poissons TE et TS sont considérés comme sédentaires

##-----------------------------MODEL ----------------------------------##

write("
      
      
      
      model {
      
      
      
      for (i in 1 : N) {
      
      
      
      # Vraisemblance

Y[i] ~ dbern(p[i])      
# logit(p[i]) <- (beta0[age[i]] + eps.river[river[i]]) + (beta1[age[i]]+betar[river[i]])*(len[i]-mean(len[]))  # MODELE COMPLET : effet rivière / intercept et pente
# logit(p[i]) <- beta0[age[i]] + (beta1[age[i]]+betar[river[i]])*(len[i]-mean(len[])) # effet rivière pente
 logit(p[i]) <- beta0[age[i]] + beta1[age[i]]*(len[i]-mean(len[])) + eps.river[river[i]] # modele avec effet rivière random sur intercept
 # logit(p[i]) <- beta0[age[i]] + beta1[age[i]]*(len[i]-mean(len[])) # modele 0, nor slope nor intercept effect
      
      } # end loop i
      
      
      
      # Prior
      
for (j in 1 : max(age[])) {

      beta1[j] ~ dnorm(0,0.001)
      beta0[j] ~ dnorm(0,0.001)
      # beta0[j] ~ dnorm(mu.beta0,tau.beta0) # common law between age classes
      # beta1[j] ~ dnorm(mu.beta1,tau.beta1)

} # end of loop j 

# mu.beta0~dnorm(0,0.001)
# tau.beta0<-1/(pow(sigma.beta0,2))
# sigma.beta0~dunif(0,100)

#mu.beta1~dnorm(0,0.001)
#tau.beta1<-1/(pow(sigma.beta1,2))
#sigma.beta1~dunif(0,100)


# # METHODE DES CONTRASTES, on fixe les effets pour la rivière n°1 et on calcule les différences par rapport à cette rivière
# eps.river[1]<-0
# eps.river[2]~ dnorm(0,0.001)
# eps.river[3]~ dnorm(0,0.001)
# betar[1] <-0
# betar[2] ~ dnorm(0,0.001)
# betar[3] ~ dnorm(0,0.001)


# CALCUL EFFETS ALEATOIRES PAR RIVIERE
 for (h in idriver){
   eps.river[h]~ dnorm(0,tau.river)
   # betar[h] ~ dnorm(0,tau.betar)
 } # end of loop h

 # tau.river~dgamma(0.1,0.1)
 # sigma<-sqrt(1/tau.river)
 # tau.betar~dgamma(0.1,0.1)
 # sigmar<-sqrt(1/tau.betar)

 tau.river<-1/(pow(sigma,2))
 sigma~dunif(0,50)
 # tau.betar<-1/(pow(sigmar,2))
 # sigmar~dunif(0,10)

# # PROBAS DES ECARTS, METHODE DES CONTRASTES
# 
#   eps.river23<-eps.river[2]-eps.river[3]
#   b23<-step(eps.river23)
#   b12<-step(eps.river[2])
#   b13<-step(eps.river[3])
# 
  # betar23<-betar[2]-betar[3]
  # br23<-step(betar23)
  # br12<-step(betar[2])
  # br13<-step(betar[3])

# PROBAS DES ECARTS, effets aleatoires par rivière

  eps.river12<-eps.river[1]-eps.river[2]
  eps.river13<-eps.river[1]-eps.river[3]
  eps.river23<-eps.river[2]-eps.river[3]
  p23<-step(eps.river23)
  p12<-step(eps.river12)
  p13<-step(eps.river13)

  # betar12<-betar[1]-betar[2]
  # betar13<-betar[1]-betar[3]
  # betar23<-betar[2]-betar[3]
  # br23<-step(betar23)
  # br12<-step(betar12)
  # br13<-step(betar13)
#       
      } # END OF THE MODEL
      
      
      
      ", "mymodel.R")





#------------------------------------------------------------------------------

# ANALYSIS

require(rjags)



nChains = 2 # Number of chains to run.

adaptSteps = 1000 # Number of steps to "tune" the samplers.

burnInSteps = 5000 # Number of steps to "burn-in" the samplers.

niter=50000 # Total number of steps in chains to save.





## PARAMETERS TO SAVE

# parameters to save for CONTRASTS METHOD
# parameters=c("beta0","beta1","eps.river", "betar", "b12","b13","b23","eps.river23", "betar23","br12","br13","br23") 

# parameters to save when considering eps.river, betar randomly
 parameters=c("beta0","beta1","eps.river", "betar", "p12","p13","p23","eps.river12", "eps.river13","eps.river23","sigma")  #"betar12", "betar13", "betar23","br12","br13","br23", "sigmar")

#------------------------------------------------------------------------------

## INITS

inits<-function(){
  # list(beta0=seq(0,-10, length=max(data$age)) , beta1=rep(0.01, max(data$age)))
# inits contrats methods  
  # list(beta0=seq(0,-10, length=max(data$age)) , beta1=rep(0.01, max(data$age))) #, sigma=1)
# inits "random" method 
   list(beta0=seq(0,-10, length=max(data$age)) , beta1=rep(0.01, max(data$age)), sigma=1) #,sigmar=1) # modele avec beta distribution normale (0, 0.001)
  #list(mu.beta0=-5, sigma.beta0=1, mu.beta1=0.01, sigma.beta1=0.001, sigma=1 ) # définition des paramètres des distributions normales de beta0 et beta1
  # list( mu.beta0=-5, sigma.beta0=1 , beta1=rep(0.01, max(data$age)), sigma=1)
}





## Compile & adapt

### Start of the run ###

## Compile & adapt

# Create, initialize, and adapt the model:

jagsfit=jags.model(
  
  'mymodel.R',
  
  data,inits,
  
  n.chains=nChains,
  
  n.adapt = adaptSteps)



# Burn-in:

cat( "Burning in the MCMC chain.\n" )

#jagsfit$recompile()

update(jagsfit, n.iter=burnInSteps,progress.bar="text")



# The saved MCMC chain:

cat( "Sampling final MCMC chain.\n" )

fit.mcmc<-coda.samples(jagsfit,variable.names=parameters, n.iter=niter, thin=1)





## BACKUP

#save(data,jagsfit,fit.mcmc,file=paste("Results_",model.name,".RData",sep=""))



#------------------------------------------------------------------------------

# EXAMINE THE RESULTS



summary(fit.mcmc)
gelman.diag(fit.mcmc,multivariate=F) # test stat de mélange chaine mcmc, convergence? multivariate psrf <1.1


library(mcmcplots)
traplot(fit.mcmc,"beta0")
traplot(fit.mcmc,"beta1")
traplot(fit.mcmc,"sigma")
traplot(fit.mcmc,"sigmar")
traplot(fit.mcmc,"eps.river")
traplot(fit.mcmc,"betar")
caterplot(fit.mcmc,c("beta0"),greek=TRUE, reorder=F)
caterplot(fit.mcmc,c("beta1"),greek=TRUE, reorder=F)
caterplot(fit.mcmc,c("eps.river"),greek=TRUE, reorder=F)
# caterplot(fit.mcmc,c("betar"),greek=TRUE, reorder=F)
#plot(fit.mcmc)

#------------ COEFFICIENTS --------------#
beta0.med<-summary(fit.mcmc)$quantiles[paste0("beta0[",1:10,"]"),"50%"] # on récupère le vecteur des beta0
beta1.med<-summary(fit.mcmc)$quantiles[paste0("beta1[",1:10,"]"),"50%"] # on récupère le vecteur des beta1
betar<-summary(fit.mcmc)$quantiles[paste0("betar[",1:3,"]"),"50%"]
eps.river<-summary(fit.mcmc)$quantiles[paste0("eps.river[",1:3,"]"),"50%"] # on récupère le vecteur des eps.river

# deltaeps<-summary(fit.mcmc)$statistics[1:3,1:2]
#deltaeps<-summary(fit.mcmc)$statistics[40:42,1:2] # faux a modifier
#-------------PREDICT--------------#

 f<-function(b0,b1,p,e){
#f<-function(b0,b1,p,ri,e){
   Lp<-(log(p/(1-p))-b0-e)/b1
  # Lp<-(log(p/(1-p))-b0)/(b1+ri)
   #Lp<-(log(p/(1-p))-b0-e)/(b1+ri) # MODELE COMPLET
  # Lp<-(log(p/(1-p))-b0)/(b1) # MODELE NUL
  return(Lp)
}

c<-c(0.05,0.25,0.5,0.75,0.95)
M<-list()
for (l in 1:nlevels(data$river)){
L<-matrix(0,5,length(beta0.med))
  for (i in 1:5){
    for (j in 1 : length(beta0.med)){
   # L[i,j]<-f(beta0.med[j],beta1.med[j],p=c[i],betar[l],eps.river[l])+mean(data$len)
       L[i,j]<-f(beta0.med[j],beta1.med[j],p=c[i],eps.river[l])+mean(data$len)
      # L[i,j]<-f(beta0.med[j],beta1.med[j],p=c[i],betar[l])+mean(data$len)
      # L[i,j]<-f(beta0.med[j],beta1.med[j],p=c[i])+mean(data$len)
    M<-append(M,L[i,j])
    } # end of loop for (j in 1 : max(data$age))
  } # end of loop i
} #end of loop l

output<-matrix(M,nrow=15,ncol=10,byrow=T) # nrow=5probas * 3 rivières
o<-as.data.frame.vector(t(output))
o<-cbind(rep(0,length(o)),rep(0,length(o)),rep(0,length(o)),o)
colnames(o)<-c("River","p","Age","Length")
o$Age<-rep(c(1:max(data$age)),15)
o$p<-rep(c(rep(c[1],10),rep(c[2],10),rep(c[3],10),rep(c[4],10),rep(c[5],10)),3)
o$River<-c(rep("Manchot",50),rep("Norvegienne",50),rep("Rohan",50))
o$p<-as.factor(o$p)
o$River<-as.factor(o$River)
o$Length<-as.numeric(o$Length)
o$Age<-as.factor(o$Age)

 
#res<-cbind(as.vector(col(t(L))), as.vector(row(t(L))), as.vector(t(L)))
#colnames(res)<-c("p","Age","Length")
#res<-as.data.frame(res)
#res$p<-as.factor(res$p)
library(ggplot2)

oR<-subset(o,o$River=="Rohan")
oM<-subset(o,o$River=="Manchot")
oN<-subset(o,o$River=="Norvegienne")

#  p0<-ggplot()+
#    geom_line(data=oM, aes(x=Age, y= Length))+ # col=p))+
#     geom_line(data=oN, aes(x=Age, y= Length) )+ #, col=p))+
#     geom_line(data=oR, aes(x=Age, y= Length)) # col=p))
# #   geom_smooth(se=F)#+
# #   #geom_point(data=subset(bdd_norm,bdd_norm$Milieu==0), aes(x=Age, y=Len_Fry), col="blue")+
# #   #geom_point(data=subset(bdd_norm,bdd_norm$Milieu==1), aes(x=Age, y=Len_Fry),col="yellow")+
# #   #geom_point(data=subset(bdd_norm,bdd_norm$Milieu==2), aes(x=Age, y=Len_Fry),col="red")
# p0

#row.names(L)<-c("0.25","0.5","0.75")
#colnames(L)<-c("1","2","3","4","5","6","7","8","9","10","11")

# plott<-ggplot(o,aes(Age,Length))+
#   + geom_point(aes(colour=River))+
#   + geom_smooth(aes(group=p))

# #ylim=c(100,550)
# p <- ggplot(data = o, aes(x = Age, y = Length),ylim=c(100,550)) +
# #     geom_line(data = oM, aes(x = Age, y = Length))+
# #     geom_line(data = oN, aes(x = Age, y = Length))+
# geom_point() +
# facet_wrap(~River, ncol = 3) +
# geom_smooth(aes(group=p,col=p)) +
#    theme_bw()
#  p

plot(NULL, ylim=c(0,500), xlim=c(1,10),xlab="Age (Year)",ylab="Length (mm)")
points(1:10,subset(oM,p=="0.5")$Length,pch=16,col="lightblue")
lines(1:10,subset(oM,p=="0.5")$Length,pch=16,col="lightblue")
segments(1:10,subset(oM,p=="0.25")$Length,1:10,subset(oM,p=="0.75")$Length,lwd = 2,col="lightblue")
segments(1:10,subset(oM,p=="0.05")$Length,1:10,subset(oM,p=="0.95")$Length,lwd=0.5,lty=2,col="lightblue")
points(c(1:10)+0.1,subset(oN,p=="0.5")$Length,pch=16,col="orange")
lines(c(1:10)+0.1,subset(oN,p=="0.5")$Length,pch=16,col="orange")
segments(c(1:10)+0.1,subset(oN,p=="0.25")$Length,c(1:10)+0.1,subset(oN,p=="0.75")$Length,lwd = 2,col="orange")
segments(c(1:10)+0.1,subset(oN,p=="0.05")$Length,c(1:10)+0.1,subset(oN,p=="0.95")$Length,col="orange",lwd=0.5,lty=2)
# points(c(1:6)+0.2,subset(oR,p=="0.5")$Length[1:6],pch=16,col="darkblue")
# lines(c(1:6)+0.2,subset(oR,p=="0.5")$Length[1:6],pch=16,col="darkblue")
# segments(c(1:6)+0.2,subset(oR,p=="0.25")$Length[1:6],c(1:6)+0.2,subset(oR,p=="0.75")$Length[1:6],lwd = 2,col="darkblue")
# segments(c(1:6)+0.2,subset(oR,p=="0.05")$Length[1:6],c(1:6)+0.2,subset(oR,p=="0.95")$Length[1:6],lwd=0.5,lty=2,col="darkblue")
legend(x=7,y=150,legend = c("MANCHOT","NORVEGIENNE","ROHAN"),lty=c(1,1,1), col=c("lightblue","orange","darkblue"), title = "p(Migration) = 50%")


gap <- c(-.1, 0, .1)
plot(NULL, ylim=c(0,700), xlim=c(1,10),xlab="Age", ylab="Length")
points(data$age+gap[data$river], data$len,col=data$river,pch=data$Y+2)

#--------------BOXPLOT DES ECARTS -------------#
deltaeps2<-t(deltaeps)
boxplot(deltaeps2)
