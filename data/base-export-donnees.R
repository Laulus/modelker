bddmodelo<-read.csv("C:/Users/laulus/Documents/Mirror/THESE/BDD/BDD-THESE-TRAVAIL.csv",sep=";")


setwd(dir = "C:/Users/laulus/Documents/Mirror/THESE/MODELE-COMPLET/modelker/data")
# We want to select ony fish that have been read
# the problem is that radius were not measured the same according to the period
# so we can basically split the dataset in two parts
# the reader X from which the measure of radii have to be converted
# all the other reader that did measure radii in micrometers

levels(bddmodelo$lecteur) # here you can see the name of the different readers
# bddmodeloX<-subset(bddmodelo,bddmodelo$lecteur=="X") # dim : 3154 fish (98 variables or columns recorded)
# bddmodeloo<-subset(bddmodelo,bddmodelo$lecteur!="" & bddmodelo$lecteur!="X") # dim: 1260 (lines read) 98 (variables or column recorded)

# we should convert the value in bddmodeloX to add them in bdd
# HERE ADD EDDY's DATA

# bdd<-bddmodeloo # to change hereafter bddmodeloo + bddmodeloX converted
bdd<-bddmodelo
# bdd<-subset(bdd, bdd$rq.lecteur=="")
bdd<-cbind(bdd[,c("id","annee","mois","Rivière","pheno","Sexe","lf","Poid.g","lecteur","ecaille.et","age.tot","age.sed")],bdd[,58:75])
# View(bdd)
bdd[,14:27]<-apply(bdd[,14:27],2, as.numeric)
bdd$mois<-as.numeric(bdd$mois)
bdd<-subset(bdd,mois!="NA")
# garder quand m?me les ?cailles et le lecteur car pour filtrer je vais en avoir besoin

# On ne garde que les poissons lus et mesur?s, ie dans la colonne lecteur il y a une nom
bdd<-subset(bdd,lecteur!="")
# bdd<-subset(bdd,lecteur!="X")
# et on ne garde que les TM
# bdd<-subset(bdd,pheno=="TM")
# probleme, il y a des age.tot que l'on ne connait pas
# il faut garder que les poissons o? on a une age.tot d?termin?
bdd<-subset(bdd,age.tot!="NA")
bdd<-subset(bdd,age.sed!="NA")
bdd<-subset(bdd,age.sed!=0)
# probleme il y a des lf manquants
bdd<-subset(bdd,bdd$lf!="NA")
# probl?me il y a des RM manquants
# bdd<-subset(bdd,RM!="NA")
# probl?me il y a des RT manquants
bdd<-subset(bdd,RT!="NA")
# probl?me il y a des ?cailles r?g?n?r?es R1
bdd<-subset(bdd,R1!="NA")
# probl?me il y a des agetot ou on a pas le radius ...

ageent<-rep(0,length=length(bdd$age.tot))
for (i in 1:length(bdd$age.tot)){
  if (bdd$mois[i]<7){ageent[i]<-bdd$age.tot[i]}
  else {ageent[i]<-bdd$age.tot[i]-1}
}

# vector<-rep(NA,nrow(bdd))
# for (i in 1:nrow(bdd)){
#   if (is.na(bdd[i,bdd$age.tot[i]+15])==TRUE){vector[i]<-i}
# }
# vivi<-na.omit(vector)
# bdd<-bdd[-vivi,]

bddi<-bdd # on enregistre la base initiale
bddi<-cbind(bddi,ageent)
bddi<-subset(bddi, ageent!="0") # on retire les 0+

# SAMPLE OR NOT
# bdd<-bddi
# bdd<-bddi[sample(1:nrow(bddi),200,replace=FALSE),] # on fait un sample de 200 lignes au hasard
# bdd<-bddi[sample(1:nrow(bddi),500,replace=FALSE),] 

# dimensions
nID<-nlevels(as.factor(bdd$id))
nScale<-nrow(bdd)
maxAge<-max(bdd$age.tot)

# definitions of variables
ID<-bdd$id
age.tot<-bdd$age.tot 
# age<-bdd$ageent

scaleMaxSize<-bdd$RT
# length(ID)
# length(scaleMaxSize)
# length(age)


radius<-matrix(data=NA, nrow=nrow(bdd), ncol=maxAge)# matrix (number of lines in the datasets (ie number of scales), maximum radius read R15)
for (i in 1:nrow(bdd)){
  for (j in 1:maxAge){
    radius[i,j]<-bdd[i,15+j]
  } # end of loop j
} # end of loop i


age<-rep(0,length=length(bdd$age.tot))
for (i in 1:nrow(bdd)){
age[i]<-sum(!is.na(radius[i,])==TRUE)}

smoltAge<-bdd$age.sed


# il y a un probleme dans les smolt Age a r?gler ... tel que coder le poissson est consid?r?
# comme dispersant r?ellement l'ann?e suivant son ?ge ? la migration (car RM est compris le plus souvent entre 2 rayons)
smolt<-matrix(data=NA, nrow=nrow(bdd), ncol=maxAge) # matrix (number of lines in the datasets (ie number of scales), maximum radius read R15)
RM<-bdd$RM
for (i in 1:nrow(bdd)){
  for (j in 1:maxAge){
    if (is.na(RM[i])==TRUE){smolt[i,j]<-0}
    if (is.na(RM[i])==FALSE & radius[i,j]<RM[i] & is.na(radius[i,j])==FALSE){smolt[i,j]<-0}
    if (is.na(RM[i])==FALSE & is.na(radius[i,j])==TRUE){smolt[i,j]<-smolt[i,j]}
    if (is.na(RM[i])==FALSE & radius[i,j]>=RM[i] & is.na(radius[i,j])==FALSE){smolt[i,j]<-1}
  } # end of loop j
} # end of loop i

fishMaxSize<-bdd$lf


# sur body size a va etre le m?me probl?me que sur smolt
bodySize<-matrix(data=NA, nrow=nrow(bdd), ncol=maxAge)
age.tot<-bdd$age.tot
for (i in 1:nrow(bdd)){
  for (j in 1:maxAge){
    if (j==age[i]){bodySize[i,j]<-bdd$lf[i]}
    else {bodySize[i,j]<-bodySize[i,j]}
  } # end of loop j
} # end of loop i



# EXPORT DATA for OPENBUGS
data<- list ("ID", "age","age.tot", "scaleMaxSize","radius","smoltAge","smolt","fishMaxSize","bodySize","nID","nScale","maxAge")
library(R2OpenBUGS)
bugs.data(data,dir=getwd(),digits=5,data.file="data-all.txt")
# bugs.data(data,dir=getwd(),digits=5,data.file=paste("data-",k,".txt"))
# } # end of loop k
