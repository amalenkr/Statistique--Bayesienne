

acci<-read.table(file="C:/Users/amale/Documents/GitHub/Statistique--Bayesienne/Lauzeral_Nokri_projet/plane_accidents.txt", head=TRUE)

require(stats4)
library(lubridate)
library(plyr)
library(LearnBayes)
library(ggplot2)
library(boot)

acci$DATE <- as.Date(acci$DAT, format='%d/%m/%y')
acci$WEEK <- week(acci$DATE)
acci$MONTH <- month(acci$DATE)
acci$YEAR <-year(acci$DATE)

attach(acci)
set.seed(1000)

x=vector()
vm=vector()
for(i in 1972:1975){
  for(j in 1:12){
    if(length(acci$MONTH[acci$MONTH==j & acci$YEAR==i])>0){
      x[j]=length(acci$MONTH[acci$MONTH==j & acci$YEAR==i])
      
    } else{
      x[j]=0
    }
  }
  vm=rbind(vm,x)
}

t=vector()
vw=vector()
for(i in 1972:1975){
  for(j in 1:53){
    if(length(acci$WEEK[acci$WEEK==j & acci$YEAR==i])>0){
      t[j]=length(acci$WEEK[acci$WEEK==j & acci$YEAR==i])
      
    } else{
      t[j]=0
    }
  }
  vw=rbind(vw,t)
}



### Bayesian estimator from uniform prior

# per week
vwv<-as.vector(vw)
max(vwv)
priorw <- rep(1/100, 100)
names(priorw) <- seq(from=0, to=max(vwv), length.out=100)
postw <- discrete.bayes(dpois, priorw, vwv)
print(postw)
summary(postw)
plot(postw)

# per month
vmv<- as.vector(vm)
max(vmv) 
priorm <- rep(1/100, 100)
names(priorm) <- seq(from=0, to=max(vmv), length.out=100)
postm <- discrete.bayes(dpois, priorm, vmv)
print(postm)
summary(postm)
plot(postm)



### Bayesian estimator from Gamma prior

# Parameters
alpha <- 142
betaW <- (212)^(-1)
betaM <- (48)^(-1)

# per week
postWG <- (sum(vwv)+alpha)/(length(vwv)+1/betaW)

# per month
postMG <- (sum(vmv)+alpha)/(length(vmv)+1/betaM)



### Graphiques per week

# Post uniform
xfitPUW=seq(min(0.3),max(1),length=500)
yfitPUW=dnorm(xfitPUW,mean=summary(postw)$mean,sd=summary(postw)$sd)
plot(xfitPUW,yfitPUW,type='l', xlim=c(0.3,1), ylim=c(0,10))

# Prior gamma
PriGw=rgamma(212, shape=alpha, scale=betaW)
xfitPriGw=seq(min(0),max(2),length=500)
yfitPriGw=dnorm(xfitPriGw,mean=mean(PriGw),sd=sd(PriGw))
lines(xfitPriGw,yfitPriGw,type='l',col="red")

# Posterior Gamma
PostGw=rgamma(212, shape=2*alpha, rate=2*212)
xfitPostGw=seq(min(0),max(2),length=500)
yfitPostGw=dnorm(xfitPostGw,mean=mean(PostGw),sd=sd(PostGw))
lines(xfitPostGw,yfitPostGw,type='l',col="blue")



### Graphiques per month

# Post uniform
xfitPUM=seq(min(1),max(5),length=500)
yfitPUM=dnorm(xfitPUM,mean=summary(postm)$mean,sd=summary(postm)$sd)
plot(xfitPUM,yfitPUM,type='l', ylim=c(0,2.7)) 

# Prior gamma
PriGm=rgamma(48, shape=alpha, scale=betaM)
xfitPriGm=seq(min(1),max(5),length=500)
yfitPriGm=dnorm(xfitPriGm,mean=mean(PriGm),sd=sd(PriGm))
lines(xfitPriGm,yfitPriGm,type='l',col="red")

# Post Gamma
PostGm=rgamma(48, shape=2*alpha, rate=2*48)
xfitPostGm=seq(min(0),max(5),length=500)
yfitPostGm=dnorm(xfitPostGm,mean=mean(PostGm),sd=sd(PostGm))
lines(xfitPostGm,yfitPostGm,type='l',col="blue")


#################################################################

rmse <- function(x,t){
  sqrt((1/length(x))*sum((x-t)^2))
}


# Estimateur bay??sien issue d'une prior uniform
BAYES <- function(data, i){
  d <- data[i]
  priorm <- rep(1/100, 100)
  names(priorm) <- seq(from=0, to=max(d), length.out=100)
  postm <- discrete.bayes(dpois, priorm, d)
  k <- summary(postm)$mean
  return(k)
}


# Estimateur bay??sien par week issue d'une prior gamma
BAYESGmw <- function(data, i){
  d=data[i]
  n=sum(d)+alpha
  d=length(d)+1/betaW
  return(n/d)
}


# Estimateur bay??sien par month issue d'une prior gamma
BAYESGm <- function(data, i){
  d <- data[i]
  n=sum(d)+alpha
  d=length(vmv)+1/betaM
  return(n/d)
}

#################################################################


### Bootstrap per week

# Cr??ation d'??chantillons bootstrap
vwv <- as.vector(vw)
length(vwv)
sum(vwv)
mleW <- (1/length(vwv))*sum(vwv)

bootstw=vector()
for(i in 1:500)
  bootstw<-rbind(bootstw, rpois(length(vwv), mleW))


# Uniform prior bayesian estimator
B1w=apply(bootstw,1,BAYES)
biasB1=mean(B1w)-summary(postw)$mean ; biasB1
rmseB1w=rmse(B1w,summary(postw)$mean) ; rmseB1w


# Gamma prior bayesian estimator
B2w=apply(bootstw,1,BAYESGmw)
biasB2w=mean(B2w)-postWG ; biasB2w
rmseB2w=rmse(B2w,postWG) ; rmseB2w



### Bootstrap per month

# Cr??ation d'??chantillons bootstrap
vmv<- as.vector(vm)
length(vmv)
sum(vmv)
mleM <- (1/length(vmv))*sum(vmv)

bootstm=vector()
for(i in 1:500)
  bootstm <-rbind(bootstm,rpois(length(vmv),mleM))


# Uniform prior bayesian estimator
B1m=apply(bootstm,1,BAYES)
biasB1m=mean(B1m)-summary(postm)$mean ; biasB1m
rmseB1m=rmse(B1m,summary(postm)$mean) ; rmseB1m


# Gamma prior bayesian estimator
B2m=apply(bootstm,1,BAYESGm)
biasB2m=mean(B2m)-postMG ; biasB2m
rmseB2m=rmse(B2m,postMG) ; rmseB2m

