####################################
############ PROJET 3 ##############
####################################


acci<-read.table(file="/Users/Sarah/Desktop/Master 1/S1/Mathematical statistic/Projet/acci.txt.txt", head=TRUE)
require(stats4)
library(lubridate)
library(plyr)

acci$DATE <- as.Date(acci$DAT, format='%d/%m/%y')
acci$WEEK <- week(acci$DATE)
acci$MONTH <- month(acci$DATE)
acci$YEAR <-year(acci$DATE)

attach(acci)

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


### Question 3

# Parameters
alpha <- 142
betaW <- (212)^(-1)
betaM <- (48)^(-1)

# a) Maximum likelihood

# we suppose data follow a poisson distribution
# in order to obtain the MLE, we need to maximize the log-likelihood function

# per week :
vwv<-as.vector(vw)
length(vwv)
sum(vwv)
mleW <- (1/length(vwv))*sum(vwv)

# per month :
vmv<- as.vector(vm)
length(vmv)
sum(vmv)
mleM <- (1/length(vmv))*sum(vmv)

# b) Bayesian estimator from discrete prior

# per week

library(LearnBayes)
max(vwv)
priorw <- rep(1/100, 100)
names(priorw) <- seq(from=0, to=max(vwv), length.out=100)
postw <- discrete.bayes(dpois, priorw, vwv)
print(postw)
summary(postw)

# per month

max(vmv) 
priorm <- rep(1/100, 100)
names(priorm) <- seq(from=0, to=max(vmv), length.out=100)
postm <- discrete.bayes(dpois, priorm, vmv)
print(postm)
summary(postm)


# c) Bayesian estimator from Gamma prior

# per week
postWG <- (sum(vwv)+alpha)/(length(vwv)+1/betaW)

# per month
postMG <- (sum(vmv)+alpha)/(length(vmv)+1/betaM)


### Question 4

# We plot several density curves in order to compare the different estimation methods used
# In order to have a better comparison of the several distribution used, we combine them onto the same graph

library(ggplot2)

# WEEK

# a) We plot the posterior distribution of theta|X1..Xn corresponding to the gamma a priori
#Post uniform
xfitPUM=seq(min(1),max(5),length=500)
yfitPUM=dnorm(xfitPUM,mean=summary(postm)$mean,sd=summary(postm)$sd)

plot(xfitPUM,yfitPUM,type='l') #post uniform

# b) An a priori Gamma density of theta that we have chosen
# We took the parameters deduced from the dataset
PriGw=rgamma(208, shape=150,rate=220)

xfitPriGw=seq(min(0),max(2),length=500)
yfitPriGw=dnorm(xfitPriGw,mean=mean(PriGw),sd=sd(PriGw))

lines(xfitPriGw,yfitPriGw,type='l',col="red") #prior gamma

# c) We plot the posterior distribution of theta|X1..Xn corresponding to the gamma a priori
#Post Gamma

#Post Gamma
PostGw=rgamma(208, shape=150+142,rate=428)

xfitPostGw=seq(min(0),max(2),length=500)
yfitPostGw=dnorm(xfitPostGw,mean=mean(PostGw),sd=sd(PostGw))

plot(xfitPUW,yfitPUW,type='l', xlim=c(0,1.3), ylim=c(0,10))
lines(xfitPostGw,yfitPostGw,type='l',col="blue") #post gamma


# MONTH

# a) We plot the posterior distribution of theta|X1..Xn corresponding to the gamma a priori
#Post uniform
xfitPUW=seq(min(0),max(2),length=500)
yfitPUW=dnorm(xfitPUW,mean=summary(postw)$mean,sd=summary(postw)$sd)

plot(xfitPUW,yfitPUW,type='l') #post uniform

# b) An a priori Gamma density of theta that we have chosen
# We took the parameters deduced from the dataset

# Prior gamma
PriGm=rgamma(48, shape=150,rate=55)

xfitPriGm=seq(min(1),max(5),length=500)
yfitPriGm=dnorm(xfitPriGm,mean=mean(PriGm),sd=sd(PriGm))

lines(xfitPriGm,yfitPriGm,type='l',col="red") #prior gamma

# c) We plot the posterior distribution of theta|X1..Xn corresponding to the gamma a priori
#Post Gamma

#Post Gamma
PostGm=rgamma(48, shape=150+142,rate=48+55)


xfitPostGm=seq(min(0),max(5),length=500)
yfitPostGm=dnorm(xfitPostGm,mean=mean(PostGm),sd=sd(PostGm))

plot(xfitPUM,yfitPUM,type='l', ylim=c(0,2.7)) #post uniform
lines(xfitPostGm,yfitPostGm,type='l',col="blue") #post gamma


### Question 5

# Now we want to estimate theta (the mean number of plane accidents) with a data-based simulation approach
# We use the bootstrap method
# We draw B = 500 random samples of X1..Xn (week/month) with :
# -> the mean set to the MLE estimates of the question 3
# -> the sample size set to the number of weeks/months in our initial sample

# We finally evaluate the bias and root mean square error (RMSE) of each estimation method
# MLE VS Discrete Prior VS Gamma Prior
# We only use the B random samples, ignoring the initial sample

library(boot)
set.seed(1000)

# WEEK

bootstw=vector()
for(i in 1:500)
  bootstw<-rbind(bootstw,rpois(length(vwv),mleW))
#500 samples with size=208 and mean=0.6826923 (MLE for week)


# Maximum likelihood
tw=apply(bootstw,1,mean) #A vector of MLE from each sample

mean(tw) #MLE mean
# Theta is predicted to be 0.6825096 by the maximum likelihood estimator

biasMLEw= mean(tw)-mleW
biasMLEw
# Low but negative bias (-0.0001826923) -> MLE underestimates theta

# We create a function calculating the root mean square error
rmse <- function(x,t){
  sqrt((1/length(x))*sum((x-t)^2))
}

rmseMLEw=rmse(tw,mleW)
rmseMLEw
# RMSE is low but positive as well (0.05601859)

# Bayesian estimation 1 : using a uniform distribution as the prior
BAYES <- function(data, i){
  d <- data[i]
  priorm <- rep(1/100, 100)
  names(priorm) <- seq(from=0, to=max(d), length.out=100)
  postm <- discrete.bayes(dpois, priorm, d)
  k<-summary(postm)$mean
  
  return(k)
}

B1w=apply(bootstw,1,BAYES) # a vector of theta estimated by the bayesian method with discrete prior

mean(B1w)
# Bayes with discrete prior predicts theta to be 0.6873173

biasB1=mean(B1w)-summary(postw)$mean
biasB1
# low but negative bias (-0.0001826975) as well as MLE

rmseB1w=rmse(B1w,summary(postw)$mean)
rmseB1w
# approximately the same RMSE as MLE

# Bayesian estimation 2 : using a Gamma as the prior
BAYESGmw <- function(data, i){
  d <- data[i]
  n=sum(d)+alpha
  d=length(d)+1/betaW
  
  return(n/d)
}

B2w=apply(bootstw,1,BAYESGmw) # a vector of theta estimated by Bayesian method with gamma prior

mean(B2w)
# theta is estimated to be 0.6821542 here

biasB2w=mean(B2w)-postWG
biasB2w
# lower but negative bias than MLE and Bayes with discrete prior
# this estimator of theta nearly unbiased (less bias than the other at least)

rmseB2w=rmse(B2w,postWG)
rmseB2w
# RMSE is lower as well
# On average, predictions are better with this estimation method
# relevant with our assumptions

# MONTH

# Generate samples
bootstm=vector()
for(i in 1:500)
  bootstm <-rbind(bootstm,rpois(length(vmv),mleM))
# 500 samples with size = 48 and mean=2.958333 (Monthly average of plane accidents over a month estimated by ML method)

# Maximum likelihood
tm=apply(bootstm,1,mean) #A vector of MLE from each sample

mean(tm) #MLE mean
# Theta is predicted to be 2.954208 by the maximum likelihood estimator

biasMLEm= mean(tm)-mleM
biasMLEm
# Low but positive bias (-0.004125) -> MLE underestimates theta

rmseMLEm=rmse(tm,mleM)
rmseMLEm

# Bayesian estimation 1 : using a uniform distribution as the prior

B1m=apply(bootstm,1,BAYES)# a vector of theta estimated by the bayesian method with discrete prior

mean(B1m)
# Bayes with discrete prior predicts theta to be 2.975042

biasB1m=mean(B1m)-summary(postm)$mean
biasB1m
# really low and negative bias (close to 0)

rmseB1m=rmse(B1m,summary(postm)$mean)
rmseB1m
# approximately same RMSE as MLE

# Bayesian estimation 2 : using a Gamma as the prior
BAYESGm <- function(data, i){
  d <- data[i]
  n=sum(d)+alpha
  d=length(vmv)+1/betaM
  
  return(n/d)
}

B2m=apply(bootstm,1,BAYESGm) # a vector of theta estimated by Bayesian method with gamma prior

mean(B2m)
# theta is estimated to be 2.833029 here

biasB2m=mean(B2m)-postMG
biasB2m

rmseB2m=rmse(B2m,postMG)
rmseB2m
# We conclude similarly (as the sample divided into weeks)
# Gamma prior used in a Bayesian method of estimation lead to better results
# Estimator is less biased

