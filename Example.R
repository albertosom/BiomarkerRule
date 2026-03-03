
## Example ##

rm(list = ls())

library(caret)
library(randomForest)
library(tidyverse)
library(nloptr)
library(ranger)
library(SuperLearner)
library(parallel)

source("PluginRule.R")
source("Standard.R")
source("IT-DOOLR.R")
source("DOOLR.R")

ncores<- detectCores() - 1 # to prevent 100% CPU usage

## Generate a cohort type data set with two biomarkers and prevalence of 0.01

beta0=-8.700
beta<-c(2.4,2.4)
prev=0.01
alpha=0.04

logistic<-function(x){
  1/(1+exp(-x))
}

data.gen<-function(n,beta0,beta,p=0.06) {
  npfloor<-floor(n*p)
  X1<- rnorm(n-npfloor,mean=0,sd=1)
  W<-rnorm(n-npfloor,mean=0,sd=1)
  p.true<-logistic(beta0 +X1*beta[1] +W*beta[2])
  D<-rbinom(n-npfloor,size=1,p=p.true)
  df1<-data.frame(D,X1,W)
  if(p==0){
    return(df1)
  }
  else{
    df2<- data.frame(D=0,X1=rep(6,npfloor),W=rep(6,npfloor))
    return(rbind(df1,df2))
  }
  
}

## data generation ##
data.test<-data.gen(n=1000000,beta0,beta,p=0.06)
data.train<-data.gen(n=1250,beta0,beta,p=0.06)
data.val<-data.gen(n=1250,beta0,beta,p=0.06)
data.dt<-rbind(data.train,data.val)

## Standard method
thres<-seq(0,0.1,length.out =500)
stand.rslt<-standard.lg(data.train=data.dt,data.test=data.test,prev=NULL,alpha=alpha,KK=2,thres,
                      data.type="cohort")

print(stand.rslt)

##DOOLR approach

n.train<-nrow(data.train)
D.train<-data.train[,1]
n1<-sum(D.train==1)
n0<-sum(D.train==0)
X.train<-subset(data.train,select=2:ncol(data.train))
X.train<-data.matrix(X.train)

beta.int<-stand.rslt$beta/norm(as.matrix(stand.rslt$beta),"F")
hval = sd(as.matrix(cbind(1,X.train)) %*% matrix(beta.int,ncol=1))/(nrow(data.train)^(1/3))
#hval = sd(as.matrix(cbind(1,X.train)) %*% matrix(beta.null,ncol=1))/(nrow(data.train)^(1/3))
h<-hval

doolr.rslt<-doolr.rule(beta.int=beta.int,data.train=data.train,
                           data.val=data.val,
                           data.test=data.test,prev=NULL,h=1,alpha=alpha
                           ,data.type="cohort")

print(doolr.rslt)

## Plug-in approach

plugin.rslt<-plugin.rule(data.train=data.train,data.val,data.test,prev=NULL,alpha=alpha,
                data.type="cohort")

print(plugin.rslt)

####################################################
#                                                  #
#    Data generation with external information     #
#                                                  #
####################################################


nonlinear_logistic <- function(x) {
  1 / (1 + exp(-x))
}

# Data-generating function

data.gen <- function(n = 5000,beta0 = -8.5,beta = c(5, - 4, 3)) {
  X1 <- rnorm(n)
  X2 <- rnorm(n)
  X3 <- rnorm(n)
  
  # Nonlinear function of biomarkers
  f <- beta0 + beta[1]*sin(X1) + beta[2] * (X2^2) + beta[3]*cos(X3)
  
  # Probability of D = 1
  p <- nonlinear_logistic(f)
  # Binary outcome
  
  err<- rlogis(n, location = 0, scale = 1)
  risk.score=f+err
  D<-ifelse(risk.score >0,1,0)
  
  df<-data.frame(D, X1, X2, X3)
  return(list("data"=df,"error"=err))
}

beta0 = -8.6; beta = c(5, - 4, 3)
alpha=0.04
betas<-c(beta0,beta)
p<-length(betas)
ncores=ncores<- detectCores() - 1 # to prevent 100% CPU usage
thres<-seq(0,0.1,length.out =500)


  data_train<- data.gen(n=1250,beta0 =beta0,beta =beta)
  data.train<-data_train$data
  error.train<- data_train$error
  data_val<- data.gen(n=1250,beta0 =beta0,beta =beta)
  data.val<-data_val$data
  error.val<- data_val$error
  data_test<- data.gen(n=1e6,beta0 =beta0,beta =beta)
  data.test<- data_test$data
  error.test<- data_test$error

  data.dt<-rbind(data.train,data.val)
  
  ## Standard approach
  
  stand.rslt<-standard.lg(data.train=data.dt,data.test=data.test,prev=NULL,alpha=alpha,KK=2,thres,
                        data.type="cohort")
  
  print(stand.rslt$rr)
  
  ## Doolr approach 
  
  beta.est.st<-stand.rslt$beta
 
  n.train<-nrow(data.dt)
  D.train<-data.dt[,1]
  X.train<-subset(data.dt,select=2:ncol(data.dt))
  X.train<-data.matrix(X.train)
  
  beta.int<-beta.est.st/norm(as.matrix(beta.est.st),"F")
  hval = sd(as.matrix(cbind(1,X.train)) %*% matrix(beta.int,ncol=1))/(n.train^(1/3))
  h<-hval
  
  doolr.rslt<-doolr.rule(beta.int=beta.int,data.train=data.train,
                         data.val=data.val,
                         data.test=data.test,prev=NULL,h=1,alpha=alpha
                         ,data.type="cohort")
  
  print(doolr.rslt)
  
## External model is the same as the true model
  
X<-data.train[,-1]

## external decisions obtained from risk scores.
risk.score<- ifelse(beta0 + beta[1]*sin(X$X1) + beta[2]*(X$X2^2) + beta[3]*cos(X$X3) +error.train>0,1,-1)
  
doolr.rule.IT.CV.rslt<-doolr.rule.IT.CV(beta.int=beta.int,data.train=data.train,
                                    data.val=data.val,
                                    data.test=data.test,prev=NULL,alpha=alpha,h=1,
                                    data.type="cohort",risk_yes_no=risk.score,etas=c(0,0.0001,0.1,0.5,1,5))
    
print(doolr.rule.IT.CV.rslt)

## Plug-in approach

plugin.rslt<-plugin.rule(data.train=data.train,data.val,data.test,prev=NULL,alpha=alpha,
                         data.type="cohort")

print(plugin.rslt)




