
## Example ##

rm(list = ls())

library(caret)
library(randomForest)
library(tidyverse)
library(nloptr)
library(ranger)
library(SuperLearner)
library(parallel)

source("pluginRule.R")
source("standard.R")
source("IT-DOOLR.R")
source("DOOLR.R")

ncores=ncores<- detectCores() - 1 # to prevent 100% CPU usage

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

data.test<-data.gen(n=100000,beta0,beta,p=0.06)
data.train<-data.gen(n=1000,beta0,beta,p=0.06)


## Standard method
thres<-seq(0,1,length.out =500)
stand.rr<-standard.lg(data.train=data.train,data.test=data.test,prev=NULL,alpha=alpha,KK=2,thres,
                      data.type="cohort")

beta.est.stand<-stand.rr$beta/norm(as.matrix(stand.rr$beta),"F")
standard.est<-stand.rr$rr
cutoff.rr<-stand.rr$cutoff

##DOOLR approach

beta.int<-stand.rr$beta/norm(as.matrix(stand.rr$beta),"F")
#hval = sd(as.matrix(cbind(1,X.train)) %*% matrix(beta.int,ncol=1))/(nrow(data.train)^(1/2))
#hval = sd(as.matrix(cbind(1,X.train)) %*% matrix(beta.null,ncol=1))/(nrow(data.train)^(0.5))
#h<-hval
ncores=7

doba.rr.new<-doolr.rule(beta.int=beta.int,data.train=data.train,
                           data.val=data.val,
                           data.test=data.test,prev=NULL,h=h,alpha=alpha
                           ,data.type="cohort")


## Plug-in approach

rr<-plugin.rule(data.train=data.train,data.val,data.test,prev=NULL,alpha=alpha,
                data.type="cohort",CV=FALSE,KK=5)




