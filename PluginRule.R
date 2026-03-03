
library(caret)
library(randomForest)
library(ranger)
library(SuperLearner)
library(parallel)

#'  Functions for constructing biomarker-based decision rule using the plug-in approach
#'  
#' @param data.train a data frame of the training set; the first column is the outcome and the rest are the covariates .
#' @param data.test a data frame of the test set; the first column is the outcome and the rest are the covariates. 
#' @param data.val a data frame of the validation set to estimate the constraint; the first column is the outcome and the rest are the covariates. 
#' @param prev Prevalence of the disease; if null, we estimate prevalence from data
#' @param alpha Pre-specified positive value (PPV) constraint
#' @param data.type a case-control or cohort type data
#' @return A data frame where each role corresponds to estimate of TPR, PPV, FPR for each plugin approach
#' \item{results.df}{Estimates of TPR, PPV, FPR on the test set for each plug-in approach}
#' @export
#' 


plugin.rule<-function(data.train,data.val,data.test,prev=NULL,alpha,data.type){
  
  lambda.hat<-NA
  
  if(is.null(prev)){
    prev<-mean(data.train[,1]==1)
  }
  
  ## training step 
  n.train<-nrow(data.train)
  D.train<-data.train[,1]
  n1<-sum(D.train==1)
  n0=sum(D.train==0)
  X.train<-subset(data.train,select=2:ncol(data.train))
  X.train<-data.matrix(X.train)
  p.hat<- prev
  gamma<-p.hat/(1-p.hat)
  
  fit.logistic<-glm(D.train~X.train,family="quasibinomial")
  theta <-fit.logistic$coefficients
  
  if(data.type=="case-control"){
    theta[1]<-theta[1] + log(prev/(1-prev)) - log(mean(data.train[,1]==1)/(1-mean(data.train[,1]==1)))
  }
  
  ## VALIDATION step ###
  D.val<-data.val[,1]
  X.val<-subset(data.val,select=2:ncol(data.train))
  X.val<-data.matrix(X.val)
  prob.est.val <-logistic(cbind(1,X.val)%*%theta)
  

    skip_to_next <- FALSE
    tryCatch(  
      lambda.hat<- uniroot(constraint,c(1e-4, 2000),prob=prob.est.val,
                           prev=prev,data=data.val)$root
      , error=function(e){skip_to_next<<-TRUE}) 
    
    if(is.na(lambda.hat)){
      
      skip_to_next <- FALSE
      tryCatch(  
        lambda.hat<-lambda_constA(f=constraint_alpha,alpha,a=1e-4,b=200,iter=1000,prob=prob.est.val,
                                  prev=prev,data=data.val)
        , error=function(e){skip_to_next<<-TRUE}) 
      constraint_alpha(lambda=100,prob=prob.est.val,
                       prev=prev,data=data.val)
    }
  
  skip_to_next <- FALSE
  
  est.lg<-pluginEst(prev,data.train=data.train,
                    data.test=data.test,alpha=alpha,lambda=lambda.hat,threshold=0,
                    method="logistic",data.type=data.type)
  
  
  ########### random forest model ########
  
  wt<-rep(c((1-prev)*(n1/(n0)),prev),times=c(n0,n1))
  
  rg.iris <- ranger( y = factor(D.train),x= X.train,  probability=TRUE)
  pred.iris <- predict(rg.iris, data = X.val)
  prob.est.val<-pred.iris$predictions[,2]
  
  if(data.type=="case-control"){
    
    prob.est1<-prob.est.val
    prob.est1<-ifelse(prob.est1==1,prob.est1-1e-14,prob.est1)
    prob.est1<-ifelse(prob.est1==0,1e-14,prob.est1)
    
    aa<-(prob.est1/(1-prob.est1))*(n0/n1)*(prev/(1-prev))
    prob.est.val<-aa/(1+aa)
    
  }
  
    skip_to_next <- FALSE
    tryCatch(  
      lambda.hat<- uniroot(constraint,c(1e-4, 2000),prob=prob.est.val,
                           prev=prev,data=data.val)$root
      , error=function(e){skip_to_next<<-TRUE}) 
    
    if(is.na(lambda.hat)){
      skip_to_next <- FALSE
      tryCatch(  
        lambda.hat<-lambda_constA(f=constraint_alpha,alpha,a=1e-4,b=200,iter=1000,prob=prob.est.val,
                                  prev=prev,data=data.val)
        , error=function(e){skip_to_next<<-TRUE}) 
    }
  
  ## evaluate on testing data
  est.rf<-pluginEst(prev,data.train=data.train,data.test=data.test,
                    alpha=alpha,lambda=lambda.hat,threshold=0,method="randomforest",data.type=data.type)
  
  ########### super learner ########
  
  SL.lib <- c("SL.mean", "SL.glm", "SL.gam")
  
  fit1 <- SuperLearner(Y = D.train, X = data.frame(X.train), 
                       SL.library = SL.lib, family = 'binomial', method = 'method.NNLS')
  
  prob.est.val<-predict(fit1, newdata = data.frame(X.val),onlySL = TRUE)$pred
  
  if(data.type=="case-control"){
    
    fit1 <- SuperLearner(Y = D.train, X = data.frame(X.train), 
                         SL.library = SL.lib, family = 'binomial', method = 'method.NNLS')
    
    prob.est1<-predict(fit1, newdata = data.frame(X.val),onlySL = TRUE)$pred
    aa<-(prob.est1/(1-prob.est1))*(n0/n1)*(prev/(1-prev))
    prob.est.val<-aa/(1+aa)
  }
  
  gamma<-p.hat/(1-p.hat)
  
    skip_to_next <- FALSE
    tryCatch(  
      lambda.hat<- uniroot(constraint,c(1e-4, 2000),prob=prob.est.val,
                           prev=prev,data=data.val)$root
      , error=function(e){skip_to_next<<-TRUE}) 
    
    if(is.na(lambda.hat)){
      skip_to_next <- FALSE
      tryCatch(  
        lambda.hat<-lambda_constA(f=constraint_alpha,alpha,a=1e-4,b=200,iter=1000,prob=prob.est.val,
                                  prev=prev,data=data.val)
        , error=function(e){skip_to_next<<-TRUE}) 
      
    }
  ## evaluate on testing data
  est.sl<-pluginEst(prev,data.train=data.train,data.test=data.test,
                    alpha=alpha,lambda=lambda.hat,threshold=0,method="superlearner",data.type=data.type)
  
  methods_names <- c("Plugin_Logistic", "Plugin_rf", "Plugin_SL")
  
  results.df <- data.frame(
    Methods = methods_names,
    TPF = c(est.lg[1], est.rf[1], est.sl[1]),
    PPV = c(est.lg[2], est.rf[2], est.sl[2]),
    FPF = c(est.lg[3], est.rf[3], est.sl[3]))
  return(results.df)
}



roots.lambda<- function(f,a,b,eps=1e-03,alpha,iter=1000){
  x=(a+b)/2
  x_new=x
  t=0
  while(t<=iter){
    
    if (f(a)*f(x)<=0) {
      b=x
      x_new=(a+b)/2
    }
    else{
      a=x
      x_new=(a+b)/2
    }
    if (abs(x_new-x)<eps|t>iter) {
      x=x_new
      break
    }
    x=x_new
    t=t+1
  }
  return (list("value"=x_new, "min"=f(x_new), "iter"= t))
}

##### solving for lambda that satisfies the constraint using bisection search###


lambda_constA<-function(f,alpha,a,b,iter=1000,prob,prev,data) {
  
  lambda_low   <-a
  lambda_up  <-b
  
  delta         <- 1
  ii            <- 1
  ALPHA         <- c()
  LAMBDA        <- c()
  PPV           <-numeric()
  
  lambdas<-seq(a,b,length.out=iter)
  
  for ( i in seq_along(lambdas) ) {
    
    PPV[i]<-f(lambdas[i],prob,prev,data)
  }
  temp              <- PPV - alpha
  temp1<-temp
  
  if(sum(temp1<0,na.rm=T)==length(lambdas)){
    index.opt<-which(abs(temp1)==min(abs(temp1),na.rm=T))[1]
    lambda  <- lambdas[index.opt]
  }
  
  else{ 
    temp1[temp1<0]<-1
    index.opts<-which(abs(temp1)==min(abs(temp1),na.rm=T))[1]
    lambda  <- lambdas[index.opts]
  }
  
  return(lambda)
  
}


lambda_const<-function(alpha,a,b,iter=200,data.train, data.val,method="logistic",p.hat=p.hat) {
  
  lambda_low   <- a
  lambda_up  <- b
  
  delta         <- 1
  ii            <- 1
  ALPHA         <- c()
  LAMBDA        <- c()
  PPV           <-numeric()
  
  
  while( delta >= 1e-5 && ii <= iter){
    
    lambda_temp  <- (lambda_low + lambda_up)/2
    LAMBDA       <- cbind(LAMBDA, lambda_temp)
    
    lambda=lambda_temp
    
    D.train<-data.train[,1]
    X.train<-subset(data.train,select=2:ncol(data.train))
    X.train<-data.matrix(X.train)
    D.val<-data.val[,1]
    X.val<-subset(data.val,select=2:ncol(data.train))
    X.val<-data.matrix(X.val)
    #p.hat<- mean(D.val==1)
    
    if(method=="logistic"){
      fit.logistic<-glm(D.train~X.train,family="quasibinomial")
      theta <-fit.logistic$coefficients
      prob.est.val <-logistic(cbind(1,X.val)%*%theta)
      
    }
    else if(method=="randomforest"){
      randomForest.fit<- randomForest(x = X.train,
                                      y = factor(D.train),
                                      xtest  = X.val,
                                      ytest  = factor(D.val),
                                      keep.forest=TRUE,
                                      importance = TRUE,mtry=ncol(X.train))
      
      prob1.est<- predict(randomForest.fit,X.val,"prob")[, 2]
      prob.est.val<-prob1.est
    }
    
    gamma<-p.hat/(1-p.hat)
    eta1<-prob.est.val/p.hat
    eta0<-(1-prob.est.val)/(1-p.hat)
    d.x<-  eta1 - gamma*alpha*lambda*eta1 - lambda*alpha*eta0 + lambda*gamma*eta1
   
    TPF<- mean((ifelse(d.x > threshold, 1, 0)*eta1))
    FPF<-mean((ifelse(d.x > threshold, 1, 0)*eta0))
    ppv<-p.hat*TPF/(p.hat*TPF + (1-p.hat)*FPF)
    
    if(is.na(ppv)){
      lambda=NA
      return(lambda)
    }
    alpha_temp   <- ppv
    
    ALPHA        <- cbind(ALPHA, alpha_temp)
    
    if (alpha_temp < alpha) {
      lambda_low   <- lambda_temp  
      
    } else {
      lambda_up   <- lambda_temp
    }
    delta        <- max(abs(lambda_up - lambda_low) )
    #delta        <- abs(f(lambda_temp)-alpha)
    ii           <- ii + 1  
  }
  
  
  if (sum(ALPHA > alpha) == length(ALPHA)) {
    lambda  <- lambda_low
    
  }
  else if (sum(ALPHA < alpha) == length(ALPHA)){
    lambda  <- lambda_up
  }
  else {
    lambda  <- min(LAMBDA[ALPHA >= alpha],na.rm=T)
  }
  
  return(lambda[1])
  
}


constraint<-function(lambda,prob,prev,data) {
  
  ### testing ##
  p.hat<- prev
  gamma<-p.hat/(1-p.hat)
  
  eta1<-prob/p.hat
  eta0<-(1-prob)/(1-p.hat)
  d.x<- eta1 - gamma*alpha*lambda*eta1 - lambda*alpha*eta0 + lambda*gamma*eta1
  data$pred <- ifelse(d.x > 0, 1, 0)
  TPF<- mean((data[data[,1]==1,])$pred==1)
  FPF<- mean((data[data[,1]==0,])$pred==1)
  PPV<-prev*TPF/(prev*TPF+ (1-prev)*FPF)
  ppvalpha<- gamma*TPF/(gamma*TPF +FPF)
  
  return(alpha-ppvalpha)

}


#### Function that returns estimate of PPV 

constraint_alpha<-function(lambda,prob,prev,data) {
 
  ### testing ##
  p.hat<- prev
  gamma<-p.hat/(1-p.hat)
  
  eta1<-prob/p.hat
  eta0<-(1-prob)/(1-p.hat)
  d.x<-  eta1 - gamma*alpha*lambda*eta1 - lambda*alpha*eta0 + lambda*gamma*eta1
  data$pred <- ifelse(d.x > 0, 1, 0)
  TPF<- mean((data[data$D==1,])$pred==1)
  FPF<- mean((data[data$D==0,])$pred==1)
  PPV<-prev*TPF/(prev*TPF+ (1-prev)*FPF)
  return(PPV)
}


## Function that returns estimate of (TPF,PPV,FPF) for each of the plug-in estimators

pluginEst<- function(prev,data.train,data.test,alpha,lambda,threshold,method,data.type) {
  
  prob.est<-numeric(nrow(data.test))
  
    if(method=="logistic"){
      
      D.train<-data.train[,1]
      X.train<-subset(data.train,select=2:ncol(data.train))
      X.train<-data.matrix(X.train)
      
      fit.logistic<-glm(D.train~X.train,family="quasibinomial")
      theta <-fit.logistic$coefficients
      
      if(data.type=="case-control"){
        theta[1]<-theta[1] + log(prev/(1-prev)) - log(mean(data.train[,1]==1)/(1-mean(data.train[,1]==1)))
      }
      
      ### Testing step ###
    
      D.test<-data.test[,1]
      X.test<-subset(data.test,select=2:ncol(data.train))
      X.test<-data.matrix(X.test)
      p.hat<- prev
      gamma<-p.hat/(1-p.hat)
      prob.est <- logistic(cbind(1,X.test)%*%theta)
      
    }
    
    else if (method == "randomforest") {
      n.train<-nrow(data.train) 
      D.train<-data.train[,1]
      n1<-sum(D.train==1)
      n0<-sum(D.train==0)
      X.train<-subset(data.train,select=2:ncol(data.train))
      X.train<-data.matrix(X.train)
      
      #### Testing step ######
      p.hat<- prev
      D.test<-data.test[,1]
      X.test<-subset(data.test,select=2:ncol(data.train))
      X.test<-data.matrix(X.test)
      
      rg.iris <- ranger( y = factor(D.train),x= X.train,  probability=TRUE)
      pred.iris <- predict(rg.iris, data = X.test)
      prob.est<-pred.iris$predictions[,2]
      
      
      if(data.type=="case-control"){
    
        prob.est1<-prob.est
        prob.est1<-ifelse(prob.est1==1,prob.est1-1e-14,prob.est1)
        prob.est1<-ifelse(prob.est1==0,1e-14,prob.est1)
        aa<-(prob.est1/(1-prob.est1))*(n0/n1)*(prev/(1-prev))
        prob.est<-aa/(1+aa)
        
      }
      
      gamma<-p.hat/(1-p.hat)
    }
    
    
    ### superlearner 
    
    else if (method == "superlearner") {
      
      n.train<-nrow(data.train) 
      D.train<-data.train[,1]
      n1<-sum(D.train==1)
      n0<-sum(D.train==0)
      X.train<-subset(data.train,select=2:ncol(data.train))
      X.train<-data.matrix(X.train)
      
      #### Testing step ######
      p.hat<- prev
      D.test<-data.test[,1]
      X.test<-subset(data.test,select=2:ncol(data.train))
      X.test<-data.matrix(X.test)
      
      SL.lib <- c("SL.mean", "SL.glm", "SL.gam")
      
      fit1 <- SuperLearner(Y = D.train, X = data.frame(X.train), 
                           SL.library = SL.lib, family = 'binomial', method = 'method.NNLS')
      
      prob.est<-predict(fit1, newdata = data.frame(X.test),onlySL = TRUE)$pred
      
      if(data.type=="case-control"){
        
        fit1 <- SuperLearner(Y = D.train, X = data.frame(X.train), 
                             SL.library = SL.lib, family = 'binomial', method = 'method.NNLS')
        
        prob.est1<-predict(fit1, newdata = data.frame(X.test),onlySL = TRUE)$pred
        
        aa<-(prob.est1/(1-prob.est1))*(n0/n1)*(prev/(1-prev))
        prob.est<-aa/(1+aa)
        
      }
      
      gamma<-p.hat/(1-p.hat)
    }
    
  #### estimating decision rule
  p.hat<- prev
  gamma<-p.hat/(1-p.hat)
  eta1<-prob.est/p.hat
  eta0<-(1-prob.est)/(1-p.hat)
  d.x<-  eta1 - gamma*alpha*lambda*eta1 - lambda*alpha*eta0 + lambda*gamma*eta1  
  
  data.test$pred <- ifelse(d.x > 0, 1, 0)
  
  TPF<- mean((data.test[data.test[,1]==1,])$pred==1)
  FPF<- mean((data.test[data.test[,1]==0,])$pred==1)
  PPV<- prev*TPF/(prev*TPF + (1-prev)*FPF)
  ## estimators ###
  if(TPF==1 & PPV==1){
    print("no active constraint") 
    return(c(NA,NA,NA))
  }
  else {
  return(c(TPF,PPV,FPF))
  }
  
}


