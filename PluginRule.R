

#'  Functions for constructing biomarker-based decision rule using the plug-in approach
#'  
#' @param data.train a data frame of the training set; the first column is the outcome and the rest are the covariates .
#' @param data.test a data frame of the test set; the first column is the outcome and the rest are the covariates. 
#'@param  data.val a data frame of the validation set; the first column is the outcome and the rest are the covariates. 
#' @param prev Prevalence of the disease; defult is null if no prevalence value is provided otherwise estimate prevalence from data
#' @param alpha Pre-specified positive value (PPV) constraint
#' @param KK Number of cross-validation (CV); we use K-fold CV to choose the cutoff value
#' @param thres a vector of values bewteen 0 and 1 to choose a cutoff value for the rule
#' @param data.type a case-control or cohort
#' @param CV whether to use cross-validation (TRUE/FALSE)
#' @return A data frame where each role corresponds to estimate of TPR, PPV, FPR for each plugin approach
#' \item{results.df}{Estimates of TPR, PPV, FPR for eahc pluin approach}
#' @export
#' 


plugin.rule<-function(data.train,data.val,data.test,prev=NULL,alpha,data.type,CV=FALSE,KK=2){
  
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
  
  if(CV==TRUE){
    lambda.hat<-kkfold.lambda(data=data.train,KK=2,method="logistic",alpha,prev)
  }
  if(CV==FALSE){
    
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
  
  if(CV==TRUE){
    lambda.hat<-kkfold.lambda(data=data.train,KK=2,method="randomforest",alpha,prev)
  }
  
  if(CV==FALSE){
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
  
  if(CV==TRUE){
    lambda.hat<-kkfold.lambda(data=data.train,KK=2,method="superlearner",alpha,prev)
  }
  
  if(CV==FALSE){
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
  }
  
  ## evaluate on testing data
  est.sl<-pluginEst(prev,data.train=data.train,data.test=data.test,
                    alpha=alpha,lambda=lambda.hat,threshold=0,method="superlearner",data.type=data.type)
  
  results.df<-rbind(est.lg,est.rf,est.sl)
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


#### CONSTRAINT THAT evaluates to alpha

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



pluginEst<- function(prev,data.train,data.test,alpha,lambda,threshold,method,data.type) {
  
  prob.est<-numeric(nrow(data.test))
  
    if(method=="logistic") {
      
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
  return(c(TPF,PPV,FPF))
}

kkfold.lambda<-function(data,KK,method="logistic",alpha,prev) {
  
  folds <- createFolds(data[,1], k = KK, list = TRUE, returnTrain = TRUE)
  
  lambdas<-numeric() 
  
  p.hat<- prev
  
  PPVs  <- rep(NA, KK)
  TPFs <- rep(NA, KK)
  FPFs<-rep(NA, KK)
  
  for (kk in 1:KK) {
    
    # Create training and test sets
    
    train_indices <- unlist(folds[-kk])
    test_indices <- folds[[kk]]
    
    data.train <- data[train_indices, ]
    data.val <- data[test_indices, ]
    lambda.hat<-NA
    
    if( sum(data.val$D==1)==0 ) {
      next
    }
    
    ### creating covariate and outcome for logistic regression ###
    n<-nrow(data.train)
    D.train<-data.train[,1]
    X.train<-subset(data.train,select=2:ncol(data.train))
    X.train<-data.matrix(X.train)
    
    if(method=="logistic") {
      
      fit.logistic<-glm(D.train~X.train,family=binomial(link = "logit"))
      theta <-fit.logistic$coefficients
      ### Validation step ###
      D.val<-data.val[,1]
      X.val<-subset(data.val,select=2:ncol(data.train))
      X.val<-data.matrix(X.val)
      # p.hat<- mean(D.val==1)
      gamma<-p.hat/(1-p.hat)
      prob.est.val <-logistic(cbind(1,X.val)%*%theta)
      
      lambda.hat<- lambda_const(alpha,a=1e-4,b=200,iter=200,data.train, data.val,method="logistic",p.hat=p.hat)
    
    }
    
    else if (method == "randomforest") {
      
      ### Validation step  ###
      D.val<-data.val[,1]
      X.val<-subset(data.val,select=2:ncol(data.frame))
      X.val<-data.matrix(X.val)
      
      gamma<-p.hat/(1-p.hat)
      
      rg.iris <- ranger( y = factor(D.train),x= X.train,  probability=TRUE)
      pred.iris <- predict(rg.iris, data = X.val)
      prob.est.val<-pred.iris$predictions[,2]
      
      lambda.hat<- lambda_const(alpha,a=1e-4,b=200,iter=200,data.train, data.val,method="randomforest",p.hat=p.hat)
    }
    
    else if(method=="superlearner") {
      
      D.val<-data.val[,1]
      X.val<-subset(data.val,select=2:ncol(data.train))
      X.val<-data.matrix(X.val)
     
      gamma<-p.hat/(1-p.hat)
      
      SL.lib <- c("SL.mean", "SL.glm", "SL.gam")
      
      fit1 <- SuperLearner(Y = D.train, X = data.frame(X.train), 
                           SL.library = SL.lib, family = 'binomial', method = 'method.NNLS')
      
      prob.est.val<-predict(fit1, newdata = data.frame(X.val),onlySL = TRUE)$pred
      lambda.hat<- lambda_const(alpha,a=1e-4,b=100,iter=200,data.train, data.val,method="randomforest",p.hat=p.hat)
      
    }
    
    if(lambda.hat==Inf|is.na(lambda.hat)){
      next
    }
    lambdas[kk]<-lambda.hat  
    lambda<-lambda.hat
    
    eta1<-prob.est.val/p.hat
    eta0<-(1-prob.est.val)/(1-p.hat)
    d.x<-  eta1 - gamma*alpha*lambda*eta1 - lambda*alpha*eta0 + lambda*gamma*eta1
    data.val$pred <- ifelse(d.x > threshold, 1, 0)
    
    TPFs[kk]<- mean((data.val[data.val[,1]==1,])$pred==1)
    PPVs[kk]<- ifelse(sum(data.val$pred==1)==0,0,mean((data.val[data.val$pred==1,])[,1]==1))
    FPFs[kk]<- mean((data.val[data.val[,1]==0,])$pred==1)
    
  }
  PPVs<-na.omit(PPVs)
  TPFs<-na.omit(TPFs)
  
  temp<-PPVs-alpha
  
  temp1<-temp
  temp1[temp1 < 0] <- 1
  min.constraint<-which(temp1==min((temp1),na.rm=T))
  lambda.opt<-lambdas[min.constraint]
  TPFscut<-TPFs[min.constraint]
  
  if(length(lambda.opt)>1)
  {
    lambda.opt=lambda.opt[which(TPFscut==max(TPFscut,na.rm=T))]
  }
  
  return(lambda.opt[1])
  
}

