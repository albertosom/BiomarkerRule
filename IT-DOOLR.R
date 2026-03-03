
library(parallel)


#'  Functions for constructing biomarker-based decision rule using the DOOLR with information transfer
#'  
#'@param beta.int a vector of initial parameter values for the linear rule
#'@param data.train a data frame of the training set; the first column is the outcome and the rest are the covariates .
#'@param data.test a data frame of the test set; the first column is the outcome and the rest are the covariates. 
#'@param  data.val a data frame of the validation set; the first column is the outcome and the rest are the covariates. 
#'@param prev Prevalence of the disease; if NULL, we estimate prevalence from data
#'@param alpha Pre-specified positive value (PPV) constraint
#'@param h scaling parameter to choose for approximating indicator function
#'@param risk_yes_no individual risk decisions from an external risk model (Yes: high risk, No: low risk)
#'@param eta tuning parameter for information transfer
#'@return A list containing estimate of parameters of the rule and estimates of TPR, PPV, and FPR
#' \item{beta}{Estimate of the parameters of the linear rule}
#' \item{measures}{ Estimate of TPR, PPV, FPR}
#' @export

doolr.rule.IT<-function(beta.int,data.train,data.val,data.test,prev,alpha,h,data.type,risk_yes_no,eta){
  
  if(is.null(prev)){
    prev<-mean(data.train[,1]==1)
  }
  
  p<-ncol(data.train) #number of parameters
  
  lambdas<-seq(1e-4,100,length.out=200)
  
  Betas<-matrix(NA,ncol=p,nrow=length(lambdas))
  ppvs<-numeric()
  
  mc.rr<-mclapply(1:length(lambdas), function(i){
    
    # for(i in 1:length(lambdas)){
    
    lambda<-lambdas[i]
    
    skip_to_next <- FALSE
    tryCatch(  
      #tMax.optim <- optim(par = beta.int,fn=object.func.doba.KL,data.train=data.train,
      # control=list(maxit=500),
      #                data.val=data.val,h=h,alpha=alpha,prev=prev,data.type=data.type,
      #               riskscore=riskscore,eta=eta,
      #              lambda.hat=lambda,method="BFGS")
      
      nlmm<-nlm(object.func.doolr.IT,beta.int,data.train = data.train,data.val=data.val,
                h=h,alpha=alpha,prev=prev,data.type=data.type,risk_yes_no=risk_yes_no,eta=eta,
                lambda.hat=lambda,steptol=1e-4,gradtol=1e-4)
      , error=function(e){ 
        skip_to_next<<-TRUE}) 
    
    if(skip_to_next) { 
      beta.opt<-matrix(rep(NA,ncol(data.train)),ncol=1)
    }  
    
    else {
      beta.opt<-matrix(nlmm$estimate,ncol=1)
      #beta.opt<-matrix(tMax.optim$par,ncol=1)
    }
    beta.opt<-beta.opt/norm(beta.opt,"F")
    soln<-unlist(doolr.est.IT(prev=prev,data.test=data.val,
                                 beta=beta.opt,h=h,alpha=alpha,data.type=data.type))
    ppvs<-soln[2]
    TPFs<-soln[1]
    # }
    return(cbind(t(beta.opt),ppvs,TPFs))
  } ,mc.cores = ncores)
  
  mc.rrr<-matrix(unlist(mc.rr), ncol=p+2, byrow=TRUE)
  Betas<-mc.rrr[,1:(ncol(data.train))]
  ppvs<- mc.rrr[,(ncol(data.train))+1]
  TPFs<- mc.rrr[,(ncol(data.train))+2]
  temp <- ppvs - alpha
  temp1<-temp
  
  if(sum(is.na(ppvs))==length(lambdas)){
    beta.opt<-matrix(rep(NA,p),ncol=1)
    cat(paste0("no feasible solution for eta = ",eta) ,"\n")
    return(list("beta"=beta.opt,"measures"=c(NA,NA,NA)))
  }
  
  if(sum(temp1<0,na.rm=T)==length(lambdas)){
    index.opts<-which(abs(temp1)==min(abs(temp1),na.rm=T))
    index.opt<-index.opts
  } else {
    temp1[temp1<0]<-1
    index.opts<-which(abs(temp1)==min(abs(temp1),na.rm=T))
    index.opt<-index.opts
  }
  
  temp1  <- temp1[index.opts]
  TPFscut<-TPFs[index.opts]
  
  if(length(index.opts)>1){
    #cutoff.opt=max(cutoff.opt)
    index.opt<-which(TPFscut==max(TPFscut,na.rm=T))[1]
    index.opt<-index.opts[index.opt]
    lambda=lambdas[index.opt]
    TPFscut<-TPFs[index.opt]
    ##ppvs[index.opt]
  }
  
  if(TPFscut==1 & ppvs[index.opt]==prev ){
    
    beta.opt<-matrix(rep(NA,p),ncol=1)
    cat("no feasible solution: change initial value","\n") 
    return(list("beta"=beta.opt,"measures"=c(NA,NA,NA)))
  }
  else{
    beta.opt<-matrix(Betas[index.opt,],ncol=1)
    beta.opt<-beta.opt/norm(beta.opt, type="F")
    aa<-unlist(doolr.est.IT(prev=prev,data.test=data.test,
                               beta=beta.opt,h=h,alpha=alpha,data.type=data.type))[1:3]
    return(list("beta"=beta.opt,"measures"=aa))
  }
  
}


doolr.rule.IT.CV<-function(beta.int,data.train,data.val,data.test,prev,alpha,h,data.type,risk_yes_no,etas){
  
  if(is.null(prev)){
    prev<-mean(data.train[,1]==1)
  }
  
  TPFs<-numeric()
  PPVs<-numeric()
  Betas<-matrix(NA,nrow=length(etas),ncol=ncol(data.train))
  
  if( length(etas)>1){
    
    for(i in seq_along(etas)){
      eta<-etas[i]
      aa<-doolr.rule.IT(beta.int=beta.int,data.train=data.train,
                       data.val=data.val,
                       data.test=data.train,prev=prev,alpha=alpha,h=h,
                       data.type=data.type,risk_yes_no=risk_yes_no,eta=eta)
      
      TPFs[i]<- aa$measure[1]
      PPVs[i]<- aa$measure[2]
      Betas[i,]<-aa$beta
    }
    
    index.less<-PPVs-alpha<0
    temp<-PPVs-alpha
    
    if (sum(index.less, na.rm = TRUE) == length(etas)) {
      # All PPVs < alpha: choose index closest to alpha
      index.opts <- which(abs(PPVs - alpha) == min(abs(PPVs - alpha), na.rm = TRUE))[1]
    } else {
      # Some PPVs >|< alpha: filter for those, and take max TPF
      PPVs.filtered <- PPVs[!index.less]
      TPFs.filtered <- TPFs[!index.less]
      temp.filtered <- temp[!index.less]
      
      index.opts <- which(TPFs.filtered == max(TPFs.filtered, na.rm = TRUE))
    }
    
    if (length(index.opts) > 1) {
      tempcut <- temp.filtered[index.opts]
      index.sub <- which(abs(tempcut) == min(abs(tempcut), na.rm = TRUE))
      index.opts <- index.opts[index.sub]
      # Map back to original index
      index.opts <- which(!index.less)[index.opts]
    }
    
    
    # Final choice (first if multiple)
    index.opt <- index.opts[1]
    
    beta.opt<-matrix(Betas[index.opt,],ncol=1)
    aa<-unlist(doolr.est.IT(prev=prev,data.test=data.test,
                               beta=beta.opt,h=h,alpha=alpha,data.type=data.type))[1:3]
    
    return(list("beta"=beta.opt,"measures"=aa))
  }
  
  else {
    
    rr<- doolr.rule.IT(beta.int=beta.int,data.train=data.train,
                      data.val=data.val,
                      data.test=data.test,prev=prev,alpha=alpha,h=h,
                      data.type=data.type,risk_yes_no=risk_yes_no,eta=etas)
    return(rr)
  }
}



doolr.est.IT<-function(prev,data.test,beta,h,alpha,data.type){
  
  if(is.null(prev)){
    prev<-mean(data.test[,1]==1)
  }
  
  p<-ncol(data.test)
  #beta<-rbind(-1,beta)
  
  #### Testing data ######
  D.test<-data.test[,1]
  X.test<-subset(data.test,select=2:p)
  X.test<-data.matrix(X.test)
  
  #f            <-  cbind(1,X.test)%*% beta
  f            <-  cbind(1,X.test)%*% beta
  
  ### difference approlambdamating methods
  
  p.hat<-  prev
  gamma<-p.hat/(1-p.hat)
  
  data.test$pred <- ifelse(f > 0, 1, 0)
  
  TPF<- mean((data.test[data.test[,1]==1,])$pred==1)
  FPF<- mean((data.test[data.test[,1]==0,])$pred==1)
  PPV<- prev*TPF/(prev*TPF+ (1-prev)*FPF)
  
  return( list(TPF=TPF,PPV=PPV, FPF = FPF) )
}


object.func.doolr.IT<-function(beta,data.train,data.val,prev,alpha,h,lambda.hat,data.type,risk_yes_no,eta){
  
  if(is.null(prev)){
    prev<-mean(data.train[,1]==1)
  }
  
  p<-ncol(data.train)
  n.train<-nrow(data.train)
  D.train<-data.train[,1]
  n1<-sum(D.train==1)
  n0<-sum(D.train==0)
  X.train<-subset(data.train,select=2:p)
  X.train<-data.matrix(X.train)
  
  X1 <- as.matrix(X.train[D.train==1,])
  X0 <- as.matrix(X.train[D.train==0,])
  
  f_1          <- cbind(1,X1)%*% beta
  f_0          <- cbind(1,X0)%*% beta
  
  X.beta<-cbind(1,X.train)%*% beta
  
  f.beta<-logistic(X.beta)
  #f.risk<-logistic(riskscore)
  risk_yes_no<-as.numeric(risk_yes_no)
  risk<-risk_yes_no[D.train==1]
  loss<- mean(pnorm(-f_1*risk/h))
  p.hat<- prev
  gamma<-p.hat/(1-p.hat)
  #eta1<-prob.est/p.hat
  #eta0<- (1-prob.est)/(1-p.hat)
  
  TPF<-mean(pnorm(f_1/h))
  FPF<-mean(pnorm(f_0/h))
  PPV<-gamma*TPF/(gamma*TPF+FPF)
  #obj<- -1*mean(pnorm(f_1/h)) + mean(pnorm(f_1/h))*lambda.hat*alpha*gamma + lambda.hat*alpha*(mean(pnorm(f_0/h))) - lambda.hat*gamma*mean(pnorm(f_1/h))
  
  if(eta==0){
    obj<- -1*TPF + lambda.hat*TPF*lambda.hat*alpha*gamma + lambda.hat*alpha*FPF - lambda.hat*gamma*TPF
  }
  else{
    obj<- -1*TPF + lambda.hat*TPF*lambda.hat*alpha*gamma + lambda.hat*alpha*FPF - lambda.hat*gamma*TPF + 
      eta*loss
  }
  return(obj)
  
}


const.func.doolr<-function(beta,data.train,data.val,prev,alpha,h,lambda.hat,data.type){
  
  if(is.null(prev)){
    prev<-mean(data.val[,1]==1)
  }
  
  p<-ncol(data.val)
  
  n.val<-nrow(data.val)
  D.val<-data.val[,1]
  n1<-sum(D.val==1)
  n0<-sum(D.val==0)
  X.val<-subset(data.val,select=2:p)
  X.val<-data.matrix(X.val)
  
  X1 <- as.matrix(X.val[D.val==1,])
  X0 <- as.matrix(X.val[D.val==0,])
  
  f_1          <- cbind(1,X1)%*% beta
  f_0          <- cbind(1,X0)%*% beta
  
  #f_1          <- cbind(X1)%*% beta
  #f_0          <- cbind(X0)%*% beta
  
  
  p.hat<- prev
  gamma<-p.hat/(1-p.hat)
  #eta1<-prob.est/p.hat
  #eta0<- (1-prob.est)/(1-p.hat)
  
  TPF<-mean(pnorm(f_1/h))
  FPF<-mean(pnorm(f_0/h))
  PPV<-gamma*TPF/(gamma*TPF+FPF)
  
  #const<- mean(pnorm(f_1/h)) + mean(pnorm(f_1/h))*alpha*gamma + alpha*(mean(pnorm(f_0/h))) - gamma*mean(pnorm(f_1/h))
  const<- TPF*alpha*gamma + alpha*FPF - gamma*TPF
  
}

constrnorm.atl <- function(beta,data.train,data.val,alpha,h,prev,lambda.hat,data.type){
  betavec <- beta
  norm(matrix(betavec), type="F") 
}
