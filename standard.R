



#'  Functions for constructing biomarker-based decision rule using standard approach; logistic regression
#'  
#' @param data.train a data frame of the training set; the first column is the outcome and the rest are the covariates .
#' @param data.test a data frame of the test set; the first column is the outcome and the rest are the covariates. 
#' @param prev Prevalence of the disease; defult is null if no prevalence value is provided otherwise estimate prevalence from data
#' @param alpha Pre-specified positive value (PPV) constraint
#' @param KK Number of cross-validation (CV); we use K-fold CV to choose the cutoff value
#' @param thres a vector of values bewteen 0 and 1 to choose a cutoff value for the rule
#' @param data.type a case-control or cohort
#' @return A list containing necessary information to compute p-value of the test
#' \item{standard.results}{Estimate of TPR, PPV, FPR}
#' \item{beta.est}{Estimates of the beta coefficients in the logistic regression model}
#' \item {cutoff.opt} {cutoff value the guarantee the pre-specified alpha value.}
#' 
#' @export
#' 

standard.lg<-function(data.train,data.test,prev=NULL,alpha,KK,thres,data.type){
  
  if(is.null(prev)){
    prev<-mean(data.train[,1]==1)
    
  }
  
  p<-ncol(data.train) # number of parameters
  rr.kfold<-kkfold.cutoff(prev,data=data.train,KK=KK,alpha,thres,data.type)
  
  cutoff.opt<-mode.fxn(rr.kfold$cutoffS)
  beta.est<-rr.kfold$beta
  
  n<-nrow(data.train)
  D.train<-data.train[,1]
  X.train<-subset(data.train,select=2:p)
  X.train<-data.matrix(X.train)
  
  beta.est<- as.matrix(glm(D.train ~ X.train, family ="quasibinomial")$coefficients)
  
  if(data.type=="case-control"){
    beta.est[1]<-beta.est[1] + log(prev/(1-prev)) - log(mean(data.train[,1]==1)/(1-mean(data.train[,1]==1)))
  }
  
  standard.results<-unlist(lg.est(prev,data=data.test,beta=beta.est, lambda=1,cutoff= cutoff.opt) )[1:3]
  
  return(list("rr"=standard.results,"beta"=beta.est,"cutoff"=cutoff.opt))
}


kkfold.cutoff<-function(prev,data,KK=5,alpha,thres,data.type) {
  
  if(is.null(prev)){
    prev<-mean(data[,1]==1)
  }
  p<-ncol(data) # number of parameters
  
  if(KK==1){
    
    n<-nrow(data)
    D.train<-data[,1]
    X.train<-subset(data,select=2:p)
    X.train<-data.matrix(X.train)
    beta.est<- as.matrix(glm(D.train ~  X.train, family ="quasibinomial")$coefficients)
    
    if(data.type=="case-control"){
      beta.est[1]<-beta.est[1] + log(prev/(1-prev)) - log(mean(data[,1]==1)/(1-mean(data[,1]==1)))
    }
    
    
    PPV  <- rep(NA, length(thres))
    TPF <- rep(NA, length(thres))
    FPF <- rep(NA, length(thres))
    
    for ( i in seq_along(thres) ) {
      
      FPF[i]<-lg.est(prev,data=data, beta= beta.est, lambda=1,cutoff=thres[i])$FPF
      TPF[i]<-lg.est(prev,data=data, beta= beta.est, lambda=1,cutoff=thres[i])$TPF
      PPV[i]<- prev*TPF[i]/(prev*TPF[i]+ (1-prev)*FPF[i])
    }
    
    temp              <- PPV - alpha
    temp1<-temp
    
    if(sum(temp1<0,na.rm=T)==length(thres)){
      index.opt<-which(abs(temp1)==min(abs(temp1),na.rm=T))[1]
    }
    
    else {
      temp1[temp1<0]<-1
      index.opt<-which(abs(temp1)==min(abs(temp1),na.rm=T))[1]
    }
    
    cutoff.opt<-thres[index.opt]
    return(list("cutoffS"=cutoff.opt,"beta"=beta.est))
  }
  
  else{

    folds <- createFolds(data[,1], k = KK, list = TRUE, returnTrain = TRUE)
    Beta<-matrix(NA,nrow=KK,ncol=p)
    cutoffS<-numeric() 
    TPFsMax<-numeric()
    
    for (kk in 1:KK) {
      
      # Create training and test sets
      train_indices <- unlist(folds[-kk])
      test_indices <- folds[[kk]]
      
      data.train <- data[train_indices, ]
      data.val <- data[test_indices, ]
      
      if(sum(data.train[,1]==1)==0 | sum(data.val[,1]==1)==0 ) {
        next
      }
      
      ### creating covariate and outcome for logistic regression ###
      n.train<-nrow(data.train)
      D.train<-data.train[,1]
      X.train<-subset(data.train,select=2:p)
      X.train<-data.matrix(X.train)
      D.val<-data.val[,1]
      X.val<-subset(data.val,select=2:p)
      X.val<-data.matrix(X.val)
      
      
      beta.est<- as.matrix(glm(D.train ~  X.train, family ="quasibinomial")$coefficients)

      if(data.type=="case-control"){
        beta.est[1]<-beta.est[1] + log(prev/(1-prev)) - log(mean(data.train[,1]==1)/(1-mean(data.train[,1]==1)))
      }
      
      
      PPV  <- rep(NA, length(thres))
      TPF <- rep(NA, length(thres))
      FPF <- rep(NA, length(thres))
      
      for ( i in seq_along(thres) ) {
        
        FPF[i]<-lg.est(prev,data=data.val, beta= beta.est, lambda=1,cutoff=thres[i])$FPF
        TPF[i]<-lg.est(prev,data=data.val, beta= beta.est, lambda=1,cutoff=thres[i])$TPF
        PPV[i]<- prev*TPF[i]/(prev*TPF[i]+ (1-prev)*FPF[i])
        
      }
      
      temp              <- PPV - alpha
      temp1<-temp
      threS<-thres[!is.na(temp1)]
      temp1<-temp1[!is.na(temp1)]
      if(sum(temp1<0,na.rm=T)==length(temp1)){
        index.opt<-which(abs(temp1)==min(abs(temp1),na.rm=T))[1]
      }
      
      else {
        temp1[temp1<0]<-1
        index.opt<-which(abs(temp1)==min(abs(temp1),na.rm=T))[1]
      }
      
      cutoff.opt<-threS[index.opt]
      TPFscut<-TPF[index.opt]
      
      if(length(cutoff.opt)>1)
      {
        cutoff.opt=cutoff.opt[which(TPFscut==max(TPFscut,na.rm=T))]
        TPFscut<-TPFscut[which(TPFscut==max(TPFscut,na.rm=T))]
      }
      
      cutoffS[kk]<-max(cutoff.opt,na.rm=T)
      Beta[kk,]<-beta.est
    }
    
    return(list("cutoffS"=mean(cutoffS),"beta"=colMeans(Beta)))
    
  }
}


mode.fxn <- function(x) {
  x<-na.omit(x)
  ux <- unique(x)
  ux[which.max(tabulate(match(x, ux)))]
}


lg.est<- function(prev,data,beta,lambda,cutoff) {
  
  if(is.null(prev)){
    prev<-mean(data[,1]==1)
  }
  
  p<-ncol(data) # number of parameters
  
  #### Testing step ######
  
  D.test<-data[,1]
  X.test<-subset(data,select=2:p)
  X.test<-data.matrix(X.test)
  
  f            <-  logistic(cbind(1,X.test)%*% beta)
  
  data$pred <- ifelse(f > cutoff, 1, 0)
  
  TPF<- mean((data[data[,1]==1,])$pred==1)
  FPF<- mean((data[data[,1]==0,])$pred==1)
  
  PPV<- prev*TPF/(prev*TPF+ (1-prev)*FPF)
  
  Diff         <- TPF - lambda*PPV
  
  return( list(TPF = TPF, PPV=PPV, FPF = FPF, cutoff = cutoff))
}


