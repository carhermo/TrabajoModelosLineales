## WRAPPER FOR ALL THE WORK
# author(mostly): Carlos Hernani Morales, 2019, carhermo@protonmail.com
library(leaps)
library(tidyverse)


predict.regsubsets = function(object,newdata,id,...){
  form = as.formula(PTS~X2PA + X3PA + FTA + AST + ORB + DRB + TOV + STL + BLK) 
  mat = model.matrix(form,newdata)   
  coefi = coef(object,id=id)         
  xvars = names(coefi)               
  mat[,xvars]%*%coefi               
}

TSS <- function(y){
  tss <- sum((y-mean(y))^2)
  return(tss)
}


RSS <- function(y,ypred){
  rss <- sum((y-ypred)^2)
  return(rss)
}


R2 <- function(y,ypred){
  r2 <- 1 - (RSS(y,ypred)/TSS(y))
  return(r2)
}


RMSE <- function(y,ypred){
  rmse <- sqrt(RSS(y,ypred)/(length(y)-2))
  return(rmse)
}


gghighlight <- function(dataframe,what="Cp",searchinto="Value",searchfor="min"){
  if(!(searchfor %in% c("min","max"))){
    stop("searchfor is min or max.")
  }
  else{
    if(searchfor=="min"){
        min <- dataframe %>% filter(Type == what) %>% filter(Value==min(Value))
        plot <- geom_point(data=min, aes(x=dims,y=Value), color='red',size=3)
    }
    else{
      max <- dataframe %>% filter(Type == what) %>% filter(Value==max(Value))
      plot <- geom_point(data=max, aes(x=dims,y=Value), color='red',size=3)
    }
  }
  return(plot)
}


tidymodel <- function(summarylist,criteria,testrmse,nvmax){
  if(criteria=="indirect"){
    Cp <- summarylist$cp
    bic <- summarylist$bic
    adjr2 <- summarylist$adjr2
    dims <- rep(seq(1,length(summarylist$rsq)),3)
    Type <- c(rep("Cp",length(summarylist$rsq)),rep("BIC",length(summarylist$rsq)),rep("R2-adj",length(summarylist$rsq)))
    Value <- c(Cp,bic,adjr2)
    dfmod <- data.frame(dims,Type,Value)
  }
  if(criteria == "validation"){ 
    r2 <- summarylist$rsq
    rss <- summarylist$rss
    n <- summarylist$obj$d[1]
    rmse.train <- sqrt(rss/(n-2))
    rmse.test <- testrmse
    dims <- rep(seq(1,length(summarylist$rsq)),2)
    Type <- c(rep("RMSE train",length(summarylist$rsq)),rep("RMSE test",length(summarylist$rsq)))
    Value <- c(rmse.train,rmse.test)
    dfmod <- data.frame(dims,Type,Value)
  }
  if(criteria=="crossval"){
    dims <- seq(1,9)
    Type <- rep("RMSE CV",9)
    Value <- unlist(testrmse)
    dfmod <- data.frame(dims,Type,Value)
    
  }
  return(dfmod)
}


model_selection <- function(vresp.test,vresp.train,trainset,testset,nvarmax=6,kfold=10,seed=1,method="exhaustive",crit="indirect"){
  #check valid method
  if(!(method %in% c("exhaustive", "backward", "forward", "seqrep"))){
    stop("Invalid Method: Valid methods:exhaustive, backward, forward, seqrep")
  }
  #check valid criteria
  if(!(crit %in% c("indirect", "validation", "crossval"))){
    stop("Invalid Criteria: Valid criteria:indirect, validation, crossval")
  }
  
  model <- regsubsets(PTS~X2PA + X3PA + FTA + AST + ORB + DRB + TOV + STL + BLK,data=trainset,nvmax=nvarmax,method = method)
  model.summary <- summary(model)
  
  
  if(crit == "indirect"){
    #model selection by Cp,Bic,R2adjusted
    #we tidy the info
    df.model <- tidymodel(model.summary,criteria = crit )
    #lets get the plot
    baselayer <- df.model %>% ggplot(aes(x=dims,y=Value,color=Type)) + geom_line() + geom_point()
    highlighted.Cp <- gghighlight(df.model) #min Cp AIC
    highlighted.bic <- gghighlight(df.model,what = "BIC") #min BIC
    highlighted.r2adj <- gghighlight(df.model,what = "R2-adj",searchfor = "max")
    layers <- baselayer + highlighted.Cp + highlighted.bic + highlighted.r2adj
    
  }
  
  if(crit== "validation"){
    if(missing(testset)){
      stop("Validation criteria needs a validation set.")
    }
    #model selection by validation set
    val_errors = rep(NA,nvarmax)
    val_rmse = rep(NA,nvarmax)

    for(i in 1:nvarmax){
      pred = predict.regsubsets(model,testset,i)
      
      val_errors[i] = mean((vresp.test-pred)^2)
      val_rmse[i] = RMSE(vresp.test,predict.regsubsets(model,testset,i))
    }
    #we tidy the info
    df.model <- tidymodel(model.summary,criteria = crit,testrmse = unlist(val_rmse) )
    #lets get the plot
    baselayer <- df.model %>% ggplot(aes(x=dims,y=Value,color=Type)) + geom_line() + geom_point()
    highlighted.rmse.train <- gghighlight(df.model,what = "RMSE train")
    highlighted.rmse.test <- gghighlight(df.model,what = "RMSE test")
    layers <- baselayer + highlighted.rmse.train + highlighted.rmse.test
  }
  
  if(crit == "crossval"){
    #model selection by crossvalidation
    set.seed(1)   
    
    folds = sample(1:kfold, nrow(trainset), replace = TRUE)
    
    cv_errors = matrix(NA, kfold, nvarmax, dimnames = list(NULL, paste(1:nvarmax)))
    cv_rmse = matrix(NA, kfold, nvarmax, dimnames = list(NULL, paste(1:nvarmax)))
    for(j in 1:kfold){
      
      model <- regsubsets(PTS~X2PA + X3PA + FTA + AST + ORB + DRB + TOV + STL + BLK,trainset[folds!=j,],nvmax = nvarmax)
      
      for(i in 1:nvarmax){
        
        pred = predict(model, trainset[folds==j,], id=i)
        
        cv_errors[j,i] = mean((vresp.train[folds==j]-pred)^2)
        cv_rmse[j,i] = RMSE(vresp.train[folds==j],pred)
      }
    }
    
    mean_cv_errors = apply(cv_errors, 2, mean)
    mean_cv_rmse = apply(cv_rmse, 2, mean)
    df.model <- tidymodel(criteria = crit,testrmse = mean_cv_rmse,nvmax=nvarmax)
    #lets get the plot
    baselayer <- df.model %>% ggplot(aes(x=dims,y=Value,color=Type)) + geom_line() + geom_point()
    highlighted.rmse.cv <- gghighlight(df.model,what = "RMSE CV")
    layers <- baselayer + highlighted.rmse.cv
  }
  
  #returnees
  retu <- list(method,crit,formula,nvarmax,df.model,layers)
  names(retu) <- c("method","criteria","formula","nvmax","df.model","layers")
  return(retu)

  
}