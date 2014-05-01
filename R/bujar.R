bujar <- function(y, cens, x, valdata = NULL, glm = TRUE, degree = 1, learner = "linear.regression", tran = "log", center=TRUE, mimpu = NULL, iter.bj = 10, max.cycle = 10, nu = 0.1, mstop = 50, tuning = TRUE, cv = FALSE, nfold = 5, method = "corrected", df = "actset", vimpint = TRUE, gamma=3, lambda=NULL, whichlambda, lamb = 0, s = 0.5, nk = 4, wt.pow = 1, theta = NULL, rel.inf = FALSE, tol = .Machine$double.eps, trace = FALSE){
  call <- match.call()
  if(learner == "acosso")
  stop("learner = 'acosso' is no longer supported, see NEWS\n")
  if(!learner%in%c("linear.regression","mars","pspline","tree","acosso","enet","MCP","SCAD")) stop(sQuote("weak learner"), " is not implemented")
  if(learner%in%c("mars","pspline","tree","acosso") && glm==TRUE) stop(sQuote("glm must be FALSE for the selected learner boosting"))
  if(learner%in%c("linear.regression","enet","MCP","SCAD") && glm==FALSE) stop(sQuote("glm must be TRUE for the selected linear.regression or enet, or MCP, or SCAD BJ method"))
  if(!is.null(valdata))
    if((dim(x)[2]+2) !=dim(valdata)[2])
      stop("check the dimension of x and valdata\n")
  if(learner=="acosso" && wt.pow < 0) stop(Quote("wt.pow should be > 0"))
  if(learner=="acosso" && cv) stop(Quote("if wt.pow is chosen by cross-validation, then BJ estimator is not stable, thus stop"))### check ?
  if(cv && nfold < 1)
    stop(sQuote("Number of CV folds"), " are less than 1")
  if(!all(unique(cens)%in%c(0,1)))
    stop(sQuote("censoring indicator"), " are not 0 or 1")
  if(!all(unique(valdata[,2]%in%c(0,1))))
    stop(sQuote("censoring indicator"), " are not 0 or 1")
  if(iter.bj < 2)
    stop(sQuote("iter.bj"), " should be greater than 1")
  if(max.cycle < 1)
    stop(sQuote("max.cycle"), " should be greater than 1")
  if(!learner %in% c("tree","mars","acosso") && degree !=1)
    stop("Not implemented for this degree with learner ",learner, "\n")
  xbar <- colMeans(x)                                    #componentwise linear least square
  f <- 0
  mse.bj.val <- nz.bj <- NA
  ynew <- y; ynew.m <- vim <- NULL; Fboost <- NA
  p <- ncol(x)
  sse <- rep(NA, p) #SSE for each univariate covariate
  res <- ystar <- matrix(NA, length(y), p) #residual matrix, one column for corresponding univariate covariate
  coef <- matrix(NA, 2, p)
  b <- matrix(NA,iter.bj+max.cycle,p)
  ybst <- matrix(NA, iter.bj+max.cycle, length(y)) #changed Apr 2, 2008
  ybstdiff <- rep(NA, iter.bj+max.cycle)       #changed Apr 2, 2008
  fnorm2 <- mseun <- ybstcon <- rep(NA, iter.bj+max.cycle)       #changed Apr 2, 2008
                                        #  res.model <- vector("list",iter.bj+max.cycle)
  mselect <- rep(NA,iter.bj+max.cycle)
  best.iter <- mstop
  tuningSwitch <- TRUE
  ydiff <- 100
  k <- 1; kt <- 1
  mse.bj <- pred.bj <- NA
  if(trace) cat("\nBJ with",learner,"\n")
  nm <- dim(x)
  n <- nm[1]
  if(trace){
    cat("\nNumber of observations:",n)
    cat("\nNumber of covariates:",nm[2],"\n")
  }
  if(learner=="linear.regression" && all(x[,1]==1)) ix <- 2:p
  else ix <- 1:p
  one <- rep(1, n)
  normx <- rep(1,dim(x)[2])
  cycleperiod <- 0
  nonconv <- FALSE
  mse <- rep(NA,iter.bj+max.cycle)
  nz.bj.iter <- rep(NA,iter.bj+max.cycle)
  mstop.mars <- mstop
  while (ydiff > tol && k <= iter.bj+max.cycle){
    oldydiff <- ydiff
                                        #when k=1, only for initializaiton purpose, find BJ estimator and predicted y-value                                      #determine ynew
### Initializtion of BJ estimator with different (3) methods
    if(is.null(mimpu) && k==1) ynew <- y #this is to say we begin with initial estimator of beta=0
    else if(mimpu==FALSE && k==1){
      for (i in 1:p){
        res.des <- try(bj(Surv(ynew,cens) ~ x[,i], link="identity",control=list(trace=FALSE)))
        ystar[,i] <- res.des$y.imputed #imputed y value
        res[,i] <- y - predict(res.des)
        sse[i] <- sum(res[,i][cens==1]^2) # SSE for uncensored
      }
      minid <- which.min(sse) #find the covariate with smallest loss
      ynew <- ystar[,minid] #construct the new 'response'
      cat("\nBJ step k=",k,"\n","\nInitial MSE for uncensored observations=", sse[minid]/sum(cens==1), "\n\n")
    }
    else{
      if(mimpu==TRUE && k==1){
        ynew <- bjboost.fit(cbind(y,cens),rep(0,length(y)))$y.imputed #imputing without covariates
                                        #        k <- iter.bj
      }
      else {
        ynew <- bjboost.fit(cbind(y,cens),Fboost)$y.imputed}
    }
### End of initialization with the imputed ynew replacing the original y values
    dat1 <- data.frame(cbind(ynew,x))
### Different methods for BJ estimation
    x <- as.matrix(x) 
    if(glm && learner=="linear.regression"){
      ctrl <- boost_control(nu=nu,mstop=mstop)
      dat1.glm <- glmboost(ynew~x,center=center, control=ctrl)
    }
    if(!glm && learner=="pspline"){
      ctrl <- boost_control(nu=nu,mstop=mstop)
      dat1.glm <- gamboost(ynew~., data=dat1,control=ctrl)
    }
    if(learner=="mars"){
      penalty <- ifelse(degree==2,3,2) #This is unneccessary since this is the default in the mars function. We follow the advice in [8] to set the cost for each basis function optimization
                                        # to be 3 in the MARS for two-way interaction models. (see Friedman 1991, and Yi Lin and Helen Zhang, AOS, 2006)
      if(tuningSwitch){
        dat1.glm <- mars(x=as.data.frame(x),y=ynew,degree=degree,penalty=penalty)
        nk <- dat1.glm$nk
      }
      else   dat1.glm <- mars(x=as.data.frame(x),y=ynew,degree=degree,penalty=penalty,nk=nk,prune=FALSE)   #should I add "prune=FALSE" here
      dat1.glm <- mars.to.earth(dat1.glm)
    }
    if(learner=="tree"){
      datfit <- as.data.frame(list(y=ynew,x=x))
      colnames(datfit) <- c("y",colnames(x));
      dat1.glm <- gbm(y~.,data=datfit,distribution="gaussian",interaction.depth=degree,shrinkage=nu,cv.folds=nfold,n.trees=best.iter,verbose=FALSE)
###############boosting trees################
      if(tuning && cv && tuningSwitch){
        best.iter <- gbm.perf(dat1.glm,plot.it=FALSE,method="cv")
        dat1.glm <- gbm(y~.,data=datfit,distribution="gaussian",interaction.depth=degree,shrinkage=nu,n.trees=best.iter,verbose=FALSE)
      }
      if(trace) cat("k=",k," number of trees utilized",best.iter,"\n")
      mselect[k] <- ifelse(!cv,mstop,best.iter)  #save the tuning parameter
    }
#    if(learner=="acosso"){
#      best.pow <- wt.pow
#      dat1.glm <- acosso(X=as.data.frame(x),y=ynew,order=degree,wt.pow=best.pow)  #use 2 from Section 5, line 4 of Adaptive COSSO; 0 for COSSO
#    }
    if(learner=="enet"){
      dat1.glm <- enet(x=as.matrix(x), y=ynew, lambda=lamb)
    }
    if(learner %in% c("MCP","SCAD")){
      if(tuning){
      cv.fit <- cv.ncvreg(X=as.matrix(x), y=ynew, penalty=learner, gamma=gamma, nfolds=nfold)
      cv.min <- which.min(cv.fit$cve)
      ### needs something with fixed lambda
      }
      else cv.min <- whichlambda
      dat1.glm <- ncvreg(X=as.matrix(x), y=ynew, penalty=learner, gamma=gamma, lambda=lambda)
    }
######################################################################################
### Choose tuning parameter for boosting methods, NOT for learner "enet", "mars","acosso" ####
###################################################################################### 
    if(learner%in% c("linear.regression","pspline") && tuning){
      if(!cv){
        aic <- AIC(dat1.glm,method=method,df=df) #model selection
        if(trace)
          cat("aic=",aic, "mstop(aic)=",mstop(aic),"\n")
        if(mstop(aic)==1)
          mselect[k] <- mstop
        else mselect[k] <- mstop(aic)
        dat1.glm <- dat1.glm[mselect[k]]
      }
      else {
        if(learner %in% c("linear.regression","pspline") ){                      
          file <- NULL
          mselect[k] <- nmstop <- cvmboost(obj=dat1.glm,nrow(dat1),nfold=nfold,figname=file);
          cat("CV selects mstop", nmstop,"\n")
          dat1.glm <- dat1.glm[nmstop]
        }
      }
    }
### Compute predicted values and convergence criteria                                    
    if(glm && learner=="linear.regression"){
      beta0bj <- coef(dat1.glm)[1] + dat1.glm$offset
      betabj <- coef(dat1.glm, which = 1:length(variable.names(dat1.glm)))[-1]
      pred.bj <- predict(dat1.glm)
      ybst[k,] <- Fboost <- pred.bj
      b[k,] <- betabj
      bdiff <- 1000
      if(k > 1)
        bdiff <- sum((b[k,] - b[k-1,])^2)
    }
    else {
      if(learner=="enet"){
        if(!is.null(valdata)) pred.bj <- predict(dat1.glm, newx=valdata[,-(1:2)], type="fit", s=s, mode="fraction")$fit
        ybst[k,] <- Fboost <- predict(dat1.glm, x, type="fit", s=s, mode="fraction")$fit
      }
      else if(learner %in% c("MCP", "SCAD")){
        if(!is.null(valdata)) pred.bj <- predict(dat1.glm, X=valdata[,-(1:2)], type="response", which= cv.min)
        ybst[k,] <- Fboost <- predict(dat1.glm, X=x, type="response", which=cv.min)
     } 
      else if(learner!="tree" && learner!="acosso"){
        if(!is.null(valdata)) pred.bj <- predict(dat1.glm, newdata=valdata[,-(1:2)])
        ybst[k,] <- Fboost <- predict(dat1.glm)
      }
      else
        if(learner=="tree"){
### predicted values on testing (validated) data
          if(!is.null(valdata))  pred.bj <- predict(dat1.glm, newdata=as.data.frame(valdata[,-(1:2)]),n.trees=dat1.glm$n.trees)
          ybst[k,] <- Fboost <- dat1.glm$fit   #predicted values on training data
        }
        else{
         # if(!is.null(valdata)) pred.bj <- predict.acosso(obj=dat1.glm, X.new=as.data.frame(valdata[,-(1:2)]))
         # ybst[k,] <- Fboost <- dat1.glm$y.hat
        }
      if(!is.null(valdata)) mse[k] <- mean((valdata[,1] - pred.bj)^2)
    }    
    if(k>1){
      ydiff <- ybstdiff[k] <- max(abs(Fboost - ybst[k-1,]))
      ybstcon[k] <- sum((ybst[k,]-ybst[k-1,])^2)/sum(ybst[k-1,]^2)
    }
    mseun[k] <- mean((Fboost-ynew)^2) #mse
    if(k >1 && trace)
      cat("    k=",k,"   ybstdiff", ybstdiff[k],"  ybstcon", ybstcon[k]," mse", mse[k],"\n")
### Check convergence status
    if(!nonconv){
      if(k > 1)
        if((glm && learner=="linear.regression" && bdiff <= tol)
           || (learner %in% c("enet","MCP","SCAD") && ybstcon[k] <= tol)
           || (!glm && ybstcon[k] <= tol)){
          contype <- 0
          break
        }
        else if(k >= iter.bj) {
          cycleydiff <- NULL
          if(glm && learner=="linear.regression") {
            cycle.coef.bj <- NULL
            firstb <- betabj
          }
          nonconv <- TRUE
          firstydiff <- ydiff
          first.ybstcon <- ybstcon[k]
          first.dat1.glm <- dat1.glm
          cycleb <- NULL
          tuningSwitch <- FALSE;
          if(learner=="tree")  best.iter <- best.iter                           #use the best number of trees of the previous model
          else if(learner=="mars") ynew <- dat1.glm$ynew
        }
    }
    else {
      if(glm && learner=="linear.regression"){
        coef.bj <- coef(dat1.glm, which = 1:length(variable.names(dat1.glm)), off2int=TRUE)
        cycle.coef.bj <- rbind(cycle.coef.bj,coef.bj)
      }
      cycleydiff <- c(cycleydiff,ydiff)
      cycleperiod <- cycleperiod + 1
      if(glm && learner=="linear.regression"){
        if(sum((firstb - coef.bj[-1])^2) < tol || ydiff <= tol){
          contype <- 1
          break
        }
        else if(cycleperiod >= max.cycle){
          contype <- 2
          break
        }
      }
      else {
        if(abs(ybstcon[k]-first.ybstcon) < tol){
          contype <- 1
          break
        }
        else if(cycleperiod >= max.cycle){
          contype <- 2
          break
        }
      }
    }
    k <- k + 1
  }
### END of BJ iteration loop
  if(trace)
    cat("\ncycle period is",cycleperiod,"\n")
  if(contype==2)
    dat1.glm <- first.dat1.glm                                    
  if(all(!is.na(Fboost))){
    tmpy <- y[cens==1]- Fboost[cens==1]
    tmpx <- (y[cens==1] + Fboost[cens==1])/2
    mse.bj <- mean(tmpy^2)
    if(glm && learner=="linear.regression"){
      beta0bj <- coef(dat1.glm)[1] + dat1.glm$offset
      betabj <- coef(dat1.glm, which = 1:length(variable.names(dat1.glm)))[-1]
    }
### prediction and MSE with censoring data in real applications
    if(!is.null(valdata)){
      if(glm && learner=="linear.regression"){
        pred.bj <- as.matrix(valdata)[,-(1:2)] %*% (c(beta0bj,betabj)/c(1,normx))
        mse.bj.val <- mean((valdata[,1][valdata[,2]==1] - pred.bj[valdata[,2]==1])^2)
      } else {
        if(learner!="tree" && learner!="acosso" && ! learner %in% c("enet","MCP","SCAD"))
          pred.bj <- predict(dat1.glm, newdata=valdata[,-(1:2)])
        else 
          if(learner=="tree"){
            pred.bj <- predict(dat1.glm, newdata=as.data.frame(valdata[,-(1:2)]),n.trees=dat1.glm$n.trees)
          }
#          else if(learner=="acosso")
#            pred.bj <- predict.acosso(X.new=valdata[,-(1:2)],obj=dat1.glm)
          else if(learner=="enet")
            pred.bj <- predict(dat1.glm, newx=valdata[,-(1:2)], s=s, type="fit", mode="fraction")$fit
          else if(learner %in%c("MCP", "SCAD"))
            pred.bj <- predict(dat1.glm, X=valdata[,-(1:2)], type="response", which=cv.min)
        mse.bj.val <- mean((valdata[,1][valdata[,2]==1] - pred.bj[valdata[,2]==1])^2)
      }
    }
### prediciton and MSE with simulations
    if(learner=="enet"){
      tmp <- predict(dat1.glm, type="coef", s=s, mode="fraction")$coef
      beta0.enet <- mean(ynew) - apply(x[,dat1.glm$allset], 2, mean) %*% tmp ### intercept
      beta.enet <- rep(0, p)
      beta.enet[dat1.glm$allset] <- tmp
    }
    else if(learner %in% c("MCP", "SCAD"))
      coef.ncv <- predict(dat1.glm, X=x, type="coefficients", which=cv.min)
    if(!glm && !is.null(valdata)){
      if(learner!="tree" && learner!="acosso" && ! learner %in% c("enet","MCP","SCAD"))
      #if(learner!="tree" && learner!="acosso" && learner!="enet")
        pred.bj <- predict(dat1.glm, newdata=valdata[,-(1:2)])
      else
        if(learner=="tree"){
          pred.bj <- predict(dat1.glm, newdata=as.data.frame(valdata[,-(1:2)]),n.trees=dat1.glm$n.trees)
        } #else if(learner=="acosso") 
#          pred.bj <- predict.acosso(X.new=valdata[,-(1:2)],obj=dat1.glm)
      mse.bj.val <- mean((pred.bj - valdata[,1])^2) 
    }
    if(trace) {
      cat("mse.bj=",mse.bj,"\n","correlation of predicted and observed times in noncensoring training data is",cor(y[cens==1],Fboost[cens==1]),"\n\n")
      cat("mse of predicted times of validate data is\n")
      cat("mse.bj.val",mse.bj.val,"\n")
    }
    coef.bj <- NA
    if(glm && learner=="linear.regression"){
      coef.bj <- c(beta0bj, betabj)
      coef.bj <- coef.bj/c(1,normx)
      nz.bj <- sum(abs(coef(dat1.glm)[-1])>0) # -1 means intercept is not counted
      if(trace) {cat("Number of Non-zero coefficients with BJ boosting excluding intercept is",nz.bj,"\n")
                 print(coef.bj[abs(coef.bj)>0])
               }
    }
  }
                                        #ynew is the final imputed response used for boosting
  cycle.coef.diff <- NA
  if(exists("cycle.coef.bj"))
    cycle.coef.diff <- max(abs(scale(cycle.coef.bj, coef.bj, FALSE))) #compute max of absolute difference between coef in the cycle to the one claimed as the solution
  vim <- 0; interactions=NULL
  d <- ncol(x)
                                        #construct a one-to-one corespondence table between input column and the ensemble index
  ind <- matrix(NA,ncol=2,nrow=d+d*(d-1)/2)
  kk <- 1
  for(i in 1:d)
    for(j in i:d){
      ind[kk,1] <- i; ind[kk,2] <- j
      kk <- kk + 1
    }
  if(rel.inf && learner=="tree" && degree > 1){
    vim <- summary(dat1.glm,plotit=FALSE,order=FALSE)[,2]
    interactions <- vim.interactions(dat1.glm,pred.data=x,x.pair=subset(ind,ind[,1]!=ind[,2]),learner="tree",verbose=trace)
  }
### Compute which variable is selected
  if(learner=="tree"){
    xselect <- summary(dat1.glm,order=FALSE,plotit=FALSE)[,2]
    xselect <- ifelse(xselect > 0, 1, 0)   #variable selected if xselect=1, o.w. 0.
  }
  else if(learner=="mars"){
    vim <- evimp(update(dat1.glm),sqrt.=TRUE,trim=FALSE)[,c(1,4)]
    #  	  vim <- evimp(update.earth(dat1.glm),sqrt=TRUE,trim=FALSE)[,c(1,4)]
    vim <- vim[order(vim[,1]),]
    vim <- vim[,2]
    xselect <- ifelse(vim > 0, 1, 0)   #variable selected if xselect=1, o.w. 0.
  }
  else if(learner=="linear.regression")
    xselect <- ifelse(abs(coef.bj[-1]) > 0, 1, 0)   #without the intercept
  else if(learner=="enet"){
    tmp <- predict(dat1.glm, type="coef", s=s, mode="fraction")$coef
    xselect <- ifelse(abs(tmp) > 0, 1, 0) 
  }
  else if(learner %in% c("MCP", "SCAD"))
    xselect <- ifelse(abs(coef.ncv) > 0, 1, 0)
  else if(learner=="pspline"){
    xselect <- rep(0,dim(x)[2])
    tmp <- which(colnames(x)%in%colnames(dat1.glm$data$input)[unique(dat1.glm$ensemble)]) 
    xselect[tmp] <- 1
  }
  else if(learner=="acosso"){
    if(dat1.glm$order==1)
      xselect <- ifelse(dat1.glm$theta > 0, 1, 0)
    else{
      ind <- gen.ind(p,learner="acosso") 
      xselect <- rep(0,dim(x)[2])
      tmp <- unique(as.vector(ind[dat1.glm$theta > 0,]))
      xselect[tmp] <- 1
    }}   
### Compute variable importance and interaction measures for MARS
  if(learner=="mars" && degree > 1 && vimpint){
    ind <- gen.ind(p) 
    interactions <- vim.interactions(dat1.glm,pred.data=x,x.pair=subset(ind,ind[,1]!=ind[,2]),learner="mars",verbose=trace)
  }
  if(learner=="tree")
    vim <- summary(dat1.glm,plotit=FALSE,order=FALSE)[,2]
  mse.tr <- NULL
  if(learner=="enet" && tuning){
    mse.tr <- sum((ynew - predict(dat1.glm, x, s=s, mode="frac", type="fit")$fit)^2)
    b <- predict(dat1.glm, type="coef", s=s, mode="frac")$coef
    if(any(abs(b) > 0)){
      b <- which(abs(b) > 0)
      x0 <- as.matrix(x[,b])
      if(lamb==0) q <- dim(x0)[2]
      else {q <- sum(diag(x0 %*% solve(t(x0) %*% x0 + diag(lamb, nrow=dim(x0)[2])) %*% t(x0))) ### trace
          }
    }
    else q <- 0
    mse.tr <- mse.tr/(length(ynew) - q)^2 
  }
  else if(learner=="enet") 
    coef.bj <- c(beta0.enet, beta.enet)
  else if(learner %in%c("MCP", "SCAD"))
    coef.bj <- coef.ncv 
  
  RET <- list(x=x,y=y,cens=cens,ynew=ynew,res.fit=dat1.glm,learner=learner,degree=degree,mse=mse,nz.bj.iter=nz.bj.iter,mse.bj=mse.bj,mse.bj.val=mse.bj.val,nz.bj=nz.bj,mse.all=mseun[1:(k-1)],yhat=Fboost,ybstdiff=c(NA,ybstdiff[1:(k-1)]),ybstcon = ybstcon,coef.bj=coef.bj,pred.bj=pred.bj,cycleperiod=cycleperiod,cycle.coef.diff = cycle.coef.diff,nonconv=nonconv,fnorm2=fnorm2,vim=100*vim/sum(vim),interactions=interactions,mselect=mselect,contype=contype,xselect=xselect,lamb=lamb, s=s, mse.tr=mse.tr,valdata=valdata)
  RET$call <- match.call()
  class(RET) <- "bujar"
  return(RET)
}

gen.ind <- function(d,learner="tree"){
  ind <- matrix(NA,ncol=2,nrow=d+d*(d-1)/2)
  if(learner=="mars"){
                                        #construct a one-to-one corespondence table between input column and the ensemble index
    kk <- 1
    for(i in 1:d)
      for(j in i:d){
        ind[kk,1] <- i; ind[kk,2] <- j
        kk <- kk + 1
      }
  }
  else if(learner=="tree" || learner=="acosso"){       #cf: get.gram function in acosso.R 
    ind[1:d,] <- cbind(1:d,1:d)
    next.ind <- d+1
    for(i in 1:(d-1))
      for(j in ((i+1): d)){
        ind[next.ind,] <- cbind(i,j)
        next.ind <- next.ind + 1
      }
  }
  ind
}

convbujar <- function(x){
  ybstdiff <- x$ybstdiff
  ybstcon <- x$ybstcon
  mseun <- x$mse.all
  mse <- x$mse
  fnorm2 <- x$fnorm2
  plot(ybstcon, type="b",xlab="Buckley-James estimator iteration",ylab="Convergence criterion",ylim=c(0,0.01))
}

###compute the number of covariates selected based on the position of x
nxselect <- function(obj, varpos) sum(obj$xselect[varpos] == 1)

print.bujar <- function(x, ...) {

  cat("\n")
  cat("\t Models Fitted with Buckley-James Regression\n")
  cat("\n")
  if (!is.null(x$call))
    cat("Call:\n", deparse(x$call), "\n\n", sep = "")
  cat("\n")
  if(x$learner%in%c("linear.regression","mars","pspline","tree"))
    cat("Base learner: ", x$learner, "\n")
  else cat("Regression methods: ", x$learner, "\n") 
  cat("\n")
  if(x$learner=="linear.regression"){
    cat("Coefficients: \n")
    cf <- x$coef.bj
    print(cf)
    cat("\n")
  }
  invisible(x)
}

### methods: coefficients
coef.bujar <- function(object, ...) {
  if(!object$learner %in% c("linear.regression","pspline","enet"))
    stop("Coefficients Not implemented for learner ",object$learner,"\n")
  if(object$learner %in%  c("linear.regression","pspline"))
    object$coef.bj
  else predict(object, type="coef")
}

plot.bujar <- function(x, ...){
  if(!x$learner %in% c("mars", "pspline", "acosso"))
    plot(x$res.fit)
  else stop("Not implemented for learner ",x$learner,"\n")
}

predict.bujar <- function(object, ...){
  if(object$learner != "acosso")
    predict(object$res.fit, ...)
#  else
    #predict.venus(object$x, object$res.fit)
}

summary.bujar <- function(object, ...)
  summary(object$res.fit, ...)
