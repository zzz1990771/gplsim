#internal support function
#function to fit plsim with iterative procedure, through non-linear least square
#input x,y,degree,nknots,alpha0(if different from default OLSe),maxiter,z
plsiml <- function(y,x,z=NULL,degree=3,nknots=12,alpha0=NULL,maxiter=2){

  if(is.null(z)==TRUE){
    temp_XZ <- x
  }else{
    temp_XZ <- cbind(x,z)
  }

  if(is.null(alpha0)){
    alpha0 <- coef(lm(y~temp_XZ))[2:(1+NCOL(x))]# OLS estimate for initial alpha0
    alpha <- sign(alpha0[1])*alpha0/sqrt(sum(alpha0^2))
  }else{
    alpha <- alpha0
  }

  xnew <- x%*%alpha
  srspl_fit <- srspl(xnew,y,degree,nknots,z)
  yhat = srspl_fit$yhat
  beta0 = srspl_fit$beta
  yhat1 = srspl_fit$yhat1
  lambda = srspl_fit$lambda

  param0 <- c(alpha,beta0)

  n <- length(y)
  m <- length(param0)
  d <- length(alpha)
  d2 <- NCOL(temp_XZ)-d

  param <- param0
  paramold <- 2 * param
  mse <- 1
  iter <- 1
  tol <- 0.000001

  while( (max(abs(paramold-param)) > tol) & (iter < maxiter)){

    paramold <- param

    alpha <- param[1:d]
    alpha <- sign(alpha[1])*alpha/sqrt(sum(alpha^2))# normalize alpha
    beta <- param[(d+1):m]

    xnew <- x%*%alpha

    srspl_fit <- srspl(xnew,y,degree,nknots,z)
    yhat = srspl_fit$yhat
    beta = srspl_fit$beta
    yhat1 = srspl_fit$yhat1
    lambda = srspl_fit$lambda

    mse <- (1/n)*sum((y-yhat)^2)

    #options = optimset('LargeScale','on','MaxFunEvals','200','TolFun',1e-6);
    #plresim(param,x=x,y=y,z=z,degree=degree,nknots=nknots,lambda=lambda)
    nls.fit <- minpack.lm::nls.lm(param,fn = plresim,x=x,y=y,z=z,degree=degree,nknots=nknots,lambda=lambda,control = list(maxiter = 500))
    param <- nls.fit$par
    # alpha is not normalized by lsqnonlin

    iter <- iter + 1
  }

  alpha <- param[1:d]
  alpha <- sign(alpha[1])*alpha/sqrt(sum(alpha^2))# normalize alpha
  xnew <- x%*%alpha

  srspl_fit <- srspl(xnew,y,degree,nknots,z)
  yhat = srspl_fit$yhat
  beta = srspl_fit$beta
  yhat1 = srspl_fit$yhat1
  lambda = srspl_fit$lambda

  mse <- 1/n*sum((y-yhat)^2)
  res <- NULL
  res$param <- param
  res$mse <- mse
  res$lambda <- lambda
  res$iter <- iter
  res$alpha <- alpha
  return(res)
}

##internal support function
#plresim is the residual vector of length n
plresim <- function(param,x,y,z,degree,nknots,lambda){

  m <- length(param)
  d <- NCOL(x)
  alpha <- param[1:d]
  alpha <- sign(alpha[1])*alpha/sqrt(sum(alpha^2))
  beta <- param[(d+1):length(param)]



  xnew <- x%*%alpha # reduce x dimentionality by this projection.
  knots <- quantileknots(xnew,nknots)
  xm <- powerbasis(xnew,degree,knots)
  exty <- c(y,  matrix(0,nrow=nknots,ncol=1))

  if(is.null(z)){
    yhat <- xm%*%beta
  }else{
    yhat <- cbind(xm,z)%*%beta
  }
  g <- exty - c(yhat, sqrt(lambda)*beta[(degree+2):(degree+nknots+1)])
  return(g)
}

#internal support function
quantileknots <- function(x,nknots){
  boundstab <- 0
  n <- length(x)
  xsort <- sort(x)
  loc <- t(n*(1:nknots+2*boundstab)) / (nknots+1+2*boundstab)
  knots <- xsort[round(loc)]
  knots <- knots[(1 + boundstab) : (nknots + boundstab)]
  return(knots)
}

#internal support function
powerbasis <- function(x,degree,knots,der=0){
  n <- NROW(x)
  nknots <- length(knots)

  if (der == 0 ){
    xm <- rep(1,n)
  }else{
    xm <- rep(0,n)
  }

  for (i in 1:degree ){
    if (i < der ){
      xm <- cbind(xm, rep(0,n))
    }else{
      xm <- cbind(xm, prod((i-der+1):i) * x^(i-der))
    }
  }

  if (nknots > 0 ){
    for (i in 1:(nknots) ){
      xm <- cbind(xm,prod((degree-der+1):degree) * (x-knots[i])^(degree-der)*(x > knots[i]))
    }
  }
  xm
}

#internal support function
pbderiv1 <- function(x,degree,knots){
  n <- NROW(x)
  nknots <- length(knots)

  xm1 <- rep(0,n)

  for (i in 1:degree ){
    xm1 <- cbind(xm1,i*x^(i-1))
  }

  if (nknots > 0 ){
    for (i in 1:(nknots) ){
      xm1 <- cbind(xm1, degree*((x-knots[i])^(degree-1))*(x > knots[i]))
    }
  }
  return(xm1)
}

#internal support function
srspl <- function(x,y,degree=3,nknots=10,z=NULL,penwt=(10^seq(-6, 7, length.out=n))){
  n = NROW(x)
  #ngrid = 200
  #ugrid = t(seq(0.4,1.4,length.out=ngrid))
  #knots = quantileknots(ugrid,nknots);
  knots = quantileknots(x,nknots)
  n = length(x) ;
  xm=rep(1,n);
  xm = powerbasis(x,degree,knots) ;
  xm = cbind(xm,z);
  xm1= pbderiv1(x,degree,knots);


  xx = t(xm)%*%xm ;
  xy = t(xm)%*%y ;
  nz = (!is.null(z))*NCOL(z)
  id = diag(c(rep(0,(degree+1)),rep(1,(nknots)),rep(0,(nz)))) ;
  m = length(penwt) ;
  beta = matrix(0,NCOL(xm),m);
  yhat = matrix(0,n,m) ;
  asr = matrix(0,m,1) ;
  gcv = asr ;
  trsd = asr ;
  ssy = t(y)%*%y ;
  for (i in 1:m ){
    #return(xx+ penwt[i]*id)
    b <- solve(xx + penwt[i]*id,tol = 0)
    xxb <- xx%*%b
    trsd[i] <- sum(diag(xxb))
    beta[,i] <- b %*% xy
    asr[i] <-  (ssy - 2*t(xy)%*%beta[,i] + t(beta[,i])%*%xx%*%beta[,i])/n
    if (i==1){
      trsdsd <- sum(diag(xxb*xxb))
      sigma2hat <-  n*asr[i]/(n-2*trsd[i]+trsdsd)
    }
    gcv[i] <- asr[i] / ((1-trsd[i]/n)^2 )
  }

  imin <- which(  (gcv==min(gcv)) )
  df <-  trsd
  dfs <- trsd[imin]
  sigma2hat <- n*asr[imin] / (n - df[imin])
  b <- solve(xx + penwt[imin]*id,tol = 0)
  beta <- beta[,imin]
  yhat <- xm%*%beta
  yhat1 <- xm1%*%beta[1:(1+degree+nknots)] # estimate of the first derivative of y
  postvarbeta <- sigma2hat*b
  postvaryhat <- sigma2hat*(xm*(xm%*%b))%*%matrix(1,length(beta),1)
  lambda <- penwt[imin]
  res <- NULL
  res$yhat = yhat
  res$beta = beta
  res$yhat1 = yhat1
  res$gcv =gcv
  res$imin = res$imin
  res$lambda = lambda
  return(res)
}













