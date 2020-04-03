generate_data <- function(n,true.alpha=c(1, 1, 1)/sqrt(3)){
  X = matrix(runif(length(true.alpha)*n), ncol=length(true.alpha))
  U = X%*%true.alpha

  #Z<-rep(c(0,1),n/2)
  Z<-rnorm(n)
  q <- NULL
  py <- NULL
  for(i in 1:n){
    q[i] = sin( (U[i]-0.3912)*pi/(1.3409 -0.3912) ) + 0.3*Z[i]
    py[i] = exp(q[i])/(1+exp(q[i]))
  }

  y= rbinom(n, size=1, prob=py)
  return(list("X" = X, "Y" = y,"Z"=Z))
}


si <- function(alpha,y,x,z,opt=TRUE,k=10,fam, fx=FALSE) {
  ## Fit single index model using gam call. Return ML is opt==TRUE
  ## and fitted gam with theta added otherwise...

  #reparameterized theta
  alpha_add1 <- c(1, alpha)
  theta <- alpha_add1/sqrt(sum(alpha_add1^2))

  a <- x%*%theta
  b <- mgcv::gam(y~s(a,fx=fx,k=k)+z-1,family=fam,method="ML")
  if (opt) return(b$gcv.ubre) else {
    b$theta <- theta  ## add theta
    return(b)
  }
}      # end of defining si;

si2 <- function(alpha,y,x,z,opt=TRUE,k=10,fam, fx=FALSE) {
  ## Fit single index model using gam call. Return ML is opt==TRUE
  ## and fitted gam with theta added otherwise...
  theta <- alpha
  theta <- theta*sign(theta[1])/sqrt(sum(theta^2))
  a <- x%*%theta
  b <- gam(y~s(a,fx=fx,k=k)+z,family=fam,method="ML")
  if (opt) return(b$gcv.ubre) else {
    b$theta <- theta  ## add theta
    return(b) }
}      # end of defining si;



gplsim<- gplsimPs <- function(Y=Y,X=X,Z=Z,family=binomial,penalty=TRUE,user.init=FALSE){
  p <- dim(X)[2]
  if(identical(FALSE,user.init)){
    if(length(user.init)!=(p-1)){
      stop("user.init length must be p-1")
    }else{
      init.alpha <- user.init
    }
  }else{
    er_np <- optim(rep(0,p-1), si, y=Y,x=X,z=Z, fam=family, hessian=TRUE, fx=TRUE)
    init.alpha <- er_np$par
  }

  er <- optim(init.alpha,si,y=Y,x=X,z=Z,fam=family,hessian=TRUE,fx=!penalty)
  b <- si(er_np$par,y=y,X,Z, fam=family, opt=FALSE)
  return(b)
}
