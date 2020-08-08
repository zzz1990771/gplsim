generate_data <- function(n,true.alpha=c(1, 1, 1)/sqrt(3),binary=TRUE){
  sigma = 0.1
  c1 = 0.3912
  c2 = 1.3409
  rho = 0.3

  X = matrix(runif(length(true.alpha)*n), ncol=length(true.alpha))
  U = X%*%true.alpha

  #Z<-rnorm(n)
  Z <- 1 - c(1:n)%%2
  q = sin( (U-c1)*pi/(c2 -c1) ) + rho*Z
  if(binary){
    py = exp(q)/(1+exp(q))
    y = rbinom(length(q), size=1, prob=py)
  }else{
    y = q + rnorm(length(q),0,sigma)
  }

  return(list("X" = X, "Y" = y,"Z"=Z))
}


si <- function(alpha,y,x,z,opt=TRUE,k=10,smooth_selection,fam, fx=FALSE) {
  ## Fit single index model using gam call. Return ML is opt==TRUE
  ## and fitted gam with theta added otherwise...

  #reparameterized theta
  theta <- sign(alpha[1])*alpha/sqrt(sum(alpha^2))

  a <- x%*%theta
  b <- mgcv::gam(y~s(a,bs="tp",fx=fx,k=k)+z,family=fam,method= smooth_selection )
  if (opt) return(b$deviance) else {
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



gplsim<- gplsimPs <- function(Y=Y,X=X,Z=Z,k=10,family=binomial,smooth_selection = "GCV.Cp",penalty=TRUE,user.init=FALSE){
  #if p <= 2
  p <- dim(X)[2]
  #op_method <- ifelse(p>2,"Nelder-Mead","Brent")
  op_method <- "Nelder-Mead"
  if(p<2){stop("SIM predictos must no less than two")}
  if(!identical(FALSE,user.init)){
    if(length(user.init)!=(p)){
      stop("user.init length must be p")
    }else{
      init.alpha <- user.init
    }
  }else{
    er_np <- suppressWarnings(optim(rep(1,p), si, y=Y,x=X,z=Z, fam=family,  smooth_selection=smooth_selection, hessian=TRUE, fx=TRUE, method=op_method,k=k))
    init.alpha <- er_np$par
  }

  er <- suppressWarnings(optim(init.alpha,si,y=Y,x=X,z=Z,fam=family,  smooth_selection=smooth_selection, hessian=TRUE,fx=!penalty, method=op_method,k=k))
  b <- si(er$par,y=y,X,Z, fam=family, smooth_selection=smooth_selection, opt=FALSE,k=k)
  return(b)
}


plot.si <- function(model_obj,index=NULL,...){
  if(is.null(index)){
    mgcv::plot.gam(model_obj,xlab = "single index",ylab = "mean",...)
  }else{
    UQ = cbind(model_obj$model$a,result$fitted.values)
    UQ0 = UQ[index==0,]
    UQ1 = UQ[index==1,]
    y_scale = range(result$fitted.values)
    plot(UQ0[order(UQ0[,1]),],type="l",ylim=y_scale,lty=2,...)
    lines(UQ1[order(UQ1[,1]),],type="l",lty=2,...)
  }
}
