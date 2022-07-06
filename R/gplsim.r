#' Data generation function for simulation and demonstration
#' A sine-bump setting has been employed.
#'
#' @param n sample size
#' @param true.theta true single-index coefficients,
#' default is c(1,1,1)/sqrt(3) for setting 1 and c(1,2)/sqrt(5) for other settings
#' @param family chose from "gaussian", "binomial" or "poisson".
#' @param ncopy generates multiple copies of data for Monte Carlo simulations
#'
#' @return X single index predictors
#' @return Y response variables, a list
#' @return Z partial linear predictor(s)
#' @return single_index_values single index term
#' @export
generate_data <- function(n,true.theta=c(1, 1, 1)/sqrt(3),family="gaussian",ncopy=1){
  #parameter setting
  sigma = 0.1
  c1 = 0.3912
  c2 = 1.3409
  rho = 0.3

  X = matrix(stats::runif(length(true.theta)*n), ncol=length(true.theta))
  U = X%*%true.theta
  fU = sin( (U-c1)*pi/(c2 -c1) )

  #Z<-rnorm(n)
  Z <- 1 - c(1:n)%%2
  q = fU + rho*Z

  if(family=="gaussian"){
    ylist <- lapply(vector(mode = "list", length = ncopy),function(x){q + rnorm(length(q),0,sigma)})
  }else if(family=="binomial"){
    py = exp(q)/(1+exp(q))
    ylist <- lapply(vector(mode = "list", length = ncopy),function(x){rbinom(length(q), size=1, prob=py)})
  }else if(family=="poisson"){
    py = exp(q)
    ylist <- lapply(vector(mode = "list", length = ncopy),function(x){rpois(length(q),py)})
  }
  if(ncopy==1){ylist = ylist[[1]]}
  return(list("X" = X, "Y" = ylist,"Z"=Z,"single_index_values"=fU))
}

#' An internal function to optimization and fitting. 
#' Don't use it solely.
#'
#' @param alpha single-index coefficients
#' @param y Response variable, should be a vector.
#' @param x Single index covariates.
#' @param z Partially linear covariates.
#' @param opt see ?gplsim
#' @param k see ?gplsim
#' @param smooth_selection see ?gplsim
#' @param fam see ?gplsim
#' @param bs see ?gplsim
#' @param fx see ?gplsim
#' @param scale see ?gplsim
#' 
#' @return b fitted gam object
si <- function(alpha,y,x,z,opt=TRUE,k=13,smooth_selection,fam,bs="ps", fx=FALSE,scale=scale) {
  ## Fit single index model using gam call. Return ML is opt==TRUE
  ## and fitted gam with theta added otherwise...

  #reparameterized theta
  theta <- sign(alpha[1])*alpha/sqrt(sum(alpha^2))

  a <- x%*%theta
  if(is.null(z)){
    b <- mgcv::gam(y~s(a,bs=bs,fx=fx,m=2,k=k),family=fam,method= smooth_selection,scale=scale )
  }else{
    b <- mgcv::gam(y~s(a,bs=bs,fx=fx,m=2,k=k)+z,family=fam,method= smooth_selection,scale=scale )
    b$gamma <- b$coefficients[c(2:(1+NCOL(z)))]
  }
  if (opt) return(b$deviance) else {
    b$theta <- theta  ## add theta
    #b$gamma <- b$coefficients[c(2:(1+NCOL(z)))]
    class(b) <- c("gplsim","gam","glm","lm")
    return(b)
  }

}

#' Function to fit generalized partially linear single-index models via penalized splines
#'
#' This function employs penalized spline (P-spline) to estimate generalized
#' partially linear single index models, which extend the generalized linear
#' models to include nonlinear effect for some predictors.
#'
#' @param Y Response variable, should be a vector.
#' @param X Single index covariates.
#' @param Z Partially linear covariates.
#' @param family A \code{family} object: a list of functions and expressions
#' for defining \code{link} and \code{variance} functions.
#' Families supported are \code{binomial}, \code{gaussian}.
#' The default family is \code{gaussian}.
#' @param penalty Whether use penalized splines or un-penalized splines to fit the model. The default is TRUE.
#' @param penalty_type The optional argument penalty_type is a character variable, which specifies the type of penalty used in the penalized splines estimation. The default penalty type is {L}_{2} penalty, while {L}_{1} is also supported.
#' @param scale The optional argument scale is a numeric indicator with a default value set to -1. Any negative value including -1 indicates that the scale of response distribution is unknown, thus need to be estimated. Another option is 0 signaling scale of 1 for Poisson and binomial distribution and unknown for others. Any positive value will be taken as the known scale parameter.
#' @param smooth_selection The optional argument smooth_selection is another character variable that specifies the criterion used in the selection of a smoothing parameter. The supported criteria include "GCV.Cp","GACV.Cp", "ML","P-ML", "P-REML" and "REML", while the default criterion is "GCV.Cp".
#' @param profile profile is a logical variable that indicates whether the algorithm with profile likelihood or algorithm with NLS procedure should be used. The default algorithm is set to algorithm with profile likelihood.
#' @param bs bs is a character variable that specifies the spline basis in the estimation of unknown univariate function of single index. Default is P-splines.
#' @param user.init The user.init is a numeric vector of the same length as the dimensionality of single index predictors. The users can use this argument to pass in any appropriate user-defined initial single-index coefficients based on prior information or domain knowledge. The default value is NULL
#' @param k k is the the dimension of the basis used to represent the smooth term. The default is set at 13.
#' 
#' @return theta  Estimation of Theta
#' @return coefficients  the coefficients of the fitted model. Parametric coefficients are first, followed by coefficients for each spline term in turn.
#' @return ... See GAM object
#' @examples
#' # parameter settings
#' n=200
#' true.theta = c(1, 1, 1)/sqrt(3)
#' # Gaussian case
#' # This function generate a plain sin bump model with gaussian response.
#' data <- generate_data(n,true.theta=true.theta,family="gaussian")
#' y=data$Y       # continous response
#' X=data$X       # single index term ;
#' Z=data$Z       # partially linear term ;
#'
#' result <- gplsim(y,X,Z,family = gaussian)
#' result$theta
#' result$coefficients
#' summary(result)
#'
#'
#' #plot the estimated single index function curve
#' plot_si(result)
#' @export
gplsim <- function(Y=Y,X=X,Z=Z,family=gaussian(),penalty=TRUE,penalty_type = "L2", scale = -1, smooth_selection = "GCV.Cp",profile = TRUE, bs="ps", user.init=NULL,k=13){

#if p <= 2
  p <- dim(X)[2]
  if(p<2){stop("SIM predictos must no less than two")}

  if(!is.null(user.init)){
    if(length(user.init)!=(p)){
      stop("user.init length must be p")
    }else{
      init.alpha <- user.init
    }
  }else{
    if(is.null(Z)){
      linear <- glm(Y~X-1,family=family)
    } else{
      temp_XZ <- cbind(X,Z)
      linear <- glm(Y~temp_XZ-1,family=family)
    }
    init.alpha <- linear$coefficients[1:(NCOL(X))]
    init.alpha <- sign(init.alpha[1])*init.alpha/sqrt(sum(init.alpha^2))
  }

  #non-profile
  #test
  if(profile==FALSE){
    non.profile.fit <- plsiml(x=X,y=Y ,z=Z,degree=3,nknots=13,maxiter = 2,alpha0 = NULL)
    a <- X%*%non.profile.fit$alpha
    if(is.null(Z)){
      b <- mgcv::gam(y~s(a,bs=bs,fx=!penalty,m=2,k=13),family=family,method= smooth_selection,scale=scale )

    }else{
      b <- mgcv::gam(y~s(a,bs=bs,fx=!penalty,m=2,k=13)+Z,family=family,method= smooth_selection,scale=scale )
      b$gamma <- b$coefficients[c(2:(1+NCOL(Z)))]
    }
    b$theta <- non.profile.fit$alpha  ## add theta

    class(b) <- c("gplsim","gam","glm","lm")


  }else{
    er <- suppressWarnings(optim(init.alpha,si,y=Y,x=X,z=Z,fam=family,  smooth_selection=smooth_selection, hessian=TRUE,fx=!penalty,k=k,bs=bs,scale=scale))
    b <- si(er$par,y=Y,X,Z, fam=family, smooth_selection=smooth_selection, opt=FALSE,k=k,bs=bs,scale=scale)
  }
  b$Xnames <- if(is.null(colnames(X)))  paste0("X.",seq(1,NCOL(X),by = 1)) else colnames(X)
  b$Znames <- if(is.null(colnames(Z)) && !is.null(Z))  paste0("Z.",seq(1,NCOL(Z),by = 1)) else colnames(Z)

  return(b)
}

#utils::globalVariables(c("yscale"))

#' Function that plot fitted curve for the unknown univariate function for single index term
#' 
#' @param x the gam/gplism fitted object
#' @param family default is gaussian()
#' @param ylab y label
#' @param yscale scale of y
#' @param plot_data controls whether to plot the data as points
#'
#' @return NULL single-index plot
#' @export
plot_si <- function(x,family=gaussian(),ylab="mean",yscale=NULL,plot_data=FALSE){
  #check the args
  #family=gaussian()
  #plot_data= ifelse(!exists("plot_data"),FALSE,plot_data)
  #ylab= ifelse(!exists("ylab"),"mean",ylab)
  
  #x$family = faussian()
  if(is.null(x$Znames)){
    offset_fit <- family$linkinv((x$linear.predictors))
  }else{
    offset_fit <- family$linkinv((x$linear.predictors-x$gamma%*%t(x$model$z)))
  }

  UQ = cbind(x$model$a,matrix(offset_fit,ncol=1))
  if(is.null(yscale)){ylim = range(offset_fit)}else{ylim = yscale}
  plot(UQ[order(UQ[,1]),],col="blue",type="l",ylim=ylim,lty=1,xlab = "single index",ylab = ylab)
  if(plot_data){graphics::points(UQ[order(UQ[,1]),1],(x$model$y)[order(UQ[,1])],pch=20)}
}

#' Summary function of gplsim object
#' 
#' @param object the gam/gplism fitted object
#' @param ... optional arguments
#'
#' @return gplsim_obj a list of summary information for a fitted gplsim object, which extends on gam object.
#' @export
summary.gplsim <- function(object,...){
  gplsim_obj <- mgcv::summary.gam(object)
  #single index term
  p.table.sim <- matrix(object$theta,ncol=1)
  dimnames(p.table.sim) <- list(object$Xnames, c("Estimate"))
  gplsim_obj$p.coeff.sim <- p.table.sim
  #partial linear term variable names
  row.names(gplsim_obj$p.table) <- c("Intercept",object$Znames)
  class(gplsim_obj) <- "summary.gplsim"
  return(gplsim_obj)
}

#' Print Summary function of gplsim object
#' 
#' @param x the gam/gplism fitted object
#' @param digits controls number of digits printed in output.
#' @param signif.stars should significance stars be printed alongside output.
#' @param ... optional arguments
#'
#' @return summarized object with nice format
#' @export
print.summary.gplsim <- function(x, digits = max(5, getOption("digits") - 3),
                              signif.stars = getOption("show.signif.stars"), ...)
  # print method for gam summary method. Improved by Henric Nilsson
{ print(x$family)
  cat("Formula:\n")

  if (is.list(x$formula)) for (i in 1:length(x$formula)) print(x$formula[[i]]) else
    print(x$formula)

  if (length(x$p.coeff)>0)
  { cat("\npartial linear coefficients:\n")
    printCoefmat(x$p.table, digits = digits, signif.stars = signif.stars, na.print = "NA", ...)
  }
  cat("\n")

  if (length(x$p.coeff.sim)>0)
  { cat("\nsingle index coefficients:\n")
    printCoefmat(x$p.coeff.sim, digits = digits, signif.stars = signif.stars, na.print = "NA", ...)
  }

  cat("\n")
  if(x$m>0)
  { cat("Approximate significance of smooth terms:\n")
    printCoefmat(x$s.table, digits = digits, signif.stars = signif.stars, has.Pvalue = TRUE, na.print = "NA",cs.ind=1, ...)
  }
  cat("\n")
  if (!is.null(x$rank) && x$rank< x$np) cat("Rank: ",x$rank,"/",x$np,"\n",sep="")
  if (!is.null(x$r.sq)) cat("R-sq.(adj) = ",formatC(x$r.sq,digits=3,width=5),"  ")
  if (length(x$dev.expl)>0) cat("Deviance explained = ",formatC(x$dev.expl*100,digits=3,width=4),"%",sep="")
  cat("\n")
  if (!is.null(x$method)&&!(x$method%in%c("PQL","lme.ML","lme.REML")))
    cat(x$method," = ",formatC(x$sp.criterion,digits=5),sep="")

  cat("  Scale est. = ",formatC(x$scale,digits=5,width=8,flag="-"),"  n = ",x$n,"\n",sep="")
  invisible(x)
}

#' function dedicated to add simulation standard error bound, in development
#' draw the bound to current plot
#' 
#' @param data a list of simulated data
#' @param family default is gaussian()
#' @param M number of simulations 
#' @param n sample size
#' @param true.theta the true coefficients
#'
#' @return NULL
#' @export
add_sim_bound <- function(data,family = gaussian(),M=200,n=1000,true.theta=c(1, 1, 1)/sqrt(3)){
  offset_fit_matrix <- matrix(0, M, n)
  for (i in 1:M){
    y=(data$Y)[[i]]       # continous response
    X=data$X       # single index term ;
    Z=data$Z       # partially linear term ;
    model_obj <- gplsim(y,X,Z,user.init=NULL,family = family)
    offset_fit <- family$linkinv(model_obj$linear.predictors-model_obj$gamma%*%t(model_obj$model$z))
    offset_fit_matrix[i,] <-  offset_fit
  }
  quan2.5=apply(offset_fit_matrix, 2, stats::quantile, prob=0.025)
  quan97.5=apply(offset_fit_matrix, 2, stats::quantile,prob=0.975)
  fit_mean=apply(offset_fit_matrix, 2, mean)
  U = data$X%*%true.theta
  graphics::lines(U[order(U)],quan2.5[order(U)],type="l",lty=2)
  graphics::lines(U[order(U)],quan97.5[order(U)],type="l",lty=2)
  graphics::lines(U[order(U)],fit_mean[order(U)],type="l",lty=1,col="blue")
  #result_table
}

#' supporting function to make tr smooth
#'
#' @param object smooth object for gam class
#' @param data the new data to predict on
#' @param knots knots
#'
#' @return tr smooth object
#' @export
smooth.construct.tr.smooth.spec<-function(object,data,knots)
  ## a truncated power spline constructor method function
  ## object$p.order = null space dimension
{ requireNamespace("splines")
  m <- object$p.order[1]
  if (is.na(m)) m <- 2 ## default
  if (m<1) stop("silly m supplied")
  if (object$bs.dim<0) object$bs.dim <- 10 ## default
  nk<-object$bs.dim-m-1 ## number of knots
  if (nk<=0) stop("k too small for m")
  x <- data[[object$term]]  ## the data
  x.shift <- mean(x) # shift used to enhance stability
  k <- knots[[object$term]] ## will be NULL if none supplied
  if (is.null(k)) # space knots through data
  { n<-length(x)
  k<-quantile(x[2:(n-1)],seq(0,1,length=nk+2))[2:(nk+1)]
  }
  if (length(k)!=nk) # right number of knots?
    stop(paste("there should be ",nk," supplied knots"))
  x <- x - x.shift # basis stabilizing shift
  k <- k - x.shift # knots treated the same!
  X<-matrix(0,length(x),object$bs.dim)
  for (i in 1:(m+1)) X[,i] <- x^(i-1)
  for (i in 1:nk) X[,i+m+1]<-(x-k[i])^m*as.numeric(x>k[i])
  object$X<-X # the finished model matrix
  if (!object$fixed) # create the penalty matrix
  { object$S[[1]]<-diag(c(rep(0,m+1),rep(1,nk)))
  }
  object$rank<-nk  # penalty rank
  object$null.space.dim <- m+1  # dim. of unpenalized space
  ## store "tr" specific stuff ...
  object$knots<-k;object$m<-m;object$x.shift <- x.shift
  object$df<-ncol(object$X)     # maximum DoF (if unconstrained)
  class(object)<-"tr.smooth"  # Give object a class
  object 
  }


#' prediction method function for the tr smooth class
#'
#' @param object smooth object for gam class
#' @param data the new data to predict on
#''
#' @return X the prediction matrix
Predict.matrix.tr.smooth<-function(object,data)
  ## prediction method function for the tr smooth class
{ requireNamespace("splines")
  x <- data[[object$term]]
  x <- x - object$x.shift # stabilizing shift
  # spline order (3=cubic)
  # knot locations
  # number of knots
  m <- object$m;
  k<-object$knots
  nk<-length(k)
  X<-matrix(0,length(x),object$bs.dim)
  for (i in 1:(m+1)) X[,i] <- x^(i-1)
  for (i in 1:nk) X[,i+m+1] <- (x-k[i])^m*as.numeric(x>k[i])
  X # return the prediction matrix
}