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
  q = as.vector(fU + rho*Z)

  if(family=="gaussian"){
    ylist <- lapply(vector(mode = "list", length = ncopy),function(x){q + rnorm(length(q),0,sigma)})
  }else if(family=="binomial"){
    py = exp(q)/(1+exp(q))
    ylist <- lapply(vector(mode = "list", length = ncopy),function(x){rbinom(length(q), size=1, prob=py)})
  }else if(family=="poisson"){
    py = exp(q)
    ylist <- lapply(vector(mode = "list", length = ncopy),function(x){rpois(length(q),py)})
  }else if(family=="gamma"){
    ylist <- lapply(vector(mode = "list", length = ncopy),function(x){
      mu <- q +1
      alpha <- exp(q)
      rgamma(length(q), shape=alpha, scale=mu/alpha)
    })

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
#' @param smooth_selection see ?gplsim
#' @param fam see ?gplsim
#' @param bs see ?gplsim
#' @param fx see ?gplsim
#' @param scale see ?gplsim
#' @param ... includes optional arguments user can pass to \code{mgcv::gam} or \code{glm}, such as \code{k}, which is the dimension of the basis of the smooth term and \code{m}, which is the order of the penalty for the smooth term
#'
#' @return b fitted gam object
si <- function(alpha,y,x,z,opt=TRUE,smooth_selection,fam,bs="ps", fx=FALSE,scale=scale,...) {
  ## from ..., expose to environment
  #print((...))
  list2env(list(...),environment())

  ## Fit single index model using gam call. Return ML is opt==TRUE
  ## and fitted gam with theta added otherwise...

  #reparameterized theta
  theta <- sign(alpha[1])*alpha/sqrt(sum(alpha^2))

  a <- x%*%theta
  if(is.null(z)){
    b <- mgcv::gam(y~s(a,bs=bs,fx=fx,m=m,k=k),family=fam,method= smooth_selection,scale=scale )
  }else{
    b <- mgcv::gam(y~s(a,bs=bs,fx=fx,m=m,k=k)+z,family=fam,method= smooth_selection,scale=scale )
    b$gamma <- b$coefficients[c(2:(1+NCOL(z)))]
  }
  if (opt) return(b$deviance) else {
    b$theta <- theta  ## add theta
    #b$gamma <- b$coefficients[c(2:(1+NCOL(z)))]
    class(b) <- c("gplsim","gam","glm","lm")
    return(b)
  }

}


#' @name gplsim
#' @export
gplsim <- function( ...)
  UseMethod("gplsim")





#' Function to fit generalized partially linear single-index models via penalized splines
#'
#' This function employs penalized spline (P-spline) to estimate generalized
#' partially linear single index models, which extend the generalized linear
#' models to include nonlinear effect for some predictors.
#'
#' For formula, method, see ?gplsim.formula
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
#' @param profile profile is a logical variable that indicates whether the algorithm with profile likelihood or algorithm with NLS procedure should be used. The default algorithm is set to algorithm with profile likelihood.
#' @param bs bs is a character variable that specifies the spline basis in the estimation of unknown univariate function of single index. Default is P-splines.
#' @param user.init The user.init is a numeric vector of the same length as the dimensionality of single index predictors. The users can use this argument to pass in any appropriate user-defined initial single-index coefficients based on prior information or domain knowledge. The default value is NULL.
#' @param ... includes optional arguments user can pass to \code{mgcv::gam} or \code{glm}, such as \code{k}, which is the dimension of the basis of the smooth term and \code{m}, which is the order of the penalty for the smooth term. Others include:
#' \code{scale} The optional argument scale is a numeric indicator with a default value set to -1. Any negative value including -1 indicates that the scale of response distribution is unknown, thus need to be estimated. Another option is 0 signaling scale of 1 for Poisson and binomial distribution and unknown for others. Any positive value will be taken as the known scale parameter.
#' \code{smooth_selection} The optional argument smooth_selection is another character variable that specifies the criterion used in the selection of a smoothing parameter. The supported criteria include "GCV.Cp","GACV.Cp", "ML","P-ML", "P-REML" and "REML", while the default criterion is "GCV.Cp".
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
#' @rdname gplsim
#' @method gplsim default
#' @export
gplsim.default <- function(Y=Y,X=X,Z=Z,family=gaussian(),penalty=TRUE,
                   penalty_type = "L2",profile = TRUE, bs="ps",user.init= NULL,...){
  #validate inputs of X,Y,Z
  stopifnot(is.numeric(Y),is.vector(Y))
  stopifnot(is.numeric(X),is.matrix(X))
  stopifnot("SIM predictos must no less than two"= dim(X)[2] >=2)
  stopifnot(is.numeric(X),is.matrix(X))
  stopifnot(is.null(Z) || (is.numeric(Z)) )

  #deal with ...
  args_list <- list(...)
  if (is.null(args_list$m)) args_list$m <- m <- 2
  if (is.null(args_list$k)) args_list$k <- k <- 13
  if (is.null(args_list$scale)) args_list$scale <- scale <- -1
  if (is.null(args_list$smooth_selection)) args_list$smooth_selection <- smooth_selection <- "GCV.Cp"
  list2env(list(...),environment())

  #check other arguments
  if(!is.null(user.init)){
      stopifnot("user.init length must be the same with X" = length(user.init) ==NCOL(X) )
      init.alpha <- user.init
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
    non.profile.fit <- plsiml(x=X,y=Y ,z=Z,degree=m,nknots=k,maxiter = 10,alpha0 = init.alpha)
    a <- X%*%non.profile.fit$alpha
    if(is.null(Z)){
      b <- mgcv::gam(Y~s(a,bs=bs,fx=!penalty,m=m,k=k),family=family,method= smooth_selection,scale=scale)
    }else{
      b <- mgcv::gam(Y~s(a,bs=bs,fx=!penalty,m=m,k=k)+Z,family=family,method= smooth_selection,scale=scale)
      b$gamma <- b$coefficients[c(2:(1+NCOL(Z)))]
    }
    b$theta <- non.profile.fit$alpha  ## add theta

    class(b) <- c("gplsim","gam","glm","lm")

  }else{
    er <- suppressWarnings(optim(init.alpha,si,y=Y,x=X,z=Z,fam=family, smooth_selection=smooth_selection, hessian=TRUE,fx=!penalty,bs=bs,scale=scale,m=m,k=k))
    b <- si(er$par,y=Y,X,Z,fam=family,smooth_selection=smooth_selection, opt=FALSE,bs=bs,scale=scale,m=m,k=k)
  }
  b$Xnames <- if(is.null(colnames(X)))  paste0("X.",seq(1,NCOL(X),by = 1)) else colnames(X)
  b$Znames <- if(is.null(colnames(Z)) && !is.null(Z))  paste0("Z.",seq(1,NCOL(Z),by = 1)) else colnames(Z)

  return(b)
}

#utils::globalVariables(c("yscale"))


#' Formula interface for gplsim
#'
#' This function add formula interface to gplsim function
#'
#' @param formula A model formula;
#' @param data A data matrix containing the variables in the formula.
#' @param family A \code{family} object: a list of functions and expressions
#' for defining \code{link} and \code{variance} functions.
#' Families supported are \code{binomial}, \code{gaussian}.
#' The default family is \code{gaussian}.
#' @param penalty Whether use penalized splines or un-penalized splines to fit the model. The default is TRUE.
#' @param penalty_type The optional argument penalty_type is a character variable, which specifies the type of penalty used in the penalized splines estimation. The default penalty type is {L}_{2} penalty, while {L}_{1} is also supported.
#' @param profile profile is a logical variable that indicates whether the algorithm with profile likelihood or algorithm with NLS procedure should be used. The default algorithm is set to algorithm with profile likelihood.
#' @param bs bs is a character variable that specifies the spline basis in the estimation of unknown univariate function of single index. Default is P-splines.
#' @param user.init The user.init is a numeric vector of the same length as the dimensionality of single index predictors. The users can use this argument to pass in any appropriate user-defined initial single-index coefficients based on prior information or domain knowledge. The default value is NULL.
#' @param ... includes optional arguments user can pass to \code{mgcv::gam} or \code{glm}, such as \code{k}, which is the dimension of the basis of the smooth term and \code{m}, which is the order of the penalty for the smooth term. Others include:
#' \code{scale} The optional argument scale is a numeric indicator with a default value set to -1. Any negative value including -1 indicates that the scale of response distribution is unknown, thus need to be estimated. Another option is 0 signaling scale of 1 for Poisson and binomial distribution and unknown for others. Any positive value will be taken as the known scale parameter.
#' \code{smooth_selection} The optional argument smooth_selection is another character variable that specifies the criterion used in the selection of a smoothing parameter. The supported criteria include "GCV.Cp","GACV.Cp", "ML","P-ML", "P-REML" and "REML", while the default criterion is "GCV.Cp".
#'
#' @return theta  Estimation of Theta
#' @return coefficients  the coefficients of the fitted model. Parametric coefficients are first, followed by coefficients for each spline term in turn.
#' @return ... See GAM object
#' @rdname gplsim
#' @method gplsim formula
#' @export
gplsim.formula <- function(formula, data, family=gaussian(), penalty=TRUE,
                   penalty_type = "L2",profile = TRUE, bs="ps",user.init= NULL,...){
  #parse formula
  lhs <- if(length(formula) == 3) formula[[2]] else NULL
  lhsVars <- all.vars(lhs)
  Y <- data[,lhsVars]

  rhs <- paste(deparse(formula[[length(formula)]], width.cutoff = 500), collapse="")
  rhs <- gsub("[[:blank:]]", "", rhs)
  si_vars <- strsplit(gsub("si\\(([^\\)]+)\\).*","\\1", rhs),split = "\\+")
  si_vars <- unlist(si_vars)
  X <- data[,si_vars]

  pl_vars <- strsplit(gsub("si\\(([^\\)]+)\\)(.*)","\\2", rhs),split = "\\+")
  pl_vars <- unlist(pl_vars)[-1]
  Z <- data[,pl_vars]

  cl <- match.call()
  cl[[1]] <- quote(gplsim)
  obj <- gplsim(Y=Y,X=X,Z=Z, ...)
  obj$call <- cl
  obj

}
