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
