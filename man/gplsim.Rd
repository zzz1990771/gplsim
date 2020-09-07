\name{gplsim}
\alias{gplsimPs}
\alias{gplsim}

\title{
Function to fit generalised partially linear single-index models via penalised splines
}
\description{
This function employs penalised spline (P-spline) to estimate generalised partially linear single index models, which extend the generalised linear models to include nonlinear effect for some predictors.
}
\usage{
gplsim(Y, X, Z, family = binomial, penalty = TRUE,
penalty_type = "L2", scale = -1,smoothing_select = "GCV.Cp",
profile = TRUE, bs="ps", user.init=NULL)
}

\arguments{
  \item{Y}{Response varible, should be a vector.}
  \item{X}{Single index covariates.}
  \item{Z}{Partially linear covariates.}
  \item{family}{A \code{family} object: a list of functions and expressions for defining \code{link} and
\code{variance} functions. Families supported are \code{binomial}, \code{gaussian}. The default family
is \code{binomial}.}
  \item{penalty}{Whether use penalised splines or un-penalised splines to fit the model. The default is TRUE.}
  \item{penalty_type}{The optional argument penalty_type is a character variable, which specifies the type of penalty used in the penalized splines estimation. The default penalty type is {L}_{2} penalty, while {L}_{1} is also supported.}
  \item{scale}{The optional argument scale is a numeric indicator with a default value set to -1. Any negative value including -1 indicates that the scale of response distribution is unknown, thus need to be estimated. Another option is 0 signaling scale of 1 for Poisson and binomial distribution and unknown for others. Any positive value will be taken as the known scale parameter.}
  \item{smoothing_select}{The optional argument smoothing_select is another character variable that specifies the criterion used in the selection of a smoothing parameter. The supported criteria include ``GCV.Cp'',``GACV.Cp'', ``ML'',``P-ML'', ``P-REML'' and ``REML'', while the default criterion is ``GCV.Cp''.}
  \item{profile}{ The optional argument profile is a logical variable that indicates whether the algorithm with profile likelihood or algorithm with NLS procedure should be used. The default algorithm is set to algorithm with profile likelihood.}
  \item{bs}{The optional argument bs is a character variable that specifies the spline basis in the estimation of unknown univariate function of single index. The default has been set to ``ps''(P-splines with B-spline basis) while other choices are ``tr'' (truncated power basis), ``tp'' (Thin plate regression splines) and others.}
  \item{user.init}{The last optional argument user.init is a numeric vector of the same length as the dimensionality of single index predictors. The users can use this argument to pass in any appropriate user-defined initial single-index coefficients based on prior information or domain knowledge. The default value is NULL, which instructs the function to estimate the initial single index coefﬁcients by a generalized linear model.}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
The function gplsim returns an object class of gplsim, which extends the gam object and glm object.
\item{theta}{Estimation of Theta}
\item{coefficients}{the coefficients of the fitted model. Parametric coefficients are first, followed by coefficients for each spline term in turn.}
\item{...}{see GAM object}
}
\references{
Yu, Yan, Chaojiang Wu, and Yuankun Zhang. (2017). Penalised Spline Estimation for Generalised Partially Linear Single-Index Models. \emph{Statistics and Computing}. https://dl.acm.org/doi/abs/10.1007/s11222-016-9639-0

Yu, Yan, and David Ruppert. (2002). Penalized Spline Estimation for Partially Linear Single-Index Models. \emph{Journal of the American Statistical Association} 97 (460): 1042-54.https://www.jstor.org/stable/3085829

}
\author{
Zu, Tianhai.
Yu, Yan.
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link{si}}
}
\examples{
# parameter settings
n=1000
true.theta = c(1, 1, 1)/sqrt(3)

# Generating data (binary case)
# This function generates a sin bump model as in Yu et al. (2007).
# You may use your data instead.
data <- generate_data(n,true.theta=true.theta,family="Binomial")
y=data$Y       # binary response
X=data$X       # single index term ;
Z=data$Z       # partially linear term ;

# Fit the generalised partially linear single-index models
result <- gplsim(y,X,Z,family = binomial)

# Estimation of Theta
result$theta

# the coefficients of the fitted model. Parametric coefficients are first, followed by coefficients for each spline term in turn.
result$coefficients

# summary of the fitted model
summary(result)


#plot the estimated single index function curve
plot.si(result)
#par(new=T)
#plot.si(result,index=Z,xaxt="n", yaxt="n",col="red")


# Gaussian case
# This function generate a plain sin bump model with gaussian response.
data <- generate_data(n,true.theta=true.theta,family="Gaussian")
y=data$Y       # continous response
X=data$X       # single index term ;
Z=data$Z       # partially linear term ;

result <- gplsim(y,X,Z,family = gaussian)
result$theta
result$coefficients
summary(result)


#plot the estimated single index function curve
plot.si(result)
#par(new=T)
#plot.si(result,index=Z,xaxt="n", yaxt="n",col="red")


# A real data example
data(air)
y=air$ozone               # response
X=as.matrix(air[,3:4])    # single index term ;
Z=air[,2]                 # partially linear term ;

result <- gplsim(y,X,Z=Z,family = gaussian,bs="tr",k=10)
result$theta
result$coefficients
summary(result)

# Or you can try different spline basis
result <- gplsim(y,X,Z=Z,family = gaussian,bs="tp",k=10)
result$theta
result$coefficients
summary(result)

# to know more about air data
?air

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ Generalised linear model }
\keyword{ gam }
\keyword{ Penalised splines }
