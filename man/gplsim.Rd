\name{gplsim}
\alias{gplsimPs}
\alias{gplsim}

\title{
Function to fit generalised partially linear single-index models via penalised splines
}
\description{
This function employs penalised spline (P-spline) to estimate generalised partially linear singleindexmodels,which extend the generalised linear models to include nonlinear effect for some predictors.
}
\usage{
gplsim(Y, X, Z, family=binomial, penalty=TRUE, user.init=FALSE)
}

\arguments{
  \item{Y}{Response varible, should be a vector.}
  \item{X}{Single index covariates.}
  \item{Z}{Partially linear covariates.}
  \item{family}{A \code{family} object: a list of functions and expressions for defining \code{link} and
\code{variance} functions. Families supported are \code{binomial}, \code{gaussian}. The default family
is \code{binomial}.}
  \item{penalty}{Whether use penalised splines or un-penalised splines to fit the model. The default is TRUE.}
  \item{user.init}{User provided initial value for reparameterized theta. Note that the length of initial values should be p-1, instead of p.}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
\item{theta}{Estimation of Theta}
\item{coefficients}{the coefficients of the fitted model. Parametric coefficients are first, followed by coefficients for each spline term in turn.}
\item{...}{see GAM object}
}
\references{
Yu, Yan, Chaojiang Wu, and Yuankun Zhang. (2017). Penalised Spline Estimation for Generalised Partially Linear Single-Index Models. \emph{Statistics and Computing}. https://dl.acm.org/doi/abs/10.1007/s11222-016-9639-0
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
n=10000
true.alpha = c(1, 1, 1)/sqrt(3)

# Generating data (binary case)
# This function generates a sin bump model as in Yu et al. (2007).
# You may use your data instead.
data <- generate_data(n,true.alpha=true.alpha,binary=TRUE)
y=data$Y       # binary response
X=data$X   # single index term ;
Z=data$Z       # partially linear term ;

# Fit the generalised partially linear single-index models
result <- gplsim(y,X,Z,user.init=c(0,0),family = binomial)

# Estimation of Theta
result$theta

# the coefficients of the fitted model. Parametric coefficients are first, followed by coefficients for each spline term in turn.
result$coefficients

# summary of the fitted model
summary(result)


# Gaussian case
# This function generate a plain sin bump model with gaussian response.
data <- generate_data(n,true.alpha=true.alpha,binary=FALSE)
y=data$Y       # binary response
X=data$X   # single index term ;
Z=data$Z       # partially linear term ;

result <- gplsim(y,X,Z,user.init=c(0,0),family = gaussian)
result$theta
result$coefficients
summary(result)


}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ Generalised linear model }
\keyword{ gam }
\keyword{ Penalised splines }
