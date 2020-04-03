\name{si}
\alias{si}

\title{
An internal function to optimization and fitting. Don't use it solely.
}
\description{
A function to reparameterize theta and fit gam. also can be used as optimization objective function.
}
\usage{
si(alpha, y, x, z, opt = TRUE, k = 10, fam, fx = FALSE)
}
\arguments{
  \item{alpha}{
reparameterized theta
}
  \item{Y}{Response varible, should be a vector.}
  \item{X}{Single index covariates.}
  \item{Z}{Partially linear covariates.}
  \item{opt}{Whether it is used to model fitting or optimization objective function.}
  \item{k}{The dimension of the basis of splines}
  \item{fam}{A \code{family} object: a list of functions and expressions for defining \code{link} and
\code{variance} functions. Families supported are \code{binomial}, \code{gaussian}. The default family
is \code{binomial}.}
  \item{fixed}{Whether use penalised splines or un-penalised splines to fit the model. The default is FALSE, which use penalised splines.}
}
\details{

}
\value{
}
\references{

}
\author{

}
\note{

}



\seealso{
\code{\link{gplsim}}
}
\examples{

}

