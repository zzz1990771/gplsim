\name{air}
\alias{air}
\docType{data}
\title{
dataset from an environmental study.
}
\description{
This dataset contains four variables: The concentration of the air pollutant ozone, wind speed, temperature and radiation. All of them are daily measurements for 111 days. Usually the concentration of the air pollutant ozone serves as the response variable while the other three are predictors.
}
\usage{data("air")}
\format{
  A data frame with 111 observations on the following 4 variables.
  \describe{
    \item{\code{ozone}}{a numeric vector in cube root ppb}
    \item{\code{radiation}}{a numeric vector in langley}
    \item{\code{temperature}}{a numeric vector in degrees F}
    \item{\code{wind_speed}}{a numeric vector in mph}

  }
}
\examples{
data(air)
y=air$ozone               # response
X=as.matrix(air[,3:4])    # single index term ;
Z=air[,2]                 # partially linear term ;

result <- gplsim(y,X,Z=Z,family = gaussian,k=10)
result$theta
result$coefficients
summary(result)

# Or you can try different spline basis
result <- gplsim(y,X,Z=Z,family = gaussian,bs="tp",k=10)
result$theta
result$coefficients
summary(result)

}
\keyword{datasets}
