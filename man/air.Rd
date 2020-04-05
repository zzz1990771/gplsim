\name{air}
\alias{air}
\docType{data}
\title{
dataset from an environmental study.
}
\description{
This dataset contains four variables: The concentration of the air pollutant ozone, wind speed, temperature and radiation. All of them are daily measurements for 111 days. Usually the concentration of the air pollutant ozone serves as the reponse variable while the other three are predictors.
}
\usage{data("air")}
\format{
  A data frame with 111 observations on the following 4 variables.
  \describe{
    \item{\code{ozone}}{a numeric vector in cube root ppb}
    \item{\code{wind_speed}}{a numeric vector in mph}
    \item{\code{temperature}}{a numeric vector in degrees F}
    \item{\code{radiation}}{a numeric vector in langleys}
  }
}
\details{

}
\source{
need to confirm this
}
\references{
%%  ~~ possibly secondary sources and usages ~~
}
\examples{
data(air)
y=air$ozone               # response
X=as.matrix(air[,2:3])    # single index term ;
Z=air[,4]                 # partially linear term ;

result <- gplsim(y,X,Z,user.init=c(0),family = gaussian)
result$theta
result$coefficients
summary(result)

# Or you can try without user defined initials
result <- gplsim(y,X,Z,family = gaussian)
result$theta
result$coefficients
summary(result)

}
\keyword{datasets}
