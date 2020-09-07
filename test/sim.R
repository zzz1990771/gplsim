library(doParallel)
library(foreach)
cl <- makeCluster(12)
registerDoParallel(cl)

# Gaussian case
# parameter settings
M=200
gaussian_coefs <- foreach(i=1:M, .combine=rbind,.packages = c("gplsim","splines","mgcv")) %dopar%{
  n=1000
  true.theta = c(1, 1, 1)/sqrt(3)
  # This function generate a plain sin bump model with gaussian response.
  data <- generate_data(n,true.theta=true.theta,family="Gaussian")
  y=data$Y       # continous response
  X=data$X       # single index term ;
  Z=data$Z       # partially linear term ;
  result <- gplsim(y,X,Z,user.init=NULL,family = gaussian)
  c(result$theta,result$gamma)
}


apply(gaussian_coefs,2,mean)
apply(gaussian_coefs,2,sd)



# Binomial case
# parameter settings
M=200
binomial_coefs <- foreach(i=1:M, .combine=rbind,.packages = c("gplsim","splines","mgcv")) %dopar%{
  n=1000
  true.theta = c(1, 1, 1)/sqrt(3)
  # This function generate a plain sin bump model with Binomial response.
  data <- generate_data(n,true.theta=true.theta,family="Binomial")
  y=data$Y       # Binomial response
  X=data$X       # single index term ;
  Z=data$Z       # partially linear term ;
  result <- gplsim(y,X,Z,user.init=true.theta,family = binomial)
  c(result$theta,result$gamma)
}



apply(binomial_coefs,2,mean)
apply(binomial_coefs,2,sd)



# poisson case
# parameter settings
M=200
poisson_coefs<- foreach(i=1:M, .combine=rbind,.packages = c("splines","mgcv","gplsim")) %dopar%{
  n=1000
  true.theta = c(1, 1, 1)/sqrt(3)
  # This function generate a plain sin bump model with gaussian response.
  data <- generate_data(n,true.theta=true.theta,family="Poisson")
  y=data$Y       # poisson response
  X=data$X       # single index term ;
  Z=data$Z       # partially linear term ;
  result <- gplsim(y,X,Z,user.init=NULL,family = poisson(link = "log"),k=13)
  c(result$theta,result$gamma)
}


apply(poisson_coefs,2,mean)
apply(poisson_coefs,2,sd)






stopCluster(cl)


# Gaussian case for add percentile vurce
# parameter settings
M=200
gaussian_models <- lsit()
n=500
true.theta = c(1, 1, 1)/sqrt(3)
# This function generate a plain sin bump model with gaussian response.

