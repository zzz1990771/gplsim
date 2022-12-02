test_that("test the main functions (gplsim) with binary case", {
  library(gplsim)
  # parameter settings
  n=5000
  true.theta = c(1, 1, 1)/sqrt(3)
  # Generating data (binary case)
  # This function generates a sin bump model as in Yu et al. (2007).
  # You may use your data instead.
  data <- generate_data(n,true.theta=true.theta,family="binomial")
  y=data$Y       # binary response
  X=data$X       # single index term ;
  Z=data$Z       # partially linear term ;

  # Fit the generalised partially linear single-index models
  result <- gplsim(y,X,Z,family = binomial)

  # summary of the fitted model
  est_beta <- summary(result)$p.coeff.sim

  summary(result)

  expect_equal(as.vector(est_beta), true.theta, tolerance = 0.2)
})


test_that("test the main functions (gplsim) with gaussian case, user init", {
  library(gplsim)
  # parameter settings
  n=2000
  true.theta = c(1, 1, 1)/sqrt(3)
  # Gaussian case
  # This function generate a plain sin bump model with gaussian response.
  data <- generate_data(n,true.theta=true.theta,family="gaussian")
  y=data$Y       # continuous response
  X=data$X       # single index term ;
  Z=data$Z       # partially linear term ;

  # Fit the generalised partially linear single-index models
  result <- gplsim(y,X,Z,family = gaussian,user.init = c(1,1,1))

  # summary of the fitted model
  est_beta <- summary(result)$p.coeff.sim

  expect_equal(as.vector(est_beta), true.theta, tolerance = 0.2)
})

test_that("test the main functions (gplsim) with gaussian case", {
  library(gplsim)
  # parameter settings
  n=2000
  true.theta = c(1, 1, 1)/sqrt(3)
  # Gaussian case
  # This function generate a plain sin bump model with gaussian response.
  data <- generate_data(n,true.theta=true.theta,family="gaussian")
  y=data$Y       # continuous response
  X=data$X       # single index term ;
  Z=data$Z       # partially linear term ;

  # Fit the generalised partially linear single-index models
  result <- gplsim(y,X,Z,family = gaussian)

  # summary of the fitted model
  est_beta <- summary(result)$p.coeff.sim

  expect_equal(as.vector(est_beta), true.theta, tolerance = 0.2)
})

test_that("test the main functions (gplsim) with gaussian case", {
  library(gplsim)
  # parameter settings
  n=2000
  true.theta = c(1, 1, 1)/sqrt(3)
  # Gaussian case
  # This function generate a plain sin bump model with gaussian response.
  data <- generate_data(n,true.theta=true.theta,family="poisson")
  y=data$Y       # continuous response
  X=data$X       # single index term ;
  Z=data$Z       # partially linear term ;

  # Fit the generalised partially linear single-index models
  result <- gplsim(y,X,Z,family = poisson)

  # summary of the fitted model
  est_beta <- summary(result)$p.coeff.sim

  expect_equal(as.vector(est_beta), true.theta, tolerance = 0.2)
})


test_that("test the main functions (gplsim) with gaussian case, profile", {
  library(gplsim)
  # parameter settings
  n=2000
  true.theta = c(1, 1, 1)/sqrt(3)
  # Gaussian case
  # This function generate a plain sin bump model with gaussian response.
  set.seed(123456)
  data <- generate_data(n,true.theta=true.theta,family="gaussian")
  set.seed(NULL)
  y=data$Y       # continuous response
  X=data$X       # single index term ;
  Z=data$Z       # partially linear term ;

  # Fit the generalised partially linear single-index models
  result <- gplsim(y,X,Z,family = gaussian, profile = FALSE)

  # summary of the fitted model
  est_beta <- summary(result)$p.coeff.sim

  expect_equal(as.vector(est_beta), true.theta, tolerance = 0.2)
})

test_that("test the main functions (gplsim) with gaussian case, profile, missing Z", {
  library(gplsim)
  # parameter settings
  n=2000
  true.theta = c(1, 1, 1)/sqrt(3)
  # Gaussian case
  # This function generate a plain sin bump model with gaussian response.
  set.seed(123456)
  data <- generate_data(n,true.theta=true.theta,family="gaussian")
  set.seed(NULL)
  y=data$Y       # continuous response
  X=data$X       # single index term ;
  Z=data$Z       # partially linear term ;

  # Fit the generalised partially linear single-index models
  result <- gplsim(y,X,Z=NULL,family = gaussian, profile = FALSE)

  # summary of the fitted model
  est_beta <- summary(result)$p.coeff.sim

  expect_equal(as.vector(est_beta), true.theta, tolerance = 0.2)
})


test_that("test the main functions (gplsim) with gaussian case, missing Z", {
  library(gplsim)
  # parameter settings
  n=2000
  true.theta = c(1, 1, 1)/sqrt(3)
  # Gaussian case
  # This function generate a plain sin bump model with gaussian response.
  set.seed(123456)
  data <- generate_data(n,true.theta=true.theta,family="gaussian")
  set.seed(NULL)
  y=data$Y       # continuous response
  X=data$X       # single index term ;
  Z=data$Z       # partially linear term ;

  # Fit the generalised partially linear single-index models
  result <- gplsim(y,X,Z=NULL,family = gaussian)

  # summary of the fitted model
  est_beta <- summary(result)$p.coeff.sim

  expect_equal(as.vector(est_beta), true.theta, tolerance = 0.2)
})

test_that("test the main functions (gplsim) with gaussian case, with p<2", {
  library(gplsim)
  # parameter settings
  n=2000
  true.theta = c(1, 1, 1)/sqrt(3)
  # Gaussian case
  # This function generate a plain sin bump model with gaussian response.
  set.seed(123456)
  data <- generate_data(n,true.theta=true.theta,family="gaussian")
  set.seed(NULL)
  y=data$Y       # continuous response
  X=data$X       # single index term ;
  Z=data$Z       # partially linear term ;

  # Fit the generalised partially linear single-index models
  expect_error(result <- gplsim(y,X[1,],NULL,family = gaussian, profile = FALSE))
})

test_that("test the main functions (gplsim) with gaussian case, with p<2", {
  # parameter settings
  n=8000
  true.theta = c(1, 1, 1)/sqrt(3)
  # Gaussian case
  # This function generate a plain sin bump model with gaussian response.
  set.seed(123456)
  data <- generate_data(n,true.theta=true.theta,family="gamma")
  set.seed(NULL)
  y=data$Y       # continuous response
  X=data$X       # single index term ;
  Z=data$Z       # partially linear term ;

  # Fit the generalised partially linear single-index models
  result <- gplsim(y,X,Z,family = Gamma, profile = FALSE)

  # summary of the fitted model
  est_beta <- summary(result)$p.coeff.sim

  expect_equal(as.vector(est_beta), true.theta, tolerance = 0.2)
})


test_that("test the formula functions (gplsim) with gaussian case, user init", {
  library(gplsim)
  # parameter settings
  n=2000
  true.theta = c(1, 1, 1)/sqrt(3)
  # Gaussian case
  # This function generate a plain sin bump model with gaussian response.
  data <- generate_data(n,true.theta=true.theta,family="gaussian")
  #formula method
  y=data$Y       # continuous response
  X=data$X       # single index term ;
  Z=data$Z       # partially linear term ;
  colnames(X) = c("x1","x2","x3")
  data = cbind(y,X,Z)
  formula = formula(y~si(x1+x2+x3)+Z)
  result <- gplsim(formula=formula,data=data)

  # summary of the fitted model
  est_beta <- summary(result)$p.coeff.sim

  expect_equal(as.vector(est_beta), true.theta, tolerance = 0.2)
})

