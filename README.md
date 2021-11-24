# gplsim

Generalized partially linear single-index models (GPLSIM) are important tools in nonparametric regression. They extend popular generalized linear models to allow flexible nonlinear dependence on some predictors while overcoming the “curse of dimensionality”.

In gplsim package, we provides functions that employ penalised spline (P-spline) to estimate generalised partially linear single index models, which extend the generalised linear models to include nonlinear effect for some predictors.

We are trying to post this package to CRAN for more general users.

Or you can install it simply with the following code in R:
```
devtools::install_github("zzz1990771/gplsim")
library(gplsim)
?gplsim
```

If you run into any dependency issues, you can install the dependent package manually:

```
install.packages("mgcv","minpack.lm")
```

Quick example:

```
# parameter settings
n=1000
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

# Estimation of Theta
result$theta

# The coefficients of the fitted model. 
# Parametric coefficients are first, followed by coefficients for each spline term in turn.
result$coefficients

# summary of the fitted model
summary(result)

#plot the estimated single index function curve
plot.si(result)
#par(new=T)
#plot.si(result,index=Z,xaxt="n", yaxt="n",col="red")


# Gaussian case
# This function generate a plain sin bump model with gaussian response.
data <- generate_data(n,true.theta=true.theta,family="gaussian")
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

result <- gplsim(y,X,Z=Z,family = gaussian,k=10)
result$theta
result$coefficients
summary(result)

# Or you can try different spline basis
result <- gplsim(y,X,Z=Z,family = gaussian,bs="tp",k=10)
result$theta
result$coefficients
summary(result)

```
