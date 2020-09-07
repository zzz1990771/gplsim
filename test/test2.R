data(air)
y=air$ozone               # response
X=as.matrix(air[,3:4])    # single index term ;
Z=air[,2]                 # partially linear term ;
#Z=rep(0,length(Z))

result <- gplsim(y,X,Z,user.init=NULL,k=20,family = gaussian)
result$theta
result$coefficients
summary(result)

sort_index <- order(X%*%true.theta)
lines((X%*%true.theta)[sort_index],data$single_index_values[sort_index])

lines((X%*%result$theta)[sort_index],result$fitted.values[sort_index])

plot((X%*%result$theta)[sort_index],result$fitted.values[sort_index])


result$model$a
attr(result$smooth[[1]],"offset") <-  result$coefficients[1]
result$smooth[[1]]@offset
mgcv::plot.gam(result,xlab = "single index",ylab = "mean",select=1)

result$offset <- rep(result$coefficients[1],n)
mgcv::plot.gam(result)
result$coefficients[1]
# parameter settings
n=200
true.alpha = c(1, 1, 1)/sqrt(3)

# Gaussian case
# This function generate a plain sin bump model with gaussian response.
data <- generate_data(n,true.alpha=true.alpha,family="Gaussian")
y=data$Y       # continous response
X=data$X       # single index term ;
Z=data$Z       # partially linear term ;

result <- gplsim(y,X,Z,user.init=FALSE,family = gaussian)
result$theta
result$coefficients
summary(result)


#plot the estimated single index function curve
plot.si(result)
#par(new=T)
#plot.si(result,index=Z,xaxt="n", yaxt="n",col="red")

