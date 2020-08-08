data(air)
y=air$ozone               # response
X=as.matrix(air[,3:4])    # single index term ;
Z=air[,2]                 # partially linear term ;
#Z=rep(0,length(Z))

result <- gplsim(y,X,Z,user.init=c(0),k=20,family = gaussian)
result$theta
result$coefficients
summary(result)

sort_index <- order(X%*%result$theta)
plot((X%*%result$theta)[sort_index],y[sort_index])
lines((X%*%result$theta)[sort_index],result$fitted.values[sort_index])

plot((X%*%result$theta)[sort_index],result$fitted.values[sort_index])





# parameter settings
n=200
true.alpha = c(1, 1, 1)/sqrt(3)

# Gaussian case
# This function generate a plain sin bump model with gaussian response.
data <- generate_data(n,true.alpha=true.alpha,binary=FALSE)
y=data$Y       # continous response
X=data$X       # single index term ;
Z=data$Z       # partially linear term ;

result <- gplsim(y,X,Z,user.init=FALSE,family = gaussian)
result$theta
result$coefficients
summary(result)
#plot the estimated single index function curve
plot.si(result)
par(new=T)
plot.si(result,index=Z,xaxt="n", yaxt="n",col="red")

