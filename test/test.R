#######################################
############     gplsim    ############


#######################################
n=10000
m=1
true.alpha = c(1, 1, 1)/sqrt(3)


# Generating data
data <- generate_data(n,true.alpha=true.alpha)
y=data$Y       # binary response
X=data$X   # single index term ;
Z=data$Z       # partially linear term ;




require(mgcv)
require(quantreg)
require(splines)

gplsimPs <- function(Y=Y,X=X,Z=Z,family=binomial,penalty=TRUE,user.init=FALSE){
  p <- dim(X)[2]
  if(user.init){
    if(length(user.init)!=(p-1)){
      stop("user.init length must be p-1")
    }else{
      init.alpha <- user.init
    }
  }else{
    er_np <- optim(rep(0,p-1), si, y=Y,x=X,z=Z, fam=family, hessian=TRUE, fx=TRUE)
    init.alpha <- er_np$par
  }
  
  er <- optim(init.alpha,si,y=Y,x=X,z=Z,fam=family,hessian=TRUE,fx=!penalty)
  b <- si(er_np$par,y=y,X,Z, fam=family, opt=FALSE) 
  return(b)
}

result <- gplsimPs(y,X,Z,user.init=c(0,0))


#binomial
alpha0 <- c(0,0)
## get initial alpha, using no penalization...
er_np <- optim(alpha0, si, y=y,x=X,z=Z, fam=binomial, hessian=TRUE, fx=TRUE,k=5)
## now get alpha with smoothing parameter selection...

b <- si(er_np$par,y=y,X,Z, fam=binomial, opt=FALSE) ## best fit model
#b
b$theta
b$coefficients[1:2]



#normal
alpha0 <- c(1,1)
## get initial alpha, using no penalization...
er <- optim(alpha0, si, y=y2,x=X,z=Z, fam=gaussian, hessian=TRUE, fx=TRUE,k=5)
## now get alpha with smoothing parameter selection...
er <- optim(er$par,si,y=y2,x=X,z=Z,fam=gaussian,hessian=TRUE,k=10)
b <- si(er$par,y=y2,X,Z, fam=gaussian, opt=FALSE) ## best fit model
#b
b$theta
b$coefficients[1]




# with iteration;
para.old = rep(0, 3)
ptm <- proc.time()
while((delta>=epsilon)&step<100){
  step=step+1

  alpha.old=para.old[1:3]
  er <- optim(alpha.old, si, y=y2,x=X,z=Z,fam=gaussian,hessian=TRUE,fx=TRUE,k=10)
  er <- optim(er$par, si, y=y2,x=X,z=Z,fam=gaussian,hessian=TRUE,k=10)
  b <- si(er$par, y2, X, Z, fam=gaussian, opt=FALSE)    ## best fit model
  alpha=b$theta                 # /sqrt(sum((b$theta)^2)) ;
  gamma=b$coefficients[1]
  para=c(alpha, gamma)
  delta = max(abs(para-para.old))
  para.old = para
}
proc.time() - ptm
#######################


###  Modified ;
### use standardized coefficient from gam
coeff = gam(y~X+Z-1, family=gaussian)$coefficients
alpha0 <- coeff[1:3]*sign(coeff[1])/sqrt(sum(coeff[1:3]^2))
## get initial alpha, using no penalization...
er <- optim(alpha0, si2, y=y,x=X,z=Z,fam=gaussian, hessian=TRUE,fx=TRUE,k=10)
#er$par
## now get alpha with smoothing parameter selection...
er <- optim(er$par, si2, y=y, x=X,z=Z,fam=gaussian,hessian=TRUE,k=10)
b <- si2(er$par, y, X, Z, fam=gaussian, opt=FALSE) ## best fit model
b$theta
b$coefficients[1:2]


#binomial
coeff = gam(y~X+Z-1, family=binomial )$coefficients
alpha0 <- coeff[1:3]*sign(coeff[1])/sqrt(sum(coeff[1:3]^2))
## get initial alpha, using no penalization...
er <- optim(alpha0, si2, y=y,x=X,z=Z,fam=binomial , hessian=TRUE,fx=TRUE,k=10)
#er$par
## now get alpha with smoothing parameter selection...
er <- optim(er$par, si2, y=y, x=X,z=Z,fam=binomial ,hessian=TRUE,k=10)
b <- si2(er$par, y, X, Z, fam=binomial , opt=FALSE) ## best fit model
b$theta
b$coefficients[1:2]


# try iteration;
coeff = gam(y~X+Z-1, family=binomial)$coefficients
alpha0 <- coeff[1:3]*sign(coeff[1])/sqrt(sum(coeff[1:3]^2))
epsilon = 10^(-4)
delta =1
step=0

para.old = c(alpha0, coeff[4])
ptm <- proc.time()
while((delta>=epsilon)&step<100){
  step=step+1

  alpha.old=para.old[1:3]
  er <- optim(alpha.old, si2, y=y,x=X,z=Z,fam=gaussian, hessian=TRUE,fx=TRUE,k=10)
  er <- optim(er$par, si2, y=y,x=X,z=Z,fam=gaussian, hessian=TRUE,k=10)
  b <- si2(er$par, y,X,Z, fam=gaussian, opt=FALSE)    ## best fit model
  alpha=b$theta                 # /sqrt(sum((b$theta)^2)) ;
  gamma=b$coefficients[1]
  para=c(alpha, gamma)
  delta = max(abs(para-para.old))
  para.old = para
}
proc.time() - ptm
step
para





##########  case 1    ##########
tau0 = 0.5

nknots.2n205=NULL; min_bic.2n205=NULL; nd.2=NULL;
bic.2n205=NULL
an205 <-matrix(rep(0, 15), ncol=3)
bn205<-NULL
iter2.2n205 <- matrix(rep(0, 500), ncol=5);   step2n205<-NULL
fitn205 <- matrix(rep(0, 20), ncol=4)
para_estn205<-matrix(rep(0, 400), ncol=4)
epsilon=10^(-5);

ptm <- proc.time()
for (j in 1:5){

  # Generating data
  data <- generate_data(n,true.alpha=true.alpha)
  y=data$Y       # binary response
  X=data$X   # single index term ;
  Z=data$Z       # partially linear term ;

  for (i in 1:5){
    nd.2[i] = i

    para.old <- quantreg::rq(y~ X+Z-1, tau=tau0)$coefficients
    delta=1
    step=0
    while ((delta>=epsilon)&(step<1000)){
      step=step+1;
      alpha.old <-para.old[1:3]
      alpha.old=sign(alpha.old[1])*alpha.old/sqrt(crossprod(alpha.old))
      beta.old <-para.old[4]

      index = X%*%alpha.old
      mm=min(index)
      MM=max(index)
      seq=seq(mm, MM, by=(MM-mm)/nd.2[i])
      inter= seq[-(nd.2[i]+1)]
      kts=c(mm, mm, mm, inter, MM, MM, MM, MM)
      bs2 = splines::splineDesign(knots=kts, ord=4, x=index, derivs=rep(0, n), outer.ok = F)
      ys <- y - Z*beta.old
      theta2 <- quantreg::rq(ys ~ bs2-1, tau=tau0)$coefficients
      eta = bs2%*%theta2
      r = y - eta - Z*beta.old;
      bs.dev2 = splines::splineDesign(knots=kts, ord=4, x=index, derivs=rep(1, n), outer.ok=F)
      etada = (as.vector(bs.dev2%*%theta2))*X

      # etadb = -(bs2%*%solve(t(bs2)%*%bs2)%*%t(bs2))%*%Z
      # etad = cbind(etada, etadb)
      # f1 = etad;
      # f1[, 4] = f1[, 4] + Z;
      ##f1=f1*diag(std(f1).^(-1));
      #y_new = r + f1%*%para.old;
      #x_new = f1
      #para = rq(y_new ~ x_new-1, tau=tau0)$coefficients
      ys_simple = y - eta + etada%*%alpha.old
      xs_simple = cbind(etada, Z)
      para = rq(ys_simple ~ xs_simple - 1, tau=tau0)$coefficients

      alpha =sign(para[1])*para[1:3]/ifelse(sqrt(crossprod(para[1:3]))==0,0.1^8,sqrt(crossprod(para[1:3])))
      beta = para[4]
      para=c(alpha, beta)
      delta = max(abs(para-para.old))
      para.old = para
    }     # end of while loop

    step2n205[i]=step
    fitn205[i, ]=para
    an205[i, ]=fitn205[i, 1:3]
    bn205[i]=fitn205[i, 4]
    index = X%*%an205[i, ]
    mm=min(index)
    MM=max(index)
    seq=seq(mm, MM, by=(MM-mm)/nd.2[i])
    inter= seq[-(nd.2[i]+1)]
    kts=c(mm, mm, mm, inter, MM, MM, MM, MM)
    bs2 = splines::splineDesign(knots=kts, ord=4, x=index, derivs=rep(0, n), outer.ok = F)
    ys = y - Z*bn205[i]
    theta2 <- rq(ys~bs2-1, tau=tau0)$coefficients

    s = y - bs2%*%theta2 - Z*bn205[i]
    rho= tau0*s - s*(s<0)
    loss = mean(rho)
    bic.2n205[i] = log(loss)+log(n)*(nd.2[i]-1+4+1)/(2*n)
  }    # end of iteration over num of knots

  iter2.2n205[j, ] = step2n205
  nknots.2n205[j] = nd.2[which(bic.2n205==min(bic.2n205))]-1
  min_bic.2n205[j] = bic.2n205[which(bic.2n205==min(bic.2n205))]
  para_estn205[j, ] = fitn205[which(bic.2n205==min(bic.2n205)), ]
}
proc.time() - ptm

length(iter2.2t205[iter2.2t205>50])
summary(min_bic.2t205)
table(nknots.2t205)
apply(para_estt205, 2, mean)
apply(para_estt205, 2, var)
c(3, 2, 1)/sqrt(14)
loss = exp(min_bic.2n205 - log(n)*(nd.2[i]-1+4+1)/(2*n) )
