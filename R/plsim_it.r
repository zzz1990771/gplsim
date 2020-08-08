#function to fit plsim with iterative procedure, through non-linear least square
#input x,y,degree,nknots,alpha0(if different from default OLSe),maxiter,z


# normalize alpha0
alpha0 = sign(alpha0[1])*alpha0/sqrt(sum(alpha0^2))
# reduce X dimensionality by this projection.
xnew = x*alpha0;    
# get predicted yhat  
[yhat, beta0,yhat1,lambda]  = srspl(xnew,y,degree,nknots,z); 

param0=[alpha0;beta0]; % coefficient of parameter associated with z is imbedded in beta;

n = length(y);
m = length(param0);
d = length(alpha0); 
d2 = length(z);

param = param0;
paramold = 2 * param;
mse = 1;
iter=1;
tol=1.0E-6;