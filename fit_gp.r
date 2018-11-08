setwd('~/Documents/STAT-Project/')
source("dace.sim.r")
source("sampler.r")
set.seed(80)

filename = 'data.csv'
N = 40 # no. of samples

maxIter = 20000 #max. iterations of Gibbs/MCMC
burn = 10000    #no. of initial samples to throw away
se   = 10     #save samples after this many iterations

df = read.csv(filename,header = TRUE)
colnames(df) = c('X1','X2','X3','X4','X5','Y','S')
X = df[1:N,1:5] # inputs
Y = df[1:N,6]   # response (stress)

# nomalize the inputs to [0,1]
xmax = apply(X,2,max)
xmin = apply(X,2,min)
X = input.normalize(X,xmax,xmin) 
F = input.transform(X) #non-linear transformation

#F = X
F = cbind(rep(1,nrow(F)),F) 

# Choose prior parameters
# Do a linear fitting of Y=F*beta 
fitbeta = lm(Y~.,data=data.frame(Y=Y,F=F))
beta.mu = c(fitbeta$coefficients[1],fitbeta$coefficients[3:7])
beta.prec = 1/(0.5*5)^2
Ye = as.vector(F%*%beta.mu)
# precision (lambda) ~ Gamma(a,b)
lambda.a = 2   
lambda.b = 0.5
# correlation (rho) ~ Beta(a,b)
rho.a = c(5,2,5,2,2)
rho.b = c(1,5,1,2,2)

# Grid of input points
l1=list(m1=abs(outer(X[,1],X[,1],"-")))
l2=list(m2=abs(outer(X[,2],X[,2],"-")))
l3=list(m3=abs(outer(X[,3],X[,3],"-")))
l4=list(m4=abs(outer(X[,4],X[,4],"-")))
l5=list(m5=abs(outer(X[,5],X[,5],"-")))
l.design=list(l1=l1,l2=l2,l3=l3,l4=l4,l5=l5)    


##-------------Fit the data----------------------
pi = list(az=lambda.a,bz=lambda.b,rhoa=rho.a,rhob=rho.b,
          beta_mean=beta.mu,beta_prec=beta.prec,F=F) #prior parameters
init=list(rho0=rho.a/(rho.a+rho.b),lambda0=1/var(Ye),beta0=beta.mu)#initial guess     
nit= list(maxIter=maxIter, burn=burn, se=se)
mh = list(rr=0.07)
est = regress(Y,l.design,nit,pi,mh,init)


##---------------Predict-------------------------
Xpred = df[(N+1):50,1:5]
Xpred = input.normalize(Xpred,xmax,xmin)

F0 = input.transform(Xpred)
Xall = rbind(X,Xpred)

F0 = cbind(rep(1,nrow(F0)),F0)

l1=list(m1=abs(outer(Xall[,1],Xall[,1],"-")))
l2=list(m2=abs(outer(Xall[,2],Xall[,2],"-")))
l3=list(m3=abs(outer(Xall[,3],Xall[,3],"-")))
l4=list(m4=abs(outer(Xall[,4],Xall[,4],"-")))
l5=list(m5=abs(outer(Xall[,5],Xall[,5],"-")))
l.v=list(l1=l1,l2=l2,l3=l3,l4=l4,l5=l5) 


Ypred = df[(N+1):50,6]
pred = predict(l.v,est,F,F0,eps=1e-6)

pred.error = sum((pred$mean-Ypred)^2)/length(Ypred)
cat('\nPrediction error =',pred.error)


##---------------Plot results-------------------------
# predictions and 95% confidence interval
par(mfrow=c(3,2))
for(k in 1:5){
  ind = order(Xpred[,k])
  plot(Xpred[,k],Ypred,col='red',pch=20,xlab=sprintf('X%d',k),ylab='Y')
  for(i in 1:nrow(pred$Y)){
    lines(Xpred[ind,k],pred$Y[i,ind],col='gray',lwd=0.25)
  }
  lines(Xpred[ind,k],pred$mean[ind],col='blue',lwd=1)
  lines(Xpred[ind,k],pred$mean[ind]-1.96*pred$sd[ind],lwd=0.75,col="black")
  lines(Xpred[ind,k],pred$mean[ind]+1.96*pred$sd[ind],lwd=0.75,col="black")
  points(Xpred[,k],Ypred,col='red',pch=20)
}


# draws of Gibbs/MCMC sampler
# plot 1000 random draws (else the plot gets messy)
ind.draws = sample(length(est$lambdaz),1000)

# rho
# display the posterior mean
est.rhoz.mean = apply(est$rhoz,2,mean)
cat('\nestimated rho =',est.rhoz.mean)
# plot the draws
par(mfrow=c(3,2))
for(k in 1:5){
  plot(est$rhoz[ind.draws,k],type='l',xlab='draws',ylab=sprintf('rho%d',k))
}
# plot the auto-correlation of samples of rho_i
par(mfrow=c(3,2))
for(k in 1:5){
  acf(est$rhoz[,k],main=sprintf('rho%d',k))
}


# beta
# display the posterior mean
est.betaz.mean = apply(est$betaz,2,mean)
cat('\nestimated beta =',est.betaz.mean)
# plot the draws
par(mfrow=c(3,2))
for(k in 1:6){
  plot(est$betaz[ind.draws,k],type='l',xlab='draws',ylab=sprintf('beta%d',k))
}
# plot the auto-correlation of samples of beta_i
par(mfrow=c(3,2))
for(k in 1:6){
  acf(est$betaz[,k],main=sprintf('beta%d',k))
}

# lambda
# display the posterior mean
est.lambdaz.mean = mean(est$lambdaz)
cat('\nestimated lambda =',est.lambdaz.mean)
# plot the posterior samples and the autocorrelation of lambda
par(mfrow=c(2,1))
plot(est$lambdaz[ind.draws],type='l',xlab='draws',ylab='lambda')
acf(est$lambdaz,main='lambda')

# empirical posterior density of samples from Gibbs/MCMC
# rho
par(mfrow=c(3,2))
for(k in 1:5){
  boxplot(est$rhoz[,k],xlab=sprintf('rho%d',k))
}

# beta
par(mfrow=c(3,2))
for(k in 1:6){
  boxplot(est$betaz[,k],xlab=sprintf('beta%d',k))
}

# lambda
par(mfrow=c(1,1))
boxplot(est$lambdaz,xlab='lambda')


##---------------Sensitivity Analysis----------------------
library(sensitivity)

gp_yhat <- function(pred.loc){
  #pred.loc should be normalized
  colnames(pred.loc) = colnames(X)
  Xall = rbind(X,pred.loc)
  l1=list(m1=abs(outer(Xall[,1],Xall[,1],"-")))
  l2=list(m2=abs(outer(Xall[,2],Xall[,2],"-")))
  l3=list(m3=abs(outer(Xall[,3],Xall[,3],"-")))
  l4=list(m4=abs(outer(Xall[,4],Xall[,4],"-")))
  l5=list(m5=abs(outer(Xall[,5],Xall[,5],"-")))
  l.v=list(l1=l1,l2=l2,l3=l3,l4=l4,l5=l5) 
  
  F0 = input.transform(pred.loc)
  F0 = cbind(rep(1,nrow(F0)),F0)
  pred = predict(l.v,est,F,F0,eps=1e-6)
  return(pred$mean)
}

N=100
X1=data.frame(matrix(runif(N*5),ncol=5))
X2=data.frame(matrix(runif(N*5),ncol=5))
si.S =sobolEff(model=gp_yhat,X1=X1,X2=X2,order=1,nboot=0)
si.TS=sobolEff(model=gp_yhat,X1=X1,X2=X2,order=0,nboot=0)
#si.S2=sobolEff(model=gp_yhat,X1=X1,X2=X2,order=2,nboot=0)
cat('1-Way sensitivity')
print(si.S$S)
cat('Total sensitivity')
print(si.TS$S)

#------------------Contour Plot-------------------------
Ngrid = 10 #no. of points along each axis
x1 = seq(0,5,length=Ngrid)
x2 = seq(0,5,length=Ngrid)
x3 = seq(0,5,length=Ngrid)
x4 = mean(X[,4])
x5 = mean(X[,5])
Xc = cbind(expand.grid(x1,x2,x3),x4,x5)

Yc = gp_yhat(Xc)
plot(Yc)

# Correlation plot
col.vals = rev(topo.colors(Ngrid))
plot(1:Ngrid^2,Yc[1:Ngrid^2],pch=20,col=col.vals[1],xlim=c(1,Ngrid^3),
     xlab='index',ylab='Tensile Stress',ylim=c(0,max(Yc)))
for(i in 2:Ngrid){
  points(((i-1)*Ngrid^2+1):(i*Ngrid^2),Yc[((i-1)*Ngrid^2+1):(i*Ngrid^2)],
         pch=20,col=col.vals[i])
}

ii = 4
Ncon = (ii-1)*Ngrid^2
plot((Ncon+1):(Ncon+Ngrid),Yc[(Ncon+1):(Ncon+Ngrid)],pch=20,col=col.vals[1],xlim=c(Ncon+1,Ncon+Ngrid^2),
     xlab='index',ylab='Tensile Stress',ylim=c(0,max(Yc)))
for(i in 2:Ngrid){
  points((Ncon+1+(i-1)*Ngrid):(Ncon+i*Ngrid),Yc[(Ncon+1+(i-1)*Ngrid):(Ncon+i*Ngrid)],
       pch=20,col=col.vals[i])
}

# from plot, it's clear that X3 reduces Y and (X1,X2) increases Y
# so we are going to make a contour plot of (X1,X2) keeping the rest fixed
Ngrid = 30
x1 = seq(0,5,length=Ngrid)
x2 = seq(0,5,length=Ngrid)
x3 = 0.0001
Xcc = cbind(expand.grid(x1,x2),x3,x4,x5)
Ycc = gp_yhat(Xcc)

#reshape Ycc to sensible form
Ycc = matrix(Ycc,Ngrid,Ngrid)

#calculate Ylim
yt = 7.5*sqrt(x1*(xmax[1]-xmin[1])+xmin[1])
Ylim = c()
for(i in 1:Ngrid){
  Ylim = cbind(Ylim,yt)
}

#plot the response on X1,X2 plane (Red is the highest value)
image(x1,x2,Ycc,col=rev(heat.colors(50)),main='Ydiff')

#plot the difference (Red is the highest value)
Ydiff = Ycc-Ylim
image(x1,x2,Ydiff,col=rev(heat.colors(50)),main='Ydiff')

#plot where cracks occur
Ycrack = Ydiff>0
image(x1,x2,Ycrack,col=c(0,1),main='crack')


