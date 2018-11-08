# This script allows us to fit the data to a GP by dropping X5
setwd('~/Documents/STAT-Project/') #change it the appropriate path
set.seed(150)
source("dace.sim.r")
source("sampler.r")

#---------------------Parameters------------------------
filename = 'data.csv'
N = 40 # no. of samples

maxIter = 20000 #max. iterations of Gibbs/MCMC
burn = 10000    #no. of initial samples to throw away
se   = 10       #save samples after this many iterations

x.drop = 5      #don't change it
F.drop = c(3,5)
#-----------------------------------------------------
df = read.csv(filename,header = TRUE)
colnames(df) = c('X1','X2','X3','X4','X5','Y','S')
X = df[1:N,1:5] # inputs
Y = df[1:N,6]   # response (stress)

# nomalize the inputs to [0,1]
xmax = apply(X,2,max)
xmin = apply(X,2,min)
X = input.normalize(X,xmax,xmin)
F = input.transform(X) #non-linear transformation

# drop column
F = F[,-F.drop]
X = X[,-x.drop]

#F = X
F = cbind(rep(1,nrow(F)),F) 

# Choose prior parameters
# Do a linear fitting of Y=F*beta 
fitbeta = lm(Y~.,data=data.frame(Y=Y,F=F))
#beta.mu = c(fitbeta$coefficients[1],fitbeta$coefficients[3:7])
beta.mu = fitbeta$coefficients[-2]
beta.prec = 1/(0.5*5)^2
Ye = as.vector(F%*%beta.mu)
# precision (lambda) ~ Gamma(a,b)
lambda.a = 2   
lambda.b = 0.5
# correlation (rho) ~ Beta(a,b)
rho.a = c(5,2,5,2)
rho.b = c(1,5,1,2)

# Grid of input points
l1=list(m1=abs(outer(X[,1],X[,1],"-")))
l2=list(m2=abs(outer(X[,2],X[,2],"-")))
l3=list(m3=abs(outer(X[,3],X[,3],"-")))
l4=list(m4=abs(outer(X[,4],X[,4],"-")))
#l5=list(m5=abs(outer(X[,5],X[,5],"-")))
#l.design=list(l1=l1,l2=l2,l3=l3,l4=l4,l5=l5)
l.design=list(l1=l1,l2=l2,l3=l3,l4=l4)    


##-------------Fit the data----------------------
pi = list(az=lambda.a,bz=lambda.b,rhoa=rho.a,rhob=rho.b,
          beta_mean=beta.mu,beta_prec=beta.prec,F=F) #prior parameters
init=list(rho0=rho.a/(rho.a+rho.b),lambda0=1/var(Ye),beta0=beta.mu)#initial guess     
nit= list(maxIter=maxIter, burn=burn, se=se)
mh = list(rr=0.07)
est = regress(Y,l.design,nit,pi,mh,init)


##---------------Predict-------------------------
Xpred = df[(N+1):50,1:5]
#Xpred = df[50,1:5]
Xpred = input.normalize(Xpred,xmax,xmin)
F0 = input.transform(Xpred)
#drop columns
Xpred = Xpred[,-x.drop]
F0   = F0[,-F.drop]

Xall = rbind(X,Xpred)


F0 = cbind(rep(1,nrow(F0)),F0)

l1=list(m1=abs(outer(Xall[,1],Xall[,1],"-")))
l2=list(m2=abs(outer(Xall[,2],Xall[,2],"-")))
l3=list(m3=abs(outer(Xall[,3],Xall[,3],"-")))
l4=list(m4=abs(outer(Xall[,4],Xall[,4],"-")))
l.v=list(l1=l1,l2=l2,l3=l3,l4=l4) 


Ypred = df[(N+1):50,6]
pred = predict(l.v,est,F,F0,eps=1e-6)

pred.error = sum((pred$mean-Ypred)^2)/length(Ypred)
cat('\nPrediction error =',pred.error)


##---------------Plot results-------------------------
# predictions and 95% confidence interval
par(mfrow=c(2,2))
for(k in 1:4){
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
par(mfrow=c(2,2))
for(k in 1:4){
  plot(est$rhoz[ind.draws,k],type='l',xlab='draws',ylab=sprintf('rho%d',k))
}
# plot the auto-correlation of samples of rho_i
par(mfrow=c(2,2))
for(k in 1:ncol(est$rhoz)){
  acf(est$rhoz[,k],main=sprintf('rho%d',k))
}


# beta
# display the posterior mean
est.betaz.mean = apply(est$betaz,2,mean)
cat('\nestimated beta =',est.betaz.mean)
# plot the draws
par(mfrow=c(2,2))
for(k in 1:ncol(est$betaz)){
  plot(est$betaz[ind.draws,k],type='l',xlab='draws',ylab=sprintf('beta%d',k))
}
# plot the auto-correlation of samples of beta_i
par(mfrow=c(2,2))
for(k in 1:ncol(est$betaz)){
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
par(mfrow=c(2,2))
for(k in 1:ncol(est$rhoz)){
  boxplot(est$rhoz[,k],xlab=sprintf('rho%d',k))
}

# beta
par(mfrow=c(2,2))
for(k in 1:ncol(est$betaz)){
  boxplot(est$betaz[,k],xlab=sprintf('beta%d',k))
}

# lambda
par(mfrow=c(1,1))
boxplot(est$lambdaz,xlab='lambda')


# Sensitivity Analysis
library(sensitivity)

gp_yhat <- function(pred.loc){
  #pred.loc should be normalized
  colnames(pred.loc) = colnames(X)
  Xall = rbind(X,pred.loc)
  l1=list(m1=abs(outer(Xall[,1],Xall[,1],"-")))
  l2=list(m2=abs(outer(Xall[,2],Xall[,2],"-")))
  l3=list(m3=abs(outer(Xall[,3],Xall[,3],"-")))
  l4=list(m4=abs(outer(Xall[,4],Xall[,4],"-")))
  l.v=list(l1=l1,l2=l2,l3=l3,l4=l4)
  
  pred.loc = cbind(pred.loc,0)
  F0 = input.transform(pred.loc)
  F0 = F0[,-F.drop]
  F0 = cbind(rep(1,nrow(F0)),F0)
  pred = predict(l.v,est,F,F0,eps=1e-6)
  return(pred$mean)
}

set.seed(200)
N=20
X1=data.frame(matrix(runif(N*4),ncol=4))
X2=data.frame(matrix(runif(N*4),ncol=4))
si.S =sobolEff(model=gp_yhat,X1=X1,X2=X2,order=1,nboot=0)
si.TS=sobolEff(model=gp_yhat,X1=X1,X2=X2,order=0,nboot=0)
#si.S2=sobolEff(model=gp_yhat,X1=X1,X2=X2,order=2,nboot=0)
cat('1-Way sensitivity')
print(si.S$S)
cat('Total sensitivity')
print(si.TS$S)
