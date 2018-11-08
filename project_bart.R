setwd('~/Desktop/project')
library(parallel)
library(BART)
library('plotrix')
##---------Data processing-------##
# normalize the input points in the range [0,1]
input.normalize<-function(X,xmax,xmin){
  for(i in 1:nrow(X)){
    X[i,] = (X[i,]-xmin)/(xmax-xmin)
  }
  return(X)
}

## training data with N observations 
filename = 'data.csv'
N = 40 # no. of samples
df = read.csv(filename,header = TRUE)
colnames(df) = c('X1','X2','X3','X4','X5','Y','S')
X = df[1:N,1:5] # inputs
Y = df[1:N,6]   # response (stress)
xmax = apply(X,2,max)
xmin = apply(X,2,min)
X = input.normalize(X,xmax,xmin) 

## test data 
Xpred = df[(N+1):50,1:5]
Xpred = as.matrix(input.normalize(Xpred,xmax,xmin))
Ypred_true = df[(N+1):50,6]

##------------Fit the BART model-------##
## Priors 
shat=sd(Y)
alpha=0.95
beta=2
nc=100

nu=3
q=0.80
k=1

## MCMC settings
npost=1000 # size of posterior sample
burn=5000 # size of burnin iterations 
thin=10 # keep every how many posterior draws 

m=200 # Number of trees 
set.seed(100)
fit.bart = mc.wbart(X,Y,sigest=shat,sigdf=nu,sigquant=q,
                  k=k,power=beta, base=alpha,ntree=m,numcut=nc,
                  ndpost=npost, nskip=burn, keepevery=thin, mc.cores = 2L)

## Model diagnostics
par(mfrow=c(2,1))
par(mar = c(3, 3, 3, 3), oma = c(2, 2, 0.5, 0.5))
plot(fit.bart$sigma[seq(burn,length(fit.bart$sigma),by=thin)],type='l',xlab="Iteration", ylab=expression(sigma))
acf(fit.bart$sigma[seq(burn,length(fit.bart$sigma),by=thin)],main='')

Ypred = predict(fit.bart, Xpred)

## Prediction error 
mean((Ypred_true - apply(Ypred,2,mean))^2)
[1] 34.31467
## MSE
mean((Y - fit.bart$yhat.train.mean)^2)
[1] 1.159105

## Coverage rate 
coverage.rate.median <- function(actural,upper,lower){
  count <- 0 
  for ( i in 1: length(actural)){
    if (actural[i]>= lower[i] && actural[i]<= upper[i]
    ) {count <- count + 1}
  }
  coverage.rate <- count/length(actural)
  return(coverage.rate)
}

median=apply(Ypred,2,median)
lower=apply(Ypred,2,function(x) quantile(x,probs=0.025))
upper=apply(Ypred,2,function(x) quantile(x,probs=0.975))
coverage.rate.median(Ypred_true,upper,lower)

##-------Sensitivity analysis-----------##

## Define a predict function used for the sensitivity function
predict.bart <- function(x.test){
  pred.mat <- predict(fit.bart, x.test)
  out <- apply(pred.mat,2,mean)
  return(out)
}

library(sensitivity)
## Samples of size N from Uniform distribution to test the SI
n=100
X1=data.frame(matrix(runif(n*5),ncol=5))
X2=data.frame(matrix(runif(n*5),ncol=5))
si.S=sobolEff(predict.bart,X1=X1,X2=X2,order=1,nboot=0)
si.TS=sobolEff(predict.bart,X1=X1,X2=X2,order=0,nboot=0)

print(si.S$S$original)
print(si.TS$S$original)
print(si.S$S)
print(si.TS$S)



## ------------- Remove X5, refit the BART model ----------- ##
X = df[1:N,1:4] # inputs
Y = df[1:N,6]   # response (stress)
xmax = apply(X,2,max)
xmin = apply(X,2,min)
X = input.normalize(X,xmax,xmin) 

## test data 
Xpred = df[(N+1):50,1:4]
Xpred = as.matrix(input.normalize(Xpred,xmax,xmin))
Ypred_true = df[(N+1):50,6]

##------------Fit the BART model-------##
## Priors 
shat=sd(Y)
alpha=0.95
beta=2
nc=100

nu=3
q=0.80
k=1

## MCMC settings
npost=1000 # size of posterior sample
burn=5000 # size of burnin iterations 
thin=20 # keep every how many posterior draws 

m=200 # Number of trees 
fit.bart.2 = mc.wbart(X,Y,sigest=shat,sigdf=nu,sigquant=q,
                    k=k,power=beta, base=alpha,ntree=m,numcut=nc,
                    ndpost=npost, nskip=burn, keepevery=thin, mc.cores = 2L)

## Model diagnostics
par(mfrow=c(1,2))
plot(fit.bart$sigma[seq(burn,length(fit.bart$sigma),by=thin)],type='l',xlab="Iteration", ylab=expression(sigma), 
     main="Posterior of BART (m=100)")
acf(fit.bart$sigma,main=expression(sigma))

Ypred.2 = predict(fit.bart.2, Xpred)
## Prediction error 
mean((Ypred_true - apply(Ypred.2,2,mean))^2)
[1] 29.74628
## MSE
mean((Y - fit.bart.2$yhat.train.mean)^2)
[1] 1.554737

## Coverage rate 
coverage.rate.median <- function(actural,upper,lower){
  count <- 0 
  for ( i in 1: length(actural)){
    if (actural[i]>= lower[i] && actural[i]<= upper[i]
    ) {count <- count + 1}
  }
  coverage.rate <- count/length(actural)
  return(coverage.rate)
}

median=apply(Ypred,2,median)
lower=apply(Ypred,2,function(x) quantile(x,probs=0.025))
upper=apply(Ypred,2,function(x) quantile(x,probs=0.975))
coverage.rate.median(Ypred_true,upper,lower)

##-------Sensitivity analysis-----------##

## Define a predict function used for the sensitivity function
predict.bart <- function(x.test){
  pred.mat <- predict(fit.bart.2, x.test)
  out <- apply(pred.mat,2,mean)
  return(out)
}

library(sensitivity)
## Samples of size N from Uniform distribution to test the SI
n=100
X1=data.frame(matrix(runif(n*4),ncol=4))
X2=data.frame(matrix(runif(n*4),ncol=4))
si.S.2=sobolEff(predict.bart,X1=X1,X2=X2,order=1,nboot=0)
si.TS.2=sobolEff(predict.bart,X1=X1,X2=X2,order=0,nboot=0)


print(si.S.2$S$original)
print(si.TS.2$S$original)
print(si.S.2$S)
print(si.TS.2$S)

##-----Plot the CI of SI------##

plotCI(1:5, si.TS$S$original, 1.96*si.TS$S$`std. error`, ylim=c(0,0.6),pch=20, ylab='Total SI', xaxt='n', xlab='')
axis(side=1, at=seq(1,5,by=1), labels=c('X1','X2','X3','X4','X5'), lty=1, col='black',tck=0.03)
abline(h=0, lty=2, col='blue')
plotCI(1:5, si.S$S$original, 1.96*si.S$S$`std. error`, ylim=c(-0.1,0.8), pch=20, ylab='First-order SI', xaxt='n', xlab='')
axis(side=1, at=seq(1,5,by=1), labels=c('X1','X2','X3','X4','X5'), lty=1, col='black',tck=0.03)
abline(h=0, lty=2, col='blue')


plotCI(1:4, si.TS.2$S$original, 1.96*si.TS.2$S$`std. error`, ylim=c(0,0.6), pch=20, ylab='Total SI', xaxt='n', xlab='')
axis(side=1, at=seq(1,4,by=1), labels=c('X1','X2','X3','X4'), lty=1, col='black',tck=0.03)
abline(h=0, lty=2, col='blue')
plotCI(1:4, si.S.2$S$original, 1.96*si.S.2$S$`std. error`,ylim=c(-0.2,0.6), pch=20, ylab='First-order SI', xaxt='n', xlab='')
axis(side=1, at=seq(1,4,by=1), labels=c('X1','X2','X3','X4'), lty=1, col='black',tck=0.03)
abline(h=0, lty=2, col='blue')


