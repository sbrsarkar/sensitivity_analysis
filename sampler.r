# This file contains essential functions used in the Project

library(tictoc) # to measure the execution time

# transform the data to the form that the simulator takes
input.transform<-function(X){
  # Ec = 57000*sqrt(X[,1])
  # Eback = 2.7*(X[,4]*X[,2]^2)/386.088
  # Ebed = 2.7*(X[,5]*X[,3]^2)/386.088
  # Pback = X[,4]
  # Pbed  = X[,5]
  Ec = sqrt(X[,1])
  Eback = (X[,4]*X[,2]^2)
  Ebed = (X[,5]*X[,3]^2)
  Pback = X[,4]
  Pbed  = X[,5]
  Xnew = cbind(Ec,Eback,Ebed,Pback,Pbed)
  colnames(Xnew) = colnames(X)
  return(Xnew)
}

# normalize the input points in the range [0,1]
input.normalize<-function(X,xmax,xmin){
  for(i in 1:nrow(X)){
    X[i,] = (X[i,]-xmin)/(xmax-xmin)
  }
  return(X)
}

# returns the log full conditional density
logp<-function(y,lambdaz,rhoz,l.dez,pi,ymean)
{
  az=pi$az; bz=pi$bz
  
  if(any(rhoz<=0) || any(rhoz>=1))
    logposterior=-Inf
  else
  {
    n=length(y)
    k=length(rhoz)
    In=diag(n)
    Rz=rhogeodacecormat(l.dez,rhoz)$R
    Ez=1/lambdaz*(Rz+In*1e-14)
    cholEz=chol(Ez)
    Ezinv=chol2inv(cholEz)
    logdetEz=2*sum(log(diag(cholEz)))
    
    logpz=-1/2*logdetEz-1/2*t(y-ymean)%*%Ezinv%*%(y-ymean)
    logplambdaz=(az-1)*log(lambdaz)-bz*lambdaz
    logprhoz=0
    for(i in 1:k) logprhoz=logprhoz+(pi$rhoa[i]-1)*log(rhoz[i])+(pi$rhob[i]-1)*log(1-rhoz[i])
    
    logposterior=logpz+logplambdaz+logprhoz
  }
  logposterior
}

# Gibbs/MCMC Sampler
regress<-function(y,l.dez,nit,pi,mh,init,adapt=TRUE)
{
  
  F = pi$F
  
  # determine the length and dimensions
  n=length(y)
  k=length(l.dez)
  p=ncol(F)
  
  # iteration parameters
  N = nit$maxIter
  burn = nit$burn
  se = nit$se  
  
  Ns = floor((N-burn)/se) + 1 #total no. of samples saved
  #if(Ns<2000) stop("Ns>=2000 please\n")
  
  draw.lambdaz=rep(NA,Ns)
  draw.rhoz=matrix(NA,nrow=Ns,ncol=k)
  draw.betaz=matrix(NA,nrow=Ns,ncol=p)
  
  rr=rep(mh$rr,k)
  
  # initial guesses
  lambdaz = init$lambda0
  rhoz    = init$rho0
  betaz   = init$beta0
  
  
  #accept/reject ratio trackers
  accept.rhoz=rep(0,k)
  
  lastadapt=0
  
  In=diag(n)
  Ip=diag(p)
  one=rep(1,n)
  it = 1
  
  cat("\n ---------Bayesian Gaussian Process Interpolation model----------")
  cat("\n Total ",N," samples from the posterior will be generated")
  cat("\n Every ",se," samples will be saved after the burning initial ",burn," samples")
  cat("\n ")
  cat("\n The stepwidth for uniform corr. param proposal distn is rr=",rr[1])
  cat("\n Prior params:  lambda ~ Gamma(",pi$az,",",pi$bz,")")
  cat("\n ----------------------------------------------------------------\n\n\n")
  
  tic()
  for(i in 1:N)
  {
    ymean = F%*%betaz
    for(j in 1:k)
    {
      rhoz.new=rhoz
      rhoz.new[j]=runif(1,rhoz[j]-rr[j],rhoz[j]+rr[j])
      a=min(0,logp(y,lambdaz,rhoz.new,l.dez,pi,ymean)
            -logp(y,lambdaz,rhoz,l.dez,pi,ymean))
      if(log(runif(1))<a)
      {
        rhoz[j] = rhoz.new[j]
        accept.rhoz[j]=accept.rhoz[j]+1
      }
    }
    
    
    # Draw the marginal precision of the data (Gibbs step)
    Rz=rhogeodacecormat(l.dez,rhoz)$R
    Rz=Rz+In*1e-14  #cheat
    cholRz=chol(Rz)
    Rzinv=chol2inv(cholRz)
    lambdaz=rgamma(1,pi$az+n/2,pi$bz+0.5*t(y-ymean)%*%Rzinv%*%(y-ymean))
    
    # Draw the trend coefficients
    beta_cov = solve(lambdaz*t(F)%*%Rzinv%*%F + pi$beta_prec*Ip)
    beta_mean= beta_cov%*%(lambdaz*t(F)%*%Rzinv%*%y + pi$beta_prec*pi$beta_mean)
    beta_cov_chol = chol(beta_cov)
    zu = rnorm(p)
    betaz = beta_mean + t(beta_cov_chol)%*%zu
    
    
    
    # Adaptive MCMC stuff:
    # adapt the proposal step every N/20 iterations, 
    # but only for the first 50% of the iterations
    if(adapt && i%%(N/20)==0 && i<(N*.5+1))
    {
      rate.rhoz=accept.rhoz/(i-lastadapt)
      cat("Adapting rates from ",rate.rhoz,"\n");
      for(j in 1:k)
        if(rate.rhoz[j]>.49 || rate.rhoz[j]<.39) 
            rr[j]=rr[j]*rate.rhoz[j]/.44
      lastadapt=i
      accept.rhoz=rep(0,k)
    }

    if ((i/N*100)%%10==0){
      cat("\n",i/N*100," percent complete\n")
    }
    
    # save the samples
    if((i>=burn) && (i-burn)%%se==0){
        draw.rhoz[it,]  = rhoz
        draw.betaz[it,] = betaz
        draw.lambdaz[it]= lambdaz
        #cat('i=',i,' it=',it)
        it = it + 1
    }

  }
  
  cat("\n Complete.")
  stop = toc(quiet=TRUE)
  te = stop$toc-stop$tic
  cat('Time elapsed: ',floor(te/3600),'hrs ',floor(te/60),'mins ',round(te%%60),'secs\n')
  rate.rhoz=accept.rhoz/(i-lastadapt)
  cat("[rate.rhoz=",rate.rhoz,"]\n")
  
  return(list(y=y,lambdaz=draw.lambdaz,rhoz=draw.rhoz,betaz=draw.betaz))
}

# Predict function
predict<-function(l.v,fit,F,F0,eps=1e-6)
{
  y = fit$y
  n=length(fit$y)  # num observations
  m=nrow(l.v[[1]][[1]])-n  # num prediction points
  N=length(fit$lambdaz)
  
  
  draw.preds=matrix(0,nrow=N,ncol=m)
  
  for(i in 1:N)
  {
    Rall=rhogeodacecormat(l.v,fit$rhoz[i,])$R
    
    # Extract the sub-matrices we need
    Ryy=Rall[1:n,1:n]
    Rgg=Rall[(n+1):(n+m),(n+1):(n+m)]
    Rgy=matrix(Rall[(n+1):(n+m),1:n],m,n)
    Ryy.inv=chol2inv(chol(Ryy))
    
    # Mean of conditional distribution:
    m.cond=F0%*%fit$betaz[i,] + Rgy%*%Ryy.inv%*%(y-F%*%fit$betaz[i,])
    
    # Covariance of conditional distribution:
    E.cond=fit$lambdaz[i]^(-1)*(Rgg-Rgy%*%Ryy.inv%*%t(Rgy))
    
    # Let's generate a realization!
    L=t(chol(E.cond+diag(m)*eps))
    u=rnorm(m)
    draw.preds[i,]=m.cond + L%*%u
  }
  
  draw.mean = apply(draw.preds,2,mean)
  draw.sd   = apply(draw.preds,2,sd)
  return(list(Y=draw.preds,mean=draw.mean,sd=draw.sd))
}
