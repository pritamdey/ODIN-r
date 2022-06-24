#Functions for fitting and evaluating the model

#MATHEMATICAL FUNCTIONS
expit=function(x) ifelse(x>500,1,exp(x)/(1+exp(x)))

l2_norm=function(x) {sqrt(sum(x^2))}

logit=function(x)  {log(x/(1-x))}

quadratic_form=function(x,a) as.vector(t(x)%*%a%*%x)

#PREDICTION & LOG-LIKELIHOOD
logistic_predict=function(z,bet,x)
{
  expit(z+x%*%bet)
}

negative_log_like=function(z,bet,x,a)
{
  eta=z+x%*%bet
  -sum(a*eta-log(1+exp(eta)))/ncol(bet)
}

total_log_like=function(z,bet,x,a,lam)
{
  negative_log_like(z,bet,x,a)+(lam/2)*sum(z**2)
}

#ESTIMATION FUNCTION
mm_logistic=function(A,X,lam=0.001,tol=1e-7,maxiter=2000,init_Z=NA,init_Beta=NA,print_out=FALSE)
{
  N=dim(A)[2]    #Sample Size
  L=dim(A)[1]    #Length of the half vectors
  p=dim(X)[2]
  
  if(is.na(init_Z))
    Z=rep(0,L)
  else
    Z=init_Z
  if(is.na(init_Beta))
    Beta=matrix(0,p,N)
  else
    Beta=init_Beta
  
  
  #Iteration parameters
  diff=Inf
  iteration=0
  
  #Keep track of the likelihood
  llike=rep(0,maxiter+1)
  norms=rep(0,maxiter)
  llike[iteration+1]=total_log_like(Z,Beta,X,A,lam)
  
  XTXinv=solve(t(X)%*%X)
  IminusPX=diag(L)-(X%*%XTXinv%*%t(X))
  Gam=solve(4*lam*diag(L)+IminusPX)
  start = Sys.time()
  
  while(diff>tol && iteration<maxiter)
  {
    iteration=iteration+1
    
    #Stuff needed to compute gradient
    B=logistic_predict(Z,Beta,X)
    api=A-B
    
    delZ=-lam*Z+IminusPX%*%rowMeans(api)
    
    #MM
    deltaZ=4*Gam%*%delZ
    deltaBeta=4*N*XTXinv%*%t(X)%*%(api/N-as.numeric(deltaZ))
    Z=Z+as.numeric(deltaZ)
    Beta=Beta+deltaBeta
    
    #New Likelihood
    llike[iteration+1]=total_log_like(Z,Beta,X,A,lam)
    
    #Norm of gradient
    #Negative of Gradient
    dirZ=-lam*Z+rowMeans(api)
    dirBeta=t(X)%*%api/N
    #Norm of gradient
    norms[iteration]=l2_norm(c(l2_norm(dirZ),l2_norm(dirBeta)))
    diff=(llike[iteration]-llike[iteration+1])/llike[iteration]
    
    
    if(print_out)
    {
      message(paste("Iteration:",iteration,"Likelihood:",llike[iteration+1],"Norm:",norms[iteration],"Relative_Likelihood_Change:",diff))
    }
  }
  message(paste("While loop ended with",iteration,"iterations in",difftime(Sys.time(), start, units = "secs"),"seconds."))
  llike=llike[1:(iteration+1)]
  return(list(Z=Z,Beta=Beta,llike=llike))
}

# CALCULATION OF INFLUENCE
influence_measure=function(list1,A,X,lam)
{
  s=Sys.time()
  
  Zhat = list1[[1]]
  Betahat = list1[[2]]
  
  L=length(Zhat)
  N=dim(Betahat)[2]
  p=dim(Betahat)[1]
  
  pis = logistic_predict(Zhat,Betahat,X)
  beta_bar = rowMeans(Betahat)
  s = Betahat %*% t(Betahat) / N - beta_bar %*% t(beta_bar)
  s_inv = solve(s)
  
  T_ = diag(rowSums(pis * (1 - pis)) + N * lam)
  d_beta = rep(0,N)
  
  for(i in 1:N)
  {
    pi_i = pis[,i]
    W = pi_i * (1 - pi_i)
    B_i = W*X
    Q_i = solve(t(X)%*%B_i)
    T_ = T_ - (B_i%*%Q_i%*%t(B_i))
    d_beta[i] = quadratic_form(Betahat[,i]-beta_bar,solve(s))
  }
  
  
  Z_inf = solve(T_, A - pis - lam * Zhat)
  IM = sqrt(colSums(Z_inf^2)) * N/(N-1)
  list(IM=IM,Z_d_beta=d_beta)
}