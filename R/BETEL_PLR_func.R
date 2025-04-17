library(matrixStats)
library(invgamma)
library(LaplacesDemon)
library(momentfit)
library(gtools)



loglik.prior<-function(theta){
  sum(dnorm(theta,mean=mutheta,sd=sigtheta,log=TRUE))
}


#Set the prior parameters
mualp<-0;sigalp<-1000 
mutheta<-0;sigtheta<-1000

## Y: Outcome 
## X: Design matrix for h(x) in the OR (without intercept)
## D: Treatment variable
## X_ps: Densign matrix for the PS model (without intercept)
## Function to estimate the ATE in the PS_cov regression model via the BETEL. See
## Luo et al. (2023) for detailed description for this approach.
## Reference: Luo, Y., Graham, D. J., McCoy, E. J. (2023) Semiparametric Bayesian doubly robust causal estimation. Journal of Statistical Planning and Inference, 225, 171-187.


PS_cov_BETEL<-function(Y,X,X_ps,D){
  
  
  ps_est <- glm(D~X_ps,family = binomial(link = "logit"))$fitted.values
  
  OR_pred <- model.matrix( ~ D+X+ps_est)
  ### calculate log likelihood for prior 

  ### estimating equations
  gf<-function(theta,deX){
    gf1<-apply(t(deX),1, function(t) t*((Y-deX%*%theta)))
    return(gf1)
  }
 
  ##initial theta
  old.theta<-lm.fit(OR_pred,Y)$coef
  emp_mod<-momentModel(gf, OR_pred, old.theta, grad=NULL, vcov="iid")
  
  ## solve the empirical likelihood to update theta
  old.theta<-solveGel(update(emp_mod, gelType="ETEL"), theta0=old.theta)$theta
  
  ### MCMC initializations 
  ###starting values
 
  nburn <- 1000
  nits<- 5000 + nburn
  theta.samp<-matrix(0,nrow=nits,ncol=ncol(OR_pred))
  theta.samp[1,]<-old.theta
  
  prop.sig<-vcov(lm(Y~OR_pred-1))
 
  fit <- gelFit(emp_mod,gelType="ETEL")
  
  emp_log_lik_old<-sum(log(getImpProb(fit)$pt))
  
  log_prior_old<-loglik.prior(old.theta)
  
  for (i in 2:nits){
  
    ##symmetric proposal for theta
    para_p<-as.vector(rmvt(1,mu=old.theta,S=prop.sig))
  ## log MH prob
    log_p<-dmvt(c(old.theta),mu=para_p,S=prop.sig,log=TRUE)-
      dmvt(para_p,mu=old.theta,S=prop.sig,log=TRUE)
 ## log prior density
    log_prior_new<-loglik.prior(para_p)  
    
    new.lam<-evalGel(emp_mod,para_p,gelType="ETEL")@lambda
    ## log emp likelihood value
    emp_log_lik_new<-sum(log(getImpProb(evalGel(emp_mod,para_p,new.lam,gelType="ETEL"))$pt))
   
    jmp<-min(log_p+emp_log_lik_new+log_prior_new-emp_log_lik_old-log_prior_old,0)
    
    
    if(log(runif(1)) <= jmp){ ### accept the move
      emp_log_lik_old<-emp_log_lik_new
      log_prior_old<-log_prior_new
      theta.samp[i,]<-old.theta<-para_p
      
    } else{## reject the move
     
      theta.samp[i,]<-old.theta
    
    }

  }
  
  theta.samp_fin<-theta.samp[(nburn+1):nits,]
  post_ate<-theta.samp_fin[,2]
  return(list(post_betel_pscov_ate = post_ate))
} 


