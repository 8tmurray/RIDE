### Functions pertaining to the comparators for the JADE Design paper
### Written by Thomas A Murray on 10/23/2018


######################
### mTPI Functions ###
######################
#Function to get next dose based on mTPI
#ntot - number treated at the current dose level
#ndlt - number of DLTs at the current dose level
#target - target DLT probability
#delta - half-width of the equivalance interval
get.mtpi.next.dose = function(ntot,ndlt,target=0.3,delta=0.05){

  # calculate unit probability mass (UMP) for each decision
  UPM = rep(NA,3)
  UPM[1] = pbeta(target-delta,ndlt+1,ntot-ndlt+1)/(target-delta)
  UPM[2] = diff(pbeta(target+c(-1,1)*delta,ndlt+1,ntot-ndlt+1))/(2*delta)
  UPM[3] = pbeta(target+delta,ndlt+1,ntot-ndlt+1,lower.tail=FALSE)/(1-target-delta)
  
  # determine next dose based on the largest UPM
  decision = c("Escalate","Retain","Deescalate")[which.max(UPM)]

  return(decision)   
}
#ntot = 15; sapply(0:ntot,function(ndlt) get.mtpi.next.dose(ntot,ndlt,delta=0.06))


simulate.mtpi.trial = function(dlt.probs,start.dose=1,cohort.size=3,ncohorts=12,min.esc.size=1,target=0.3,threshold=0.975,delta=0.05){
  ndose = length(dlt.probs); n = y = rep(0,ndose); 
  
  curr.dose = start.dose; overdose = rep(FALSE,ndose)
  while(sum(n) < ncohorts*cohort.size){
    #Generate new observations at the current dose
    ntot = n[curr.dose] = n[curr.dose] + cohort.size
    ndlt = y[curr.dose] = y[curr.dose] + rbinom(1,cohort.size,dlt.probs[curr.dose])
    
    #Assess safety and decision criterion
    if(pbeta(target,ndlt+1,ntot-ndlt+1,lower.tail=FALSE) > threshold & n[curr.dose] >= 3) overdose[curr.dose:ndose] = TRUE
    if(overdose[1]) break #Stop trial when all doses are overdoses

    decision = get.mtpi.next.dose(ntot,ndlt,target=target,delta=delta)
    if((decision == "Deescalate" | overdose[curr.dose]) & curr.dose > 1) curr.dose = curr.dose-1
    if(decision == "Escalate" & !overdose[min(curr.dose+1,ndose)] & curr.dose < ndose & n[curr.dose] >= min.esc.size) curr.dose = curr.dose+1
  }
  dat = rbind(n,y,overdose); rownames(dat) = c("ntot","ndlt","overdose"); colnames(dat) = 1:ndose
  return(dat)  
}
#simulate.mtpi.trial(get.dlt.probs(6)$dlt.probs)



#######################
### Keyboard Design ###
#######################
#Function to get next dose based on keyboard
#ntot - number treated at the current dose level
#ndlt - number of DLTs at the current dose level
#target - target DLT probability
#delta - half-width of the equivalance interval
get.keyboard.next.dose = function(ntot,ndlt,target=0.3,delta=0.05){
  
  # determine key intervals based on target and delta
  getkey <- function(a1, a2)
  {
    delta=a2-a1
    lkey=NULL; rkey=NULL
    i=0; cutoff=0.3
    while(cutoff>0)
    {
      i=i+1
      cutoff = a1-i*delta
      lkey = c(cutoff, lkey)
    }
    lkey[lkey<0]=0
    
    i=0; cutoff=0.3
    while(cutoff<1)
    {
      i=i+1
      cutoff = a2+i*delta
      rkey = c(rkey, cutoff)
    }
    rkey[rkey>1]=1
    key=c(lkey, a1, a2, rkey)
    
    return(key)
  }
  keys = getkey(target-delta, target+delta)
  nkeys = length(keys)-1
  targetkey = which(keys==(target-delta))
  
  # calculate unit probability mass (UMP) for each key
  highkey = which.max((pbeta(keys[-1], ndlt+1, ntot-ndlt+1) - pbeta(keys[-length(keys)], ndlt+1, ntot-ndlt+1))/diff(keys))
  if(highkey<targetkey) decision = "Escalate"
  if(highkey==targetkey)decision = "Retain"
  if(highkey>targetkey) decision = "Deescalate"
  
  return(decision)   
}
#ntot = 15; sapply(0:ntot,function(ndlt) get.keyboard.next.dose(ntot,ndlt,delta=0.06))


simulate.keyboard.trial = function(dlt.probs,start.dose=1,cohort.size=3,ncohorts=12,min.esc.size=1,target=0.3,threshold=0.975,delta=0.05){
  ndose = length(dlt.probs); n = y = rep(0,ndose); 
  
  curr.dose = start.dose; overdose = rep(FALSE,ndose)
  while(sum(n) < ncohorts*cohort.size){
    #Generate new observations at the current dose
    ntot = n[curr.dose] = n[curr.dose] + cohort.size
    ndlt = y[curr.dose] = y[curr.dose] + rbinom(1,cohort.size,dlt.probs[curr.dose])
    
    #Assess safety and decision criterion
    if(pbeta(target,ndlt+1,ntot-ndlt+1,lower.tail=FALSE) > threshold & n[curr.dose] >= 3) overdose[curr.dose:ndose] = TRUE
    if(overdose[1]) break #Stop trial when all doses are overdoses
    
    decision = get.keyboard.next.dose(ntot,ndlt,target=target,delta=delta)
    if((decision == "Deescalate" | overdose[curr.dose]) & curr.dose > 1) curr.dose = curr.dose-1
    if(decision == "Escalate" & !overdose[min(curr.dose+1,ndose)] & curr.dose < ndose & n[curr.dose] >= min.esc.size) curr.dose = curr.dose+1
  }
  dat = rbind(n,y,overdose); rownames(dat) = c("ntot","ndlt","overdose"); colnames(dat) = 1:ndose
  return(dat)  
}
#simulate.keyboard.trial(c(0,0,0.1,0.2,0.4,0.8))





#####################
### CRM Functions ###
#####################
# Function to automatically deterimines a skeleton for the original CRM based on the method of Lee and Cheung (2009)
# ***adapted from the getprior() function in the dfcrm package
# ndose - number of dose levels
# target - target DLT toxicity probablity
# prior.best - prior guess for which dose is best
# delta - halfwidth of the indifference interval
get.skeleton <- function(ndose, target=0.3, delta=0.25*target, prior.best=ceiling((ndose+1)/2)){
  skeleton <- rep(NA, ndose)
  b <- rep(NA, ndose + 1); b[1] <- -Inf; b[(ndose + 1)] <- Inf
  skeleton[prior.best] <- target
  if(prior.best > 1){
    for(j in (prior.best-1):1){
      skeleton[j] <- exp(log(target-delta)*log(skeleton[j+1])/log(target+delta))
    }
  }
  if(prior.best < ndose){
    for(j in (prior.best+1):ndose){
      skeleton[j] <- exp(log(target+delta)*log(skeleton[j-1])/log(target-delta))
    }
  }
  return(skeleton)
}
#get.skeleton(ndose=6)

#Function to get next dose based on CRM with default skeleton
#curr.dose - current dose level
#n - number treated at each dose level
#y - number of DLTs at each dose level
#skeleton - CRM skeleton, i.e., prior DLT probabilities
#target - target DLT probability
get.crm.next.dose = function(curr.dose,n,y,skeleton,target=0.3){
  # posterior = likelihood x prior
  posterior <- function(alpha, p, y, n){
    sigma2 = 1.24^2; lik=1;
    for(j in 1:length(p)){
      pj = p[j]^(exp(alpha));
      lik = lik*pj^y[j]*(1-pj)^(n[j]-y[j]);
    }
    return(lik*exp(-0.5*alpha*alpha/sigma2));
  }
  
  # the posterior mean of ptox
  posttoxf <- function(alpha, p, y, n, j){ 
    p[j]^(exp(alpha))*posterior(alpha, p, y, n); 
  }
  
  
  marginal = integrate(posterior,lower=-Inf,upper=Inf, skeleton, y, n)$value
  ndose = length(skeleton); #number of doses in the trial
  pi.hat = numeric(ndose); # estimate of toxicity prob
  for(j in 1:ndose){
    pi.hat[j] = integrate(posttoxf,lower=-Inf,upper=Inf, skeleton, y, n, j)$value/marginal
  }
  
  # deterimine next dose
  bestdose = which.min(abs(pi.hat-target))
  if(bestdose>curr.dose) decision = "Escalate"
  if(bestdose==curr.dose) decision = "Retain"
  if(bestdose<curr.dose) decision = "Deescalate"

  return(decision)   
}
#skeleton = get.skeleton(6)
#ntot = 3; sapply(0:ntot,function(ndlt) get.crm.next.dose(curr.dose=1,n=c(ntot,rep(0,5)),y=c(ndlt,rep(0,5)),skeleton=skeleton))
#ntot = 3; sapply(0:ntot,function(ndlt) get.crm.next.dose(curr.dose=2,n=c(3,ntot,rep(0,4)),y=c(1,ndlt,rep(0,4)),skeleton=skeleton))


simulate.crm.trial = function(dlt.probs,skeleton,start.dose=1,cohort.size=3,ncohorts=12,min.esc.size=1,target=0.3,threshold=0.975){
  ndose = length(dlt.probs); n = y = rep(0,ndose); 
  
  curr.dose = start.dose; overdose = rep(FALSE,ndose)
  while(sum(n) < ncohorts*cohort.size){
    #Generate new observations at the current dose
    ntot = n[curr.dose] = n[curr.dose] + cohort.size
    ndlt = y[curr.dose] = y[curr.dose] + rbinom(1,cohort.size,dlt.probs[curr.dose])
    
    #Assess decision and safety criterion
    if(pbeta(target,ndlt+1,ntot-ndlt+1,lower.tail=FALSE) > threshold & n[curr.dose] >= 3) overdose[curr.dose:ndose] = TRUE
    if(overdose[1]) break #Stop trial when all doses are overdoses
    
    decision = get.crm.next.dose(curr.dose,n,y,skeleton,target)
    if((decision == "Deescalate" | overdose[curr.dose]) & curr.dose > 1) curr.dose = curr.dose-1
    if(decision == "Escalate" & !overdose[min(curr.dose+1,ndose)] & curr.dose < ndose & n[curr.dose] >= min.esc.size) curr.dose = curr.dose+1
  }
  dat = rbind(n,y,overdose); rownames(dat) = c("ntot","ndlt","overdose"); colnames(dat) = 1:ndose
  return(dat)  
}
#simulate.crm.trial(get.dlt.probs(6)$dlt.probs,skeleton)

#Function to get best dose based on CRM with default skeleton
#n - number treated at each dose level
#y - number of DLTs at each dose level
#skeleton - CRM skeleton, i.e., prior DLT probabilities
#target - target DLT probability
get.crm.best.dose = function(dat,skeleton,target=0.3){
  #Grab data
  n = dat[1,]; y = dat[2,]; 
  
  # posterior = likelihood x prior
  posterior <- function(alpha, p, y, n){
    sigma2 = 1.24^2; lik=1;
    for(j in 1:length(p)){
      pj = p[j]^(exp(alpha));
      lik = lik*pj^y[j]*(1-pj)^(n[j]-y[j]);
    }
    return(lik*exp(-0.5*alpha*alpha/sigma2));
  }
  
  # the posterior mean of ptox
  posttoxf <- function(alpha, p, y, n, j){ 
    p[j]^(exp(alpha))*posterior(alpha, p, y, n); 
  }
  
  marginal = integrate(posterior,lower=-Inf,upper=Inf, skeleton, y, n)$value
  ndose = length(skeleton); #number of doses in the trial
  pi.hat = numeric(ndose); # estimate of toxicity prob
  for(j in 1:ndose){
    pi.hat[j] = integrate(posttoxf,lower=-Inf,upper=Inf, skeleton, y, n, j)$value/marginal
  }
  
  #Determine Best Dose
  overdose = dat[3,]; #underdose = dat[4,]
  if(overdose[1]==1) best.dose = 0 #All overdoses
  if(overdose[1]==0) best.dose = which.min(abs(pi.hat[overdose==0]-target)) #Closest to the target among safe doses
  
  return(best.dose)   
}


###BMA-CRM

#Function to get next dose based on BMA-CRM with default skeletons
#curr.dose - current dose level
#n - number treated at each dose level
#y - number of DLTs at each dose level
#skeletons - matrix of BMA-CRM skeletons, i.e., prior DLT probabilities
#target - target DLT probability
#skeletons = t(sapply(c(1,4,6),function(j) get.skeleton(ndose=6,prior.best=j)))
get.bmacrm.next.dose = function(curr.dose,n,y,skeletons,target=0.3){
  # posterior = likelihood x prior
  posterior <- function(alpha, p, y, n){
    sigma2 = 1.24^2; lik=1;
    for(j in 1:length(p)){
      pj = p[j]^(exp(alpha));
      lik = lik*pj^y[j]*(1-pj)^(n[j]-y[j]);
    }
    return(lik*exp(-0.5*alpha*alpha/sigma2));
  }
  
  # the posterior mean of ptox
  posttoxf <- function(alpha, p, y, n, j){ 
    p[j]^(exp(alpha))*posterior(alpha, p, y, n); 
  }
  
  # calculate BMA pi.hat
  marginals = apply(skeletons,1,function(sk) integrate(posterior,lower=-Inf,upper=Inf, sk, y, n)$value)
  pi.hat = sapply(1:ncol(skeletons),function(j) sum(marginals/sum(marginals)*t(apply(skeletons,1,function(sk) integrate(posttoxf,lower=-Inf,upper=Inf, sk, y, n, j)$value))/marginals))
  
  # deterimine next dose
  bestdose = which.min(abs(pi.hat-target))
  if(bestdose>curr.dose) decision = "Escalate"
  if(bestdose==curr.dose) decision = "Retain"
  if(bestdose<curr.dose) decision = "Deescalate"
  
  return(decision)   
}


simulate.bmacrm.trial = function(dlt.probs,skeletons,start.dose=1,cohort.size=3,ncohorts=12,min.esc.size=1,target=0.3,threshold=0.975){
  ndose = length(dlt.probs); n = y = rep(0,ndose); 
  
  curr.dose = start.dose; overdose = rep(FALSE,ndose)
  while(sum(n) < ncohorts*cohort.size){
    #Generate new observations at the current dose
    ntot = n[curr.dose] = n[curr.dose] + cohort.size
    ndlt = y[curr.dose] = y[curr.dose] + rbinom(1,cohort.size,dlt.probs[curr.dose])
    
    #Assess decision and safety criterion
    if(pbeta(target,ndlt+1,ntot-ndlt+1,lower.tail=FALSE) > threshold & n[curr.dose] >= 3) overdose[curr.dose:ndose] = TRUE
    if(overdose[1]) break #Stop trial when all doses are overdoses
    
    decision = get.bmacrm.next.dose(curr.dose,n,y,skeletons,target)
    if((decision == "Deescalate" | overdose[curr.dose]) & curr.dose > 1) curr.dose = curr.dose-1
    if(decision == "Escalate" & !overdose[min(curr.dose+1,ndose)] & curr.dose < ndose & n[curr.dose] >= min.esc.size) curr.dose = curr.dose+1
  }
  dat = rbind(n,y,overdose); rownames(dat) = c("ntot","ndlt","overdose"); colnames(dat) = 1:ndose
  return(dat)  
}
#simulate.bmacrm.trial(c(0,0,0.1,0.2,0.5,0.8),skeletons)


#Function to get best dose based on BMA-CRM with default skeletons
#n - number treated at each dose level
#y - number of DLTs at each dose level
#skeletons - matrix of BMA-CRM skeleton, i.e., prior DLT probabilities
#target - target DLT probability
get.bmacrm.best.dose = function(dat,skeletons,target=0.3){
  #Grab data
  n = dat[1,]; y = dat[2,]; 
  
  # posterior = likelihood x prior
  posterior <- function(alpha, p, y, n){
    sigma2 = 1.24^2; lik=1;
    for(j in 1:length(p)){
      pj = p[j]^(exp(alpha));
      lik = lik*pj^y[j]*(1-pj)^(n[j]-y[j]);
    }
    return(lik*exp(-0.5*alpha*alpha/sigma2));
  }
  
  # the posterior mean of ptox
  posttoxf <- function(alpha, p, y, n, j){ 
    p[j]^(exp(alpha))*posterior(alpha, p, y, n); 
  }
  
  # calculate BMA pi.hat
  marginals = apply(skeletons,1,function(sk) integrate(posterior,lower=-Inf,upper=Inf, sk, y, n)$value)
  pi.hat = sapply(1:ncol(skeletons),function(j) sum(marginals/sum(marginals)*t(apply(skeletons,1,function(sk) integrate(posttoxf,lower=-Inf,upper=Inf, sk, y, n, j)$value))/marginals))
  
  #Determine Best Dose
  overdose = dat[3,]; #underdose = dat[4,]
  if(overdose[1]==1) best.dose = 0 #All overdoses
  if(overdose[1]==0) best.dose = which.min(abs(pi.hat[overdose==0]-target)) #Closest to the target among safe doses
  
  return(best.dose)   
}
#get.bmacrm.best.dose(dat,skeletons)



######################
### BOIN Functions ###
######################
library(BOIN)
#int - BOIN dose escalation interval
get.boin.next.dose = function(ntot,ndlt,boin.int){

  decision = "Retain"
  if(ndlt/ntot<boin.int[1]) decision = "Escalate"
  if(ndlt/ntot>boin.int[2]) decision = "Deescalate"
  
  return(decision)   
}
#boin.int = as.numeric(get.boundary(target=0.3,ncohort=12,cohortsize=3)[1:2]) #Dose escalation interval
#ntot = 6; sapply(0:ntot,function(ndlt) get.boin.next.dose(ntot,ndlt,boin.int))


simulate.boin.trial = function(dlt.probs,boin.int,start.dose=1,cohort.size=3,ncohorts=12,min.esc.size=1,target=0.3,threshold=0.975){
  ndose = length(dlt.probs); n = y = rep(0,ndose); 
  
  curr.dose = start.dose; overdose = rep(FALSE,ndose)
  while(sum(n) < ncohorts*cohort.size){
    #Generate new observations at the currdent dose
    ntot = n[curr.dose] = n[curr.dose] + cohort.size
    ndlt = y[curr.dose] = y[curr.dose] + rbinom(1,cohort.size,dlt.probs[curr.dose])
    
    #Assess safety and decision criterion
    if(pbeta(target,ndlt+1,ntot-ndlt+1,lower.tail=FALSE) > threshold & n[curr.dose] >= 3) overdose[curr.dose:ndose] = TRUE
    if(overdose[1]) break #Stop trial when all doses are overdoses
    
    decision = get.boin.next.dose(ntot,ndlt,boin.int)
    if((decision == "Deescalate" | overdose[curr.dose]) & curr.dose > 1) curr.dose = curr.dose-1
    if(decision == "Escalate" & !overdose[min(curr.dose+1,ndose)] & curr.dose < ndose & n[curr.dose] >= min.esc.size) curr.dose = curr.dose+1
  }
  dat = rbind(n,y,overdose); rownames(dat) = c("ntot","ndlt","overdose"); colnames(dat) = 1:ndose
  return(dat)  
}
#simulate.boin.trial(get.dlt.probs(6)$dlt.probs,boin.int)
