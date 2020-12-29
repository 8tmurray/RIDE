### Functions pertaining to the Retainment Interval Dose Escalation (RIDE) Design
### Written by Thomas A Murray on 08/16/2019 during manuscript revision

### Function to get next dose of loRIDE design
# ntot - number of patients treated at the current dose
# ndlt - number of DLTs at the current dose
# target - target DLT rate
# or - odds ratio of DLT for adjacent doses (default value is or = orp = orm = 2)
# orp - odds ratio of DLT with adjacent higher dose relative to current dose (supplied value overrides or)
# orm - odds ratio of DLT with current dose relative to adjacent lower dose (supplied value overrides or)
# nu - degrees of freedom for t prior (default is 7)
# pess - prior effective sample size (default is 2)
# sig - scale for t prior, overrides pess (defualt is NULL which leads to the value that equates to default pess)
### library(cubature) hcubature(,vectorize=TRUE) and library(distrEx) distrExIntegrate() are other options
### Secret is to calculate log-likelihood and then exponentiate it, otherwise you get numerical intregration errors
### I checked the values against a Gibbs sampler and they match!
### I also considered L2 loss which is slightly more aggressive in terms of escalation/de-escalation
get.next.dose <- function(ntot,ndlt,target,or=2,orp=NULL,orm=NULL,nu=7,pess=2,sig=NULL){
  
  #Default OR+ and OR-
  if(is.null(orp)) orp=or
  if(is.null(orm)) orm=or
  
  #Calculate scale (sig) of t-distribution based on pess when sig is not specified
  if(is.null(sig)) sig = sqrt((nu+1)/nu/target/(1-target)/pess)

  # Calculate posterior loss associated with each decision
  cst = integrate(function(x) exp(lchoose(ntot,ndlt)+ndlt*x-ntot*log(1+exp(x)))*dt((x-log(target/(1-target)))/sig,nu)/sig,-Inf,Inf)$value
  ple = integrate(function(x) abs(x+log(orp)-log(target/(1-target)))*exp(lchoose(ntot,ndlt)+ndlt*x-ntot*log(1+exp(x)))*dt((x-log(target/(1-target)))/sig,nu)/sig,-Inf,Inf)$value/cst
  plr = integrate(function(x) abs(x-log(target/(1-target)))*exp(lchoose(ntot,ndlt)+ndlt*x-ntot*log(1+exp(x)))*dt((x-log(target/(1-target)))/sig,nu)/sig,-Inf,Inf)$value/cst
  pld = integrate(function(x) abs(x-log(orm)-log(target/(1-target)))*exp(lchoose(ntot,ndlt)+ndlt*x-ntot*log(1+exp(x)))*dt((x-log(target/(1-target)))/sig,nu)/sig,-Inf,Inf)$value/cst
  
  # Determine which decision minimizes posterior loss
  decision = c("Escalate","Retain","Deescalate")[which.min(c(ple,plr,pld))]
  
  return(decision)  
}
#Check out decisions under various prior configurations, retainment interval halfwidths and outcomes at current dose
#ntot=3; target = 0.3; or = 2; sapply(0:ntot,function(x) get.next.dose(ntot,x,target=target,or=or,pess=2))


### Function to get the dose exploration/boundary guide for loRIDE
# ns - sample sizes at which you want the boundary
# target - target DLT rate
# or - odds ratio of DLT for adjacent doses (default value is or = orp = orm = 2)
# orp - odds ratio of DLT with adjacent higher dose relative to current dose (supplied value overrides or)
# orm - odds ratio of DLT with current dose relative to adjacent lower dose (supplied value overrides or)
# nu - degrees of freedom for t prior (default is 7)
# pess - prior effective sample size (default is 2)
# sig - scale for t prior, overrides pess (defualt is NULL which leads to the value that equates to default pess)
#Returns a|b, where y <= a -> Escalate and y >= b -> De-escalate, so a < y < b -> Retain 
get.boundary = function(ns,target,or=2,orp=NULL,orm=NULL,pess=2,sig=NULL){
  bnd = rbind(sapply(ns,function(n){
    foo = sapply(0:n,function(y) get.next.dose(ntot = n, ndlt = y, target = target, or = or, orp = orp, orm = orm, pess = pess, sig = sig))
    a = sum(foo=="Escalate")
    b = n-sum(foo=="Deescalate")
    if(a==0 & b==n) return("-|-")
    if(a==0 & b<n) return(paste("-|",b+1,sep=""))
    if(a>0 & b==n) return(paste(a-1,"|-",sep=""))
    if(a>0 & b<n) return(paste(a-1,"|",b+1,sep=""))
  }))
  return(bnd)
}




### Function to simulate an entire loRIDE trial
# dlt.probs - True DLT probabilities
# start.dose - Starting dose level
# cohort.size - cohort size
# ncohorts - maximum number of cohorts
# target - target DLT rate
# or - odds ratio of DLT for adjacent doses (default value is or = orp = orm = 2)
# orp - odds ratio of DLT with adjacent higher dose relative to current dose (supplied value overrides or)
# orm - odds ratio of DLT with current dose relative to adjacent lower dose (supplied value overrides or)
# nu - degrees of freedom for t prior (default is 7)
# pess - prior effective sample size (default is 2)
# sig - scale for t prior, overrides pess (defualt is NULL which leads to the value that equates to default pess)
# threshold - posterior probability threshold to declare over-dose (default value is 0.95)
simulate.trial = function(dlt.probs,target,start.dose=1,cohort.size=3,ncohorts=12,min.esc.size=1,or=2,orp=NULL,orm=NULL,nu=7,pess=2,sig=NULL,threshold=0.975){

  ndose = length(dlt.probs); n = y = rep(0,ndose); 
  
  curr.dose = start.dose; overdose = rep(FALSE,ndose); # underdose = rep(FALSE,ndose)
  while(sum(n) < ncohorts*cohort.size){
    #Generate new observations at the current dose
    ntot = n[curr.dose] = n[curr.dose] + cohort.size
    ndlt = y[curr.dose] = y[curr.dose] + rbinom(1,cohort.size,dlt.probs[curr.dose])
    
    #Assess safety criterion then decision rule
    if(pbeta(target,ndlt+1,ntot-ndlt+1,lower.tail=FALSE) > threshold & n[curr.dose] >= 3) overdose[curr.dose:ndose] = TRUE
    #if(pbeta(target,ndlt+2*target,ntot-ndlt+2*(1-target),lower.tail=FALSE) > threshold) overdose[curr.dose:ndose] = TRUE
    if(overdose[1]) break #Stop trial when all doses are overdoses
    
    decision = get.next.dose(ntot,ndlt,target=target,or=or,orp=orp,orm=orm,nu=nu,pess=pess,sig=sig)
    if((decision == "Deescalate" | overdose[curr.dose]) & curr.dose > 1) curr.dose = curr.dose-1
    if(decision == "Escalate" & !overdose[min(curr.dose+1,ndose)] & curr.dose < ndose & n[curr.dose] >= min.esc.size) curr.dose = curr.dose+1
  }
  dat = rbind(n,y,overdose); rownames(dat) = c("ntot","ndlt","overdose"); colnames(dat) = 1:ndose
  return(dat)  
}
#simulate.trial(dlt.probs=c(0.05,0.10,0.25,0.45,0.60,0.85),target=0.2)
#system.time(sapply(1:1000,function(iter) simulate.trial(dlt.probs=c(0.05,0.10,0.25,0.45,0.60,0.85),target=0.2)))

# Function to run PAVA algorithm which is used in the get.best.dose() function defined next
pava <- function(x, wt = rep(1, length(x))){
  n <- length(x)
  if (n <= 1) return(x)
  if (any(is.na(x)) || any(is.na(wt))) stop("Missing values in 'x' or 'wt' not allowed")
  lvlsets <- (1:n)
  repeat {
    viol <- (as.vector(diff(x)) < 0)
    if (!(any(viol))) break
    i <- min((1:(n - 1))[viol])
    lvl1 <- lvlsets[i]
    lvl2 <- lvlsets[i + 1]
    ilvl <- (lvlsets == lvl1 | lvlsets == lvl2)
    x[ilvl] <- sum(x[ilvl] * wt[ilvl])/sum(wt[ilvl])
    lvlsets[ilvl] <- lvl1
  }
  x
}


### Function to select best dose amongst safe and tested doses based on final data
#dat is a 3xndose matrix, row 1 is ntot, row 2 is ndlt, row 3 indicates overdoses
#target - target DLT probability
get.best.dose = function(dat,target){
  #Grab data
  ntot = dat[1,]; ndlt = dat[2,]; 
  
  #Run PAVA to estimate DLT probabilities
  pi.hat <- pava((ndlt[ntot>0]+0.005)/(ntot[ntot>0]+0.01), (ntot[ntot>0]+0.01)**2*(ntot[ntot>0]+1.01)/((ndlt[ntot>0]+0.05)*(ntot[ntot>0]-ndlt[ntot>0]+0.05)))
  if(sum(ntot>0)>1) pi.hat[2:sum(ntot>0)] <- pi.hat[2:sum(ntot>0)] + (2:sum(ntot>0))*1E-10 # break ties so that largest dose tied below MTD is selected, and lowest dose tied above the MTD is selected
  
  #Determine Best Dose
  overdose = dat[3,]; #underdose = dat[4,]
  if(overdose[1]==1) best.dose = 0 #All overdoses
  if(overdose[1]==0) best.dose = which.min(abs(pi.hat[overdose==0]-target)) #Closest to the target among safe doses
  
  return(as.numeric(best.dose))
}
#get.best.dose(dat=simulate.trial(dlt.probs=c(0.05,0.10,0.25,0.45,0.60,0.85),target=0.25),target=0.25)