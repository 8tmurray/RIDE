### This program generates the dose exploration guides for Table 1 and safety guides for Table 2
### Written by Thomas A Murray on 10/03/2019
setwd("C:/Users/Thomas/Dropbox/Research/Active/RIDE/Software") #Set working directory
source("comparator-functions.R")
source("ride-functions.R")



### function to get decision boundaries at ns
#Returns a:b, where y <= a -> Escalate and y >= b -> De-escalate, so a < y < b -> Retain 
get.boundary = function(ns,target,or=1.6,orp=NULL,orm=NULL){

  #Default OR+ and OR-
  if(is.null(orp)) orp=or
  if(is.null(orm)) orm=or
  
  #Get Boundary
  bnd = rbind(sapply(ns,function(n){
    foo = sapply(0:n,function(y) get.next.dose(ntot = n, ndlt = y, target = target, orp = orp, orm = orm))
    a = sum(foo=="Escalate")
    b = n-sum(foo=="Deescalate")
    if(a==0 & b==n) return("-/-")
    if(a==0 & b<n) return(paste("-/",b+1,sep=""))
    if(a>0 & b==n) return(paste(a-1,"/-",sep=""))
    if(a>0 & b<n) return(paste(a-1,"/",b+1,sep=""))
  }))
  return(bnd)
}


#Dose Guide for default BOIN and mTPI Designs
get.boin.boundary = function(ns,target){
  boin.int = as.numeric(BOIN::get.boundary(target=target,ncohort=1,cohortsize=2)[1:2])
  bnd = sapply(ns,function(n){
    foo = sapply(0:n,function(y) get.boin.next.dose(ntot = n, ndlt = y, boin.int = boin.int))
    a = sum(foo=="Escalate")
    b = n-sum(foo=="Deescalate")
    if(a==0 & b==n) return("-/-")
    if(a==0 & b<n) return(paste("-/",b+1,sep=""))
    if(a>0 & b==n) return(paste(a-1,"/-",sep=""))
    if(a>0 & b<n) return(paste(a-1,"/",b+1,sep=""))
  })
  return(bnd)
}
boin.guide = get.boin.boundary(ns,target=0.3)

get.mtpi.boundary = function(ns,target,delta){
  bnd = sapply(ns,function(n){
    foo = sapply(0:n,function(y) get.mtpi.next.dose(ntot = n, ndlt = y, target = target, delta = delta))
    a = sum(foo=="Escalate")
    b = n-sum(foo=="Deescalate")
    if(a==0 & b==n) return("-/-")
    if(a==0 & b<n) return(paste("-/",b+1,sep=""))
    if(a>0 & b==n) return(paste(a-1,"/-",sep=""))
    if(a>0 & b<n) return(paste(a-1,"/",b+1,sep=""))
  })
  return(bnd)
}
mtpi.guide = get.mtpi.boundary(ns,target,delta=0.07)


#RIDE Dose Guide
library(xtable)
ns = c(1:6,9,12,18,24); target = 0.1
foo = rbind(get.boin.boundary(ns=ns,target=target),
      get.mtpi.boundary(ns=ns,target=target,delta=0.07),     
      get.boundary(ns=ns,target=target,or=1.6),
      get.boundary(ns=ns,target=target,or=4),
      get.boundary(ns=ns,target=target,or=1.05),
      get.boundary(ns=ns,target=target,orm=1.6,orp=4),
      get.boundary(ns=ns,target=target,orm=4,orp=1.6)
)
colnames(foo) = ns
print(xtable(foo),include.rownames=FALSE)


### Table 2, safety rule
get.unsafe.boundary = function(ns,target,delta,prior){
  bnd = sapply(ns,function(n){
    foo = sapply(0:n,function(y) assess.adequacy(ntot = n, ndlt = y, target = target, delta = delta, prior = prior))
    cut = n-sum(foo=="Overdose")+1
    return(cut)
  })
  return(bnd)
}
table2 = t(sapply(c("mean","median","mode"),function(p) get.unsafe.boundary(ns,target,delta,prior=p)))
colnames(table2) = ns
table2

library(xtable)
xtable(table2,digits=0)





### Safety Rule Table (Table 2 in manuscript)
unsafe.prob = function(ntot,ndlt,target,sig=2.5,nu=7){
  cst = integrate(function(x) exp(lchoose(ntot,ndlt)+ndlt*x-ntot*log(1+exp(x)))*dt((x-log(target/(1-target)))/sig,nu)/sig,-Inf,Inf)$value
  integrate(function(x) exp(lchoose(ntot,ndlt)+ndlt*x-ntot*log(1+exp(x)))*dt((x-log(target/(1-target)))/sig,nu)/sig/cst,log(target/(1-target)),Inf)$value  
}
target=0.1; cut = 0.95; ns = c(1:6,9,12,18,24)
cat(sapply(ns,function(n) which(sapply(0:n, function(y) unsafe.prob(n,y,target))>cut)[1]-1),sep=" & ")

