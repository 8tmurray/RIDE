### This program generates the dose exploration guides for Table 1 and safety guides for Table 2
### Written by Thomas A Murray on 10/03/2019
setwd("C:/Users/Thomas/Dropbox/Research/Active/RIDE/Software") #Set working directory
source("comparator-functions.R")
source("ride-functions.R")



### function to get decision boundaries at ns
#Returns a|b, where y <= a -> Escalate and y >= b -> De-escalate, so a < y < b -> Retain 
get.ride.boundary = function(ns,target,or=2,orp=NULL,orm=NULL,pess=2,sig=NULL){
  bnd = rbind(sapply(ns,function(n){
    foo = sapply(0:n,function(y) get.next.dose(ntot = n, ndlt = y, target = target, or = or, orp = orp, orm = orm, pess = pess, sig = sig))
    a = sum(foo=="Escalate")
    b = n-sum(foo=="Deescalate")
    if(a==0 & b==n) return("-$\\vert$-")
    if(a==0 & b<n) return(paste("-$\\vert$",b+1,sep=""))
    if(a>0 & b==n) return(paste(a-1,"$\\vert$-",sep=""))
    if(a>0 & b<n) return(paste(a-1,"$\\vert$",b+1,sep=""))
  }))
  return(bnd)
}
#get.ride.boundary(ns,target)

#Dose Guide for default BOIN
get.boin.boundary = function(ns,target){
  boin.int = as.numeric(BOIN::get.boundary(target=target,ncohort=1,cohortsize=2)[1:2])
  bnd = sapply(ns,function(n){
    foo = sapply(0:n,function(y) get.boin.next.dose(ntot = n, ndlt = y, boin.int = boin.int))
    a = sum(foo=="Escalate")
    b = n-sum(foo=="Deescalate")
    if(a==0 & b==n) return("-$\\vert$-")
    if(a==0 & b<n) return(paste("-$\\vert$",b+1,sep=""))
    if(a>0 & b==n) return(paste(a-1,"$\\vert$-",sep=""))
    if(a>0 & b<n) return(paste(a-1,"$\\vert$",b+1,sep=""))
  })
  return(bnd)
}
#get.boin.boundary(ns,target)


#Dose Guide for mTPI
get.mtpi.boundary = function(ns,target,delta=0.25*target){
  bnd = sapply(ns,function(n){
    foo = sapply(0:n,function(y) get.mtpi.next.dose(ntot = n, ndlt = y, target = target, delta = delta))
    a = sum(foo=="Escalate")
    b = n-sum(foo=="Deescalate")
    if(a==0 & b==n) return("-$\\vert$-")
    if(a==0 & b<n) return(paste("-$\\vert$",b+1,sep=""))
    if(a>0 & b==n) return(paste(a-1,"$\\vert$-",sep=""))
    if(a>0 & b<n) return(paste(a-1,"$\\vert$",b+1,sep=""))
  })
  return(bnd)
}
#get.mtpi.boundary(ns,target)


#Dose Guide for default Keyboard
get.keyboard.boundary = function(ns,target,delta=0.25*target){
  bnd = sapply(ns,function(n){
    foo = sapply(0:n,function(y) get.keyboard.next.dose(ntot = n, ndlt = y, target = target, delta = delta))
    a = sum(foo=="Escalate")
    b = n-sum(foo=="Deescalate")
    if(a==0 & b==n) return("-$\\vert$-")
    if(a==0 & b<n) return(paste("-$\\vert$",b+1,sep=""))
    if(a>0 & b==n) return(paste(a-1,"$\\vert$-",sep=""))
    if(a>0 & b<n) return(paste(a-1,"$\\vert$",b+1,sep=""))
  })
  return(bnd)
}
#get.keyboard.boundary(ns,target)

1/(1+exp(-(log(0.1/0.9)+c(-1,1)*log(3.4)/2)))

### Table 1, Dose Exploration Guides
ns = c(1:9,12,18,24)
tab3 = rbind(get.boin.boundary(ns=ns,target=0.3),
      get.mtpi.boundary(ns=ns,target=0.3),
      get.keyboard.boundary(ns=ns,target=0.3),
      get.ride.boundary(ns=ns,target=0.3),
      get.ride.boundary(ns=ns,target=0.3,or=3.4),
      get.ride.boundary(ns=ns,target=0.3,or=1.2),
      get.ride.boundary(ns=ns,target=0.3,orp=2,orm=3.4),
      get.ride.boundary(ns=ns,target=0.3,orp=3.4,orm=2),
      get.ride.boundary(ns=ns,target=0.3,orp=2,orm=1.2),
      get.ride.boundary(ns=ns,target=0.3,orp=1.2,orm=2))
rownames(tab3) = c("BOIN","mTPI","Keyboard","loRIDE","loRIDE(+,+)","loRIDE(-,-)","loRIDE(=,+)","loRIDE(+,=)","loRIDE(=,-)","loRIDE(-,=)")
tab2 = rbind(get.boin.boundary(ns=ns,target=0.2),
             get.mtpi.boundary(ns=ns,target=0.2),
             get.keyboard.boundary(ns=ns,target=0.2),
             get.ride.boundary(ns=ns,target=0.2))
rownames(tab2) = c("BOIN","mTPI","Keyboard","loRIDE")
tab1 = rbind(get.boin.boundary(ns=ns,target=0.1),
             get.mtpi.boundary(ns=ns,target=0.1),
             get.keyboard.boundary(ns=ns,target=0.1),
             get.ride.boundary(ns=ns,target=0.1))
rownames(tab1) = c("BOIN","mTPI","Keyboard","loRIDE")
colnames(tab3) = colnames(tab2) = colnames(tab1) = ns

write.table(rbind(tab3,tab2,tab1),sep=" & ",quote=F)

tab3 = rbind(get.boin.boundary(ns=ns,target=0.2),
             get.mtpi.boundary(ns=ns,target=0.2),
             get.keyboard.boundary(ns=ns,target=0.2),
             get.ride.boundary(ns=ns,target=0.2),
             get.ride.boundary(ns=ns,target=0.2,or=3.4),
             get.ride.boundary(ns=ns,target=0.2,or=1.2),
             get.ride.boundary(ns=ns,target=0.2,orp=2,orm=3.4),
             get.ride.boundary(ns=ns,target=0.2,orp=3.4,orm=2),
             get.ride.boundary(ns=ns,target=0.2,orp=2,orm=1.2),
             get.ride.boundary(ns=ns,target=0.2,orp=1.2,orm=2))
rownames(tab3) = c("BOIN","mTPI","Keyboard","loRIDE","loRIDE(+,+)","loRIDE(-,-)","loRIDE(=,+)","loRIDE(+,=)","loRIDE(=,-)","loRIDE(-,=)")
write.table(tab3,sep=" & ",quote=F)


#rbind(get.boundary(ns=ns,target=0.1,or=2),get.boundary(ns=ns,target=0.1,or=1.2),get.boundary(ns=ns,target=0.1,or=3.4))







### Safety Rule Table (Table 2 in manuscript)
unsafe.prob = function(ntot,ndlt,target,pess=2,sig=NULL,nu=7){
  #Calculate scale (sig) of t-distribution based on pess when sig is not specified
  if(is.null(sig)) sig = sqrt((nu+1)/nu/target/(1-target)/pess)
  
  cst = integrate(function(x) exp(lchoose(ntot,ndlt)+ndlt*x-ntot*log(1+exp(x)))*dt((x-log(target/(1-target)))/sig,nu)/sig,-Inf,Inf)$value
  integrate(function(x) exp(lchoose(ntot,ndlt)+ndlt*x-ntot*log(1+exp(x)))*dt((x-log(target/(1-target)))/sig,nu)/sig/cst,log(target/(1-target)),Inf)$value  
}

ns = c(1:9,12,18,24)
tab = rbind(paste(ns,collapse=" & "),  
            paste(sapply(ns,function(n) which(sapply(0:n, function(y) unsafe.prob(n,y,target=0.3))>0.975)[1]-1),collapse=" & "),
            paste(sapply(ns,function(n) which(sapply(0:n, function(y) unsafe.prob(n,y,target=0.2))>0.975)[1]-1),collapse=" & "),
            paste(sapply(ns,function(n) which(sapply(0:n, function(y) unsafe.prob(n,y,target=0.1))>0.975)[1]-1),collapse=" & "),
            paste(sapply(ns,function(n) which(sapply(0:n, function(y) unsafe.prob(n,y,target=0.3))>0.95)[1]-1),collapse=" & "),
            paste(sapply(ns,function(n) which(sapply(0:n, function(y) unsafe.prob(n,y,target=0.2))>0.95)[1]-1),collapse=" & "),
            paste(sapply(ns,function(n) which(sapply(0:n, function(y) unsafe.prob(n,y,target=0.1))>0.95)[1]-1),collapse=" & "))
tab = rbind(sapply(ns,function(n) which(sapply(0:n, function(y) unsafe.prob(n,y,target=0.3))>0.975)[1]-1),
            sapply(ns,function(n) which(sapply(0:n, function(y) unsafe.prob(n,y,target=0.2))>0.975)[1]-1),
            sapply(ns,function(n) which(sapply(0:n, function(y) unsafe.prob(n,y,target=0.1))>0.975)[1]-1),
            sapply(ns,function(n) which(sapply(0:n, function(y) unsafe.prob(n,y,target=0.3))>0.95)[1]-1),
            sapply(ns,function(n) which(sapply(0:n, function(y) unsafe.prob(n,y,target=0.2))>0.95)[1]-1),
            sapply(ns,function(n) which(sapply(0:n, function(y) unsafe.prob(n,y,target=0.1))>0.95)[1]-1))
colnames(tab) = ns; rownames(tab) = rep(c(0.30,0.20,0.10),2)
write.table(tab,sep=" & ",quote=F)
