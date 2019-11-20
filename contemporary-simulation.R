### Simulation Study based on New Scenario Generation Algorithm
### Written by Thomas A Murray on 8/16/2019 
setwd("C:/Users/Thomas/Dropbox/Research/Active/RIDE/Software") #Set working directory
source('ride-functions.R') #Load RIDE functions
source('comparator-functions.R') #Load comparator functions

### Function to generate true toxicity probabilities
# ndose - number of investigational dose levels
# scen - type of scenario 
#     NULL = random draw (default)
#     0 = all overdoses (delta above target)
#     1 to ndose = corresponding dose level is closest to target (and within delta)
#     ndose+1 = all underdoses (delta below target)
# target - target DLT rate
# delta - maximum deviation from target in cases with an acceptable dose and minimum deviation from target in cases without an acceptable dose
get.dlt.probs = function(ndose, target, delta, scen=NULL){
  if(is.null(scen)) scen = sample(0:(ndose+1),1)
  if(scen==0) dlt.probs = sort(runif(ndose,runif(1,target+delta,1),1))
  if(scen %in% 1:ndose){
    dlt.probs = rep(NA,ndose)
    dlt.probs[scen] = runif(1,target-delta,target+delta)
    if(scen > 1) dlt.probs[1:(scen-1)] = sort(runif(scen-1,0,target-abs(dlt.probs[scen]-target)))
    if(scen < ndose) dlt.probs[(scen+1):ndose] = sort(runif(ndose-scen,target+abs(dlt.probs[scen]-target),runif(1,target+abs(dlt.probs[scen]-target),1)))
  }
  if(scen==ndose+1) dlt.probs = sort(runif(ndose,0,runif(1,0,target-delta)))
  
  return(list(scen=scen,dlt.probs=dlt.probs))
}

#Function to calculate summary statistics based on data and dose selection
# ldat - list of datasets
# scen - type of scenario
# dlt.probs - DLT probabilities
# target - target DLT Probability
# delta - acceptable deviation from target for adequate dosing
get.stats <- function(ldat,scen,dlt.probs,target,delta,nmax,crm=FALSE,skeleton=NULL){

  #Scenario information
  ndose = length(dlt.probs)
  best.dose = ifelse(scen == ndose+1,ndose,scen)
  dlt.probs = c(0,dlt.probs)
  
  #Dose Selection Results
  if(!crm) psel = prop.table(table(factor(sapply(ldat,function(dat) get.best.dose(dat,target=target)),levels=0:ndose)))
  if(crm) psel = prop.table(table(factor(sapply(ldat,function(dat) get.crm.best.dose(dat,skeleton=skeleton,target=target)),levels=0:ndose)))
  sbd = psel[best.dose+1] #% Selected Best Dose
  sad = sum(psel[abs(dlt.probs-target)<=delta]) #% selected adequate dose (within \delta of target)
  sod = sum(psel[dlt.probs>target+delta]) #% selected an overdose (more than \delta above target)
  sud = sum(psel[dlt.probs<target-delta]) #% selected an underdose (more than \delta above target)
  sai = 1-length(dlt.probs)*(sum(abs(dlt.probs-target)*psel)/sum(abs(dlt.probs-target))) #Wages Accuracy Index

  #Dose Allocation Results
  mdlt = mean(sapply(ldat,function(dat) sum(dat[2,]))) #Average Number of DLTs
  ptrts = sapply(ldat,function(x) c(nmax-sum(x[1,]),x[1,]))/nmax #% given each dose (including dose 0 = stopped for safety)
  ptrt = rowMeans(ptrts) #Average % given each dose
  tbd = ptrt[best.dose+1] #Average % given best dose
  tad = sum(ptrt[abs(dlt.probs-target)<=delta]) #Average % given adequate dose (within \delta of target)
  tod = sum(ptrt[dlt.probs>target+delta]) #Average % given an overdose (more than \delta above target)
  tud = sum(ptrt[dlt.probs<target-delta]) #Average % given an underdose (more than \delta above target)
  tai = 1-length(dlt.probs)*(sum(abs(dlt.probs-target)*ptrt)/sum(abs(dlt.probs-target))) #Wages Accuracy Index
  pod40 = pod70 = pud40 = pud70 = NA
  if(sum(dlt.probs[-1]>target+delta)>0){ pod = colSums(rbind(ptrts[dlt.probs>target+delta,])); pod40 = mean(pod >= 0.4); pod70 = mean(pod >= 0.7) } #% of trials where >=X% given overdose
  if(sum(dlt.probs[-1]<target-delta)>0){ pud = colSums(rbind(ptrts[dlt.probs<target-delta,])); pud40 = mean(pud >= 0.4); pud70 = mean(pud >= 0.7) } #% of trials where >=X% given underdose

  stats = rbind(c(sbd,sad,sod,sud,sai,mdlt,tbd,tad,tod,tud,tai,pod40,pod70,pud40,pud70))
  colnames(stats) = c("sBD","sAD","sOD","sUD","sAI","mDLT","tBD","tAD","tOD","tUD","tAI","pOD40","pOD70","pUD40","pUD70")
  return(stats)
}


##################
# Run Simulation #
##################

#Simulation Parameters
targets = c(0.1,0.2,0.3) #Target DLT Probability
ndose = 6 #Number of investigational dose levels
cohort.size = 3 #Cohort Size
ncohorts.vec = c(8,12) #Maximum number of cohorts
delta = 0.07 #Deviation from target that defines overdose
nrep = 100 #Number of times each type of scenario is repeated
ntrials = 2000 #Number of simulated trials for each true scenario

#Parallely loop over every true scenario
library(foreach); library(doParallel); library(doRNG)
RNGkind("L'Ecuyer"); registerDoParallel(cores=8)
#Loop over targets
for(target in targets){

  #Comparator Design parameters
  crm.skeleton = get.skeleton(ndose,target,delta) #CRM skeleton
  boin.int = as.numeric(BOIN::get.boundary(target,1,3)[1:2]) #BOIN dose escalation interval
  
  #Generate <nscens> sets of true DLT probabilties via pseudo-uniform algorithm
  set.seed(1985) #This allows us to recover true DLT probabilities in each scenario
  scen.type = rep(0:(ndose+1),each=nrep) #Scenario Type
  scenarios = lapply(scen.type,function(scen) get.dlt.probs(ndose=ndose,target=target,delta=delta,scen=scen)) #True Randomly Generated Scenarios

  #Loop over ncohorts.vec  
  for(ncohorts in ncohorts.vec){
    
    #Loop over scenarios
    for(iter in 1:length(scen.type)){

      ###True DLT Probabilities
      scen = scenarios[[iter]]$scen
      dlt.probs = scenarios[[iter]]$dlt.probs
      
      ###Simulate <ntrials> trials for each design and get summary statistics
      #RIDE (with default working value OR=1.6)
      #ride1.dat = lapply(1:100,function(m) simulate.trial(dlt.probs,cohort.size=cohort.size,ncohorts=ncohorts,target=target,or=1.60,min.esc.size=3))
      ride1.dat = foreach(m=1:ntrials) %dorng% { simulate.trial(dlt.probs,cohort.size=cohort.size,ncohorts=ncohorts,target=target,or=1.60,min.esc.size=3) }
      ride1.stats = get.stats(ride1.dat,scen,dlt.probs,nmax=ncohorts*cohort.size,target=target,delta=delta)
      write.table(cbind(rbind(iter),round(ride1.stats,4)),sep=",",row.names=FALSE,
                  file=paste("Results/Contemporary/RIDE1-Stats-N",ncohorts*cohort.size,"-P",100*target,".txt",sep=""),
                  col.names=!file.exists(paste("Results/Contemporary/RIDE1-Stats-N",ncohorts*cohort.size,"-P",100*target,".txt",sep="")),
                  append=file.exists(paste("Results/Contemporary/RIDE1-Stats-N",ncohorts*cohort.size,"-P",100*target,".txt",sep="")))
      
      #RIDE (with smaller than default working value OR=1.05)
      ride2.dat = foreach(m=1:ntrials) %dorng% { simulate.trial(dlt.probs,cohort.size=cohort.size,ncohorts=ncohorts,target=target,or=1.05,min.esc.size=3) }
      ride2.stats = get.stats(ride2.dat,scen,dlt.probs,nmax=ncohorts*cohort.size,target=target,delta=delta)
      write.table(cbind(rbind(iter),round(ride2.stats,4)),sep=",",row.names=FALSE,
                  file=paste("Results/Contemporary/RIDE2-Stats-N",ncohorts*cohort.size,"-P",100*target,".txt",sep=""),
                  col.names=!file.exists(paste("Results/Contemporary/RIDE2-Stats-N",ncohorts*cohort.size,"-P",100*target,".txt",sep="")),
                  append=file.exists(paste("Results/Contemporary/RIDE2-Stats-N",ncohorts*cohort.size,"-P",100*target,".txt",sep="")))
      
      #RIDE (with larger than default working value OR=4)
      ride3.dat = foreach(m=1:ntrials) %dorng% { simulate.trial(dlt.probs,cohort.size=cohort.size,ncohorts=ncohorts,target=target,or=4.00,min.esc.size=3) }
      ride3.stats = get.stats(ride3.dat,scen,dlt.probs,nmax=ncohorts*cohort.size,target=target,delta=delta)
      write.table(cbind(rbind(iter),round(ride3.stats,4)),sep=",",row.names=FALSE,
                  file=paste("Results/Contemporary/RIDE3-Stats-N",ncohorts*cohort.size,"-P",100*target,".txt",sep=""),
                  col.names=!file.exists(paste("Results/Contemporary/RIDE3-Stats-N",ncohorts*cohort.size,"-P",100*target,".txt",sep="")),
                  append=file.exists(paste("Results/Contemporary/RIDE3-Stats-N",ncohorts*cohort.size,"-P",100*target,".txt",sep="")))
      
      #mTPI
      mtpi.dat = foreach(m=1:ntrials) %dorng% { simulate.mtpi.trial(dlt.probs,cohort.size=cohort.size,ncohorts=ncohorts,target=target,delta=0.07,min.esc.size=3) }
      mtpi.stats = get.stats(mtpi.dat,scen,dlt.probs,nmax=ncohorts*cohort.size,target=target,delta=delta)
      write.table(cbind(rbind(iter),round(mtpi.stats,4)),sep=",",row.names=FALSE,
                  file=paste("Results/Contemporary/mTPI-Stats-N",ncohorts*cohort.size,"-P",100*target,".txt",sep=""),
                  col.names=!file.exists(paste("Results/Contemporary/mTPI-Stats-N",ncohorts*cohort.size,"-P",100*target,".txt",sep="")),
                  append=file.exists(paste("Results/Contemporary/mTPI-Stats-N",ncohorts*cohort.size,"-P",100*target,".txt",sep="")))
      
      #CRM (with default skeleton)
      crm.dat = foreach(m=1:ntrials) %dorng% { simulate.crm.trial(dlt.probs,skeleton=crm.skeleton,cohort.size=cohort.size,ncohorts=ncohorts,target=target,min.esc.size=3) }
      crm.stats = get.stats(crm.dat,scen,dlt.probs,nmax=ncohorts*cohort.size,target=target,delta=delta,crm=TRUE,skeleton=crm.skeleton)
      write.table(cbind(rbind(iter),round(crm.stats,4)),sep=",",row.names=FALSE,
                  file=paste("Results/Contemporary/CRM-Stats-N",ncohorts*cohort.size,"-P",100*target,".txt",sep=""),
                  col.names=!file.exists(paste("Results/Contemporary/CRM-Stats-N",ncohorts*cohort.size,"-P",100*target,".txt",sep="")),
                  append=file.exists(paste("Results/Contemporary/CRM-Stats-N",ncohorts*cohort.size,"-P",100*target,".txt",sep="")))
      
      #BOIN
      boin.dat = foreach(m=1:ntrials) %dorng% { simulate.boin.trial(dlt.probs,boin.int=boin.int,cohort.size=cohort.size,ncohorts=ncohorts,target=target,min.esc.size=3) }
      boin.stats = get.stats(boin.dat,scen,dlt.probs,nmax=ncohorts*cohort.size,target=target,delta=delta)
      write.table(cbind(rbind(iter),round(boin.stats,4)),sep=",",row.names=FALSE,
                  file=paste("Results/Contemporary/BOIN-Stats-N",ncohorts*cohort.size,"-P",100*target,".txt",sep=""),
                  col.names=!file.exists(paste("Results/Contemporary/BOIN-Stats-N",ncohorts*cohort.size,"-P",100*target,".txt",sep="")),
                  append=file.exists(paste("Results/Contemporary/BOIN-Stats-N",ncohorts*cohort.size,"-P",100*target,".txt",sep="")))
      
      #Print Scen Number
      print(paste("     Finished Scenario ",iter," of ",length(scen.type),sep=""))
    }
    #Print Scen Number
    print(paste("Finished Scenarios: N=",ncohorts*cohort.size,", P=",100*target,sep=""))
  }
}
stopImplicitCluster()



###################
# Analyze Results #
###################
#Plot summarizing scenarios
foo = t(sapply(1:length(scenarios),function(i) scenarios[[i]]$dlt.probs)) #Matrix of True DLT Probabilities
rs = unlist(lapply(0:7*100,function(x) sample(x+1:100,5))) #Random subset of the scenarios (10 of each type of scenario)

pdf(file="Scenario-Boxplot.pdf")
par(lab=c(6,10,7),family="serif")
boxplot(foo,las=1,ylab=expression(paste("DLT Probability ( ",pi[j]," )"),sep=""),xlab="Dose Level (j)",main="Scenario Coverage")
abline(h=target,col="grey",lwd=2)
dev.off()


#Analyze results
library(dplyr); library(xtable)
target = 0.2; ncohorts = 12; cohort.size = 3
ride1.stats = read.table(paste("Results/Contemporary/RIDE1-Stats-N",ncohorts*cohort.size,"-P",100*target,".txt",sep=""),sep=",",header=TRUE)
ride2.stats = read.table(paste("Results/Contemporary/RIDE2-Stats-N",ncohorts*cohort.size,"-P",100*target,".txt",sep=""),sep=",",header=TRUE)
ride3.stats = read.table(paste("Results/Contemporary/RIDE3-Stats-N",ncohorts*cohort.size,"-P",100*target,".txt",sep=""),sep=",",header=TRUE)
mtpi.stats = read.table(paste("Results/Contemporary/mTPI-Stats-N",ncohorts*cohort.size,"-P",100*target,".txt",sep=""),sep=",",header=TRUE)
crm.stats = read.table(paste("Results/Contemporary/CRM-Stats-N",ncohorts*cohort.size,"-P",100*target,".txt",sep=""),sep=",",header=TRUE)
boin.stats = read.table(paste("Results/Contemporary/BOIN-Stats-N",ncohorts*cohort.size,"-P",100*target,".txt",sep=""),sep=",",header=TRUE)

#Averaging across all scenarios
tab <- rbind(colMeans(ride1.stats[,-1],na.rm = TRUE),
      colMeans(ride2.stats[,-1],na.rm = TRUE),
      colMeans(ride3.stats[,-1],na.rm = TRUE),
      colMeans(mtpi.stats[,-1],na.rm = TRUE),
      colMeans(crm.stats[,-1],na.rm = TRUE),
      colMeans(boin.stats[,-1],na.rm = TRUE))
rownames(tab) = c("RIDE","RI(-)","RI(+)","mTPI","CRM","BOIN")
print.xtable(xtable(t(tab),digits=3),include.rownames=TRUE)

#Averaging across scenarios of specified type
type = 2; scens = 1:100 + 100*type
tab <- rbind(colMeans(ride1.stats[scens,-1],na.rm = TRUE),
             colMeans(ride2.stats[scens,-1],na.rm = TRUE),
             colMeans(ride3.stats[scens,-1],na.rm = TRUE),
             colMeans(mtpi.stats[scens,-1],na.rm = TRUE),
             colMeans(crm.stats[scens,-1],na.rm = TRUE),
             colMeans(boin.stats[scens,-1],na.rm = TRUE))
rownames(tab) = c("RIDE","RI(-)","RI(+)","mTPI","CRM","BOIN")
print.xtable(xtable(t(tab),digits=3),include.rownames=TRUE)

