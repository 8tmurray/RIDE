### Simulation based on "classic" five scenarios
### Written by Thomas A Murray on 09/12/2019 
setwd("C:/Users/Thomas/Dropbox/Research/Active/RIDE/Software") #Set working directory
source('ride-functions.R') #Load JADE functions
source('comparator-functions.R') #Load comparator functions

#Safety rule differences
#target = 0.2; a=2*target; b=2*(1-target)
#sapply(1:12,function(ntot) rbind(1-pbeta(target,0:ntot+2*target,ntot-0:ntot+2*(1-target)),1-pbeta(target,0:ntot+1,ntot-0:ntot+1)))

#Simulation Specifications (scenarios from Ivanova and Kim 2009)
target = 0.20    #Target DLT Probability
cohort.size = 1  #Cohort Size
ncohorts.vec = c(25,48)    #Maximum number of cohorts
min.esc.size = 3 #Minimum number of people treated at current dose before escalation is allowed
ntrials = 10000  #Number of simulated trials for each true scenario
#True DLT Probabilities in each scenario
dlt.probs = rbind(c(0.05,0.10,0.20,0.30,0.50,0.70), 
                  c(0.30,0.40,0.52,0.61,0.76,0.87),
                  c(0.05,0.06,0.08,0.11,0.19,0.34),
                  c(0.06,0.08,0.12,0.18,0.40,0.71),
                  c(0.00,0.00,0.03,0.05,0.11,0.22))
ndose = ncol(dlt.probs)

#Comparator Design parameters
crm.skeleton = get.skeleton(ndose,target,delta=0.07) #CRM skeleton
boin.int = as.numeric(BOIN::get.boundary(target=target,ncohort=3,cohortsize=1)[1:2]) #BOIN dose escalation interval


#Function to Calculate Summary Metrics
get.stats <- function(ldat,cohort.size,ncohorts,target,crm=FALSE,skeleton=NULL){
  
  # Number of doses
  ndose = ncol(ldat[[1]])
  
  # Average numbers assigned to each dose
  ptrt = rowMeans(sapply(ldat,function(x) c(cohort.size*ncohorts-sum(x[1,]),x[1,])))/(cohort.size*ncohorts)

  # Selection Percentages
  if(!crm) psel = prop.table(table(factor(sapply(ldat,function(dat) get.best.dose(dat,target=target)),levels=0:ndose)))
  if(crm)  psel = prop.table(table(factor(sapply(ldat,function(dat) get.crm.best.dose(dat,skeleton=skeleton,target=target)),levels=0:ndose)))
  
  # Collect statistics
  stats = rbind(ptrt=ptrt,psel=psel)
  colnames(stats) = 0:ndose
  
  # Return statistics
  return(stats)
}

#Function to calculate Wages Accuracy Index
get.ai = function(dlt.probs,target,ntrials,sel.pct){
  
  # Add Dose Level 0 == No Dose Given
  dlt.probs = c(0,dlt.probs)
  
  # Calculate AI
  ai = 1-length(dlt.probs)*(sum(abs(dlt.probs-target)*sel.pct)/sum(abs(dlt.probs-target)))
  
  # Return AI
  return(ai)
}

#Loop over every true scenario and simulate trials under each design parallely
library(foreach); library(doParallel); library(doRNG)
RNGkind("L'Ecuyer"); set.seed(1985)
registerDoParallel(cores=8)
for(ncohorts in ncohorts.vec){
  for(scen in 1:nrow(dlt.probs)){
    
    #RIDE (with default working value OR=2)
    ride1.dat = foreach(m=1:ntrials) %dorng% { simulate.trial(dlt.probs[scen,],cohort.size=cohort.size,ncohorts=ncohorts,target=target,or=1.60,min.esc.size=min.esc.size) }
    ride1.stats = get.stats(ncohorts,cohort.size,ldat=ride1.dat,target=target)
    ride1.ai = round(get.ai(dlt.probs[scen,],target=target,ntrials=ntrials,sel.pct=ride1.stats[2,]),3)
    write.table(rbind(pdlt=c(NA,dlt.probs[scen,]),ride1.stats,ai=c(ride1.ai,rep(NA,length(dlt.probs[scen,])))),sep=",",
                file=paste("Results/Classic/RIDE1-Stats-N",ncohorts,".txt",sep=""),
                col.names=!file.exists(paste("Results/Classic/RIDE1-Stats-N",ncohorts,".txt",sep="")),
                row.names=FALSE,append=file.exists(paste("Results/Classic/RIDE1-Stats-N",ncohorts,".txt",sep="")))
    
    #RIDE (with smaller than default working value OR=1.4)
    ride2.dat = foreach(m=1:ntrials) %dorng% { simulate.trial(dlt.probs[scen,],cohort.size=cohort.size,ncohorts=ncohorts,target=target,or=1.05,min.esc.size=min.esc.size) }
    ride2.stats = get.stats(ncohorts,cohort.size,ldat=ride2.dat,target=target)
    ride2.ai = round(get.ai(dlt.probs[scen,],target=target,ntrials=ntrials,sel.pct=ride2.stats[2,]),3)
    write.table(rbind(pdlt=c(NA,dlt.probs[scen,]),ride2.stats,ai=c(ride2.ai,rep(NA,length(dlt.probs[scen,])))),sep=",",
                file=paste("Results/Classic/RIDE2-Stats-N",ncohorts,".txt",sep=""),
                col.names=!file.exists(paste("Results/Classic/RIDE2-Stats-N",ncohorts,".txt",sep="")),
                row.names=FALSE,append=file.exists(paste("Results/Classic/RIDE2-Stats-N",ncohorts,".txt",sep="")))
    
    #RIDE (with larger than default working value OR=4)
    ride3.dat = foreach(m=1:ntrials) %dorng% { simulate.trial(dlt.probs[scen,],cohort.size=cohort.size,ncohorts=ncohorts,target=target,or=4.00,min.esc.size=min.esc.size) }
    ride3.stats = get.stats(ncohorts,cohort.size,ldat=ride3.dat,target=target)
    ride3.ai = round(get.ai(dlt.probs[scen,],target=target,ntrials=ntrials,sel.pct=ride3.stats[2,]),3)
    write.table(rbind(pdlt=c(NA,dlt.probs[scen,]),ride3.stats,ai=c(ride3.ai,rep(NA,length(dlt.probs[scen,])))),sep=",",
                file=paste("Results/Classic/RIDE3-Stats-N",ncohorts,".txt",sep=""),
                col.names=!file.exists(paste("Results/Classic/RIDE3-Stats-N",ncohorts,".txt",sep="")),
                row.names=FALSE,append=file.exists(paste("Results/Classic/RIDE3-Stats-N",ncohorts,".txt",sep="")))
    
    #mTPI
    mtpi.dat = foreach(m=1:ntrials) %dorng% { simulate.mtpi.trial(dlt.probs[scen,],cohort.size=cohort.size,ncohorts=ncohorts,target=target,delta=0.07,min.esc.size=min.esc.size) }
    mtpi.stats = get.stats(ncohorts,cohort.size,ldat=mtpi.dat,target=target)
    mtpi.ai = round(get.ai(dlt.probs[scen,],target=target,ntrials=ntrials,sel.pct=mtpi.stats[2,]),3)
    write.table(rbind(pdlt=c(NA,dlt.probs[scen,]),mtpi.stats,ai=c(mtpi.ai,rep(NA,length(dlt.probs[scen,])))),sep=",",
                file=paste("Results/Classic/mTPI-Stats-N",ncohorts,".txt",sep=""),
                col.names=!file.exists(paste("Results/Classic/mTPI-Stats-N",ncohorts,".txt",sep="")),
                row.names=FALSE,append=file.exists(paste("Results/Classic/mTPI-Stats-N",ncohorts,".txt",sep="")))
    
    #CRM (with default skeleton)
    crm.dat = foreach(m=1:ntrials) %dorng% { simulate.crm.trial(dlt.probs[scen,],skeleton=crm.skeleton,cohort.size=cohort.size,ncohorts=ncohorts,target=target,min.esc.size=min.esc.size) }
    crm.stats = get.stats(ncohorts,cohort.size,ldat=crm.dat,target=target,crm=TRUE,skeleton=crm.skeleton)
    crm.ai = round(get.ai(dlt.probs[scen,],target=target,ntrials=ntrials,sel.pct=crm.stats[2,]),3)
    write.table(rbind(pdlt=c(NA,dlt.probs[scen,]),crm.stats,ai=c(crm.ai,rep(NA,length(dlt.probs[scen,])))),sep=",",
                file=paste("Results/Classic/CRM-Stats-N",ncohorts,".txt",sep=""),
                col.names=!file.exists(paste("Results/Classic/CRM-Stats-N",ncohorts,".txt",sep="")),
                row.names=FALSE,append=file.exists(paste("Results/Classic/CRM-Stats-N",ncohorts,".txt",sep="")))
    
    #BOIN
    boin.dat = foreach(m=1:ntrials) %dorng% { simulate.boin.trial(dlt.probs[scen,],boin.int=boin.int,cohort.size=cohort.size,ncohorts=ncohorts,target=target,min.esc.size=min.esc.size) }
    boin.stats = get.stats(ncohorts,cohort.size,ldat=boin.dat,target=target)
    boin.ai = round(get.ai(dlt.probs[scen,],target=target,ntrials=ntrials,sel.pct=boin.stats[2,]),3)
    write.table(rbind(pdlt=c(NA,dlt.probs[scen,]),boin.stats,ai=c(boin.ai,rep(NA,length(dlt.probs[scen,])))),sep=",",
                file=paste("Results/Classic/BOIN-Stats-N",ncohorts,".txt",sep=""),
                col.names=!file.exists(paste("Results/Classic/BOIN-Stats-N",ncohorts,".txt",sep="")),
                row.names=FALSE,append=file.exists(paste("Results/Classic/BOIN-Stats-N",ncohorts,".txt",sep="")))
    
    #Print Scen Number
    print(paste("Finished Scenario ",scen,sep=""))
  }
}
stopImplicitCluster()



### Summarize Results
#Load results
ncohorts = 25 #ncohorts = 48
ride1.res = read.table(file=paste("Results/Classic/RIDE1-Stats-N",ncohorts,".txt",sep=""),sep=",",header=TRUE)
ride2.res = read.table(file=paste("Results/Classic/RIDE2-Stats-N",ncohorts,".txt",sep=""),sep=",",header=TRUE)
ride3.res = read.table(file=paste("Results/Classic/RIDE3-Stats-N",ncohorts,".txt",sep=""),sep=",",header=TRUE)
mtpi.res = read.table(file=paste("Results/Classic/mTPI-Stats-N",ncohorts,".txt",sep=""),sep=",",header=TRUE)
crm.res = read.table(file=paste("Results/Classic/CRM-Stats-N",ncohorts,".txt",sep=""),sep=",",header=TRUE)
boin.res = read.table(file=paste("Results/Classic/BOIN-Stats-N",ncohorts,".txt",sep=""),sep=",",header=TRUE)
  
### Create Tables for Manuscript
#Percent Treated At Each Dose
r= 2; res1.ptrt = rbind(sprintf(fmt="%.2f",c(0,dlt.probs[1,])),sprintf(fmt="%.1f",100*ride1.res[r,]),sprintf(fmt="%.1f",100*ride2.res[r,]),sprintf(fmt="%.1f",100*ride3.res[r,]),sprintf(fmt="%.1f",100*mtpi.res[r,]),sprintf(fmt="%.1f",100*crm.res[r,]),sprintf(fmt="%.1f",100*boin.res[r,]))
r= 6; res2.ptrt = rbind(sprintf(fmt="%.2f",c(0,dlt.probs[2,])),sprintf(fmt="%.1f",100*ride1.res[r,]),sprintf(fmt="%.1f",100*ride2.res[r,]),sprintf(fmt="%.1f",100*ride3.res[r,]),sprintf(fmt="%.1f",100*mtpi.res[r,]),sprintf(fmt="%.1f",100*crm.res[r,]),sprintf(fmt="%.1f",100*boin.res[r,]))
r=10; res3.ptrt = rbind(sprintf(fmt="%.2f",c(0,dlt.probs[3,])),sprintf(fmt="%.1f",100*ride1.res[r,]),sprintf(fmt="%.1f",100*ride2.res[r,]),sprintf(fmt="%.1f",100*ride3.res[r,]),sprintf(fmt="%.1f",100*mtpi.res[r,]),sprintf(fmt="%.1f",100*crm.res[r,]),sprintf(fmt="%.1f",100*boin.res[r,]))
r=14; res4.ptrt = rbind(sprintf(fmt="%.2f",c(0,dlt.probs[4,])),sprintf(fmt="%.1f",100*ride1.res[r,]),sprintf(fmt="%.1f",100*ride2.res[r,]),sprintf(fmt="%.1f",100*ride3.res[r,]),sprintf(fmt="%.1f",100*mtpi.res[r,]),sprintf(fmt="%.1f",100*crm.res[r,]),sprintf(fmt="%.1f",100*boin.res[r,]))
r=18; res5.ptrt = rbind(sprintf(fmt="%.2f",c(0,dlt.probs[5,])),sprintf(fmt="%.1f",100*ride1.res[r,]),sprintf(fmt="%.1f",100*ride2.res[r,]),sprintf(fmt="%.1f",100*ride3.res[r,]),sprintf(fmt="%.1f",100*mtpi.res[r,]),sprintf(fmt="%.1f",100*crm.res[r,]),sprintf(fmt="%.1f",100*boin.res[r,]))
colnames(res1.ptrt) = colnames(res2.ptrt) = colnames(res3.ptrt) = colnames(res4.ptrt) = colnames(res5.ptrt) = c("None",1:ncol(dlt.probs))
rownames(res1.ptrt) = rownames(res2.ptrt) = rownames(res3.ptrt) = rownames(res4.ptrt) = rownames(res5.ptrt) = c("DLT Prob","RIDE","RI(-)","RI(+)","mTPI","CRM","BOIN")
tab.ptrt = rbind(res1.ptrt,res2.ptrt,res3.ptrt,res4.ptrt,res5.ptrt)
tab.ptrt

#Percent Selected as MTD + Accuracy Index
r= 3; res1.psel = rbind(sprintf(fmt="%.2f",c(0,dlt.probs[1,],NA)),sprintf(fmt="%.1f",100*unlist(c(ride1.res[r,],ride1.res[r+1,1]))),sprintf(fmt="%.1f",100*unlist(c(ride2.res[r,],ride2.res[r+1,1]))),sprintf(fmt="%.1f",100*unlist(c(ride3.res[r,],ride3.res[r+1,1]))),sprintf(fmt="%.1f",100*unlist(c(mtpi.res[r,],mtpi.res[r+1,1]))),sprintf(fmt="%.1f",100*unlist(c(crm.res[r,],crm.res[r+1,1]))),sprintf(fmt="%.1f",100*unlist(c(boin.res[r,],boin.res[r+1,1]))))
r= 7; res2.psel = rbind(sprintf(fmt="%.2f",c(0,dlt.probs[2,],NA)),sprintf(fmt="%.1f",100*unlist(c(ride1.res[r,],ride1.res[r+1,1]))),sprintf(fmt="%.1f",100*unlist(c(ride2.res[r,],ride2.res[r+1,1]))),sprintf(fmt="%.1f",100*unlist(c(ride3.res[r,],ride3.res[r+1,1]))),sprintf(fmt="%.1f",100*unlist(c(mtpi.res[r,],mtpi.res[r+1,1]))),sprintf(fmt="%.1f",100*unlist(c(crm.res[r,],crm.res[r+1,1]))),sprintf(fmt="%.1f",100*unlist(c(boin.res[r,],boin.res[r+1,1]))))
r=11; res3.psel = rbind(sprintf(fmt="%.2f",c(0,dlt.probs[3,],NA)),sprintf(fmt="%.1f",100*unlist(c(ride1.res[r,],ride1.res[r+1,1]))),sprintf(fmt="%.1f",100*unlist(c(ride2.res[r,],ride2.res[r+1,1]))),sprintf(fmt="%.1f",100*unlist(c(ride3.res[r,],ride3.res[r+1,1]))),sprintf(fmt="%.1f",100*unlist(c(mtpi.res[r,],mtpi.res[r+1,1]))),sprintf(fmt="%.1f",100*unlist(c(crm.res[r,],crm.res[r+1,1]))),sprintf(fmt="%.1f",100*unlist(c(boin.res[r,],boin.res[r+1,1]))))
r=15; res4.psel = rbind(sprintf(fmt="%.2f",c(0,dlt.probs[4,],NA)),sprintf(fmt="%.1f",100*unlist(c(ride1.res[r,],ride1.res[r+1,1]))),sprintf(fmt="%.1f",100*unlist(c(ride2.res[r,],ride2.res[r+1,1]))),sprintf(fmt="%.1f",100*unlist(c(ride3.res[r,],ride3.res[r+1,1]))),sprintf(fmt="%.1f",100*unlist(c(mtpi.res[r,],mtpi.res[r+1,1]))),sprintf(fmt="%.1f",100*unlist(c(crm.res[r,],crm.res[r+1,1]))),sprintf(fmt="%.1f",100*unlist(c(boin.res[r,],boin.res[r+1,1]))))
r=19; res5.psel = rbind(sprintf(fmt="%.2f",c(0,dlt.probs[5,],NA)),sprintf(fmt="%.1f",100*unlist(c(ride1.res[r,],ride1.res[r+1,1]))),sprintf(fmt="%.1f",100*unlist(c(ride2.res[r,],ride2.res[r+1,1]))),sprintf(fmt="%.1f",100*unlist(c(ride3.res[r,],ride3.res[r+1,1]))),sprintf(fmt="%.1f",100*unlist(c(mtpi.res[r,],mtpi.res[r+1,1]))),sprintf(fmt="%.1f",100*unlist(c(crm.res[r,],crm.res[r+1,1]))),sprintf(fmt="%.1f",100*unlist(c(boin.res[r,],boin.res[r+1,1]))))
colnames(res1.psel) = colnames(res2.psel) = colnames(res3.psel) = colnames(res4.psel) = colnames(res5.psel) = c("None",1:ncol(dlt.probs),"Acc Idx")
rownames(res1.psel) = rownames(res2.psel) = rownames(res3.psel) = rownames(res4.psel) = rownames(res5.psel) = c("DLT Prob","RIDE","RI(-)","RI(+)","mTPI","CRM","BOIN")
tab.psel = rbind(res1.psel,res2.psel,res3.psel,res4.psel,res5.psel)
tab.psel

write.table(cbind(tab.ptrt,tab.psel),sep=" & ",quote=FALSE)

### Wages Benchmark and Accuracy Index
# R-code: http://faculty.virginia.edu/model-based_dose-finding/nonparametric%20benchmark.R
# Manuscript: https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5630493/
