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
get.stats <- function(ldat,scen,dlt.probs,target,delta,nmax,crm=FALSE,bmacrm=FALSE,skeleton=NULL){

  #Scenario information
  ndose = length(dlt.probs)
  best.dose = ifelse(scen == ndose+1,ndose,scen)
  dlt.probs = c(0,dlt.probs)
  
  #Dose Selection Results
  if(!crm & !bmacrm) psel = prop.table(table(factor(sapply(ldat,function(dat) get.best.dose(dat,target=target)),levels=0:ndose)))
  if(crm) psel = prop.table(table(factor(sapply(ldat,function(dat) get.crm.best.dose(dat,skeleton=skeleton,target=target)),levels=0:ndose)))
  if(bmacrm) psel = prop.table(table(factor(sapply(ldat,function(dat) get.bmacrm.best.dose(dat,skeletons=skeleton,target=target)),levels=0:ndose)))
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
nrep = 100 #Number of times each type of scenario is repeated
ntrials = 2000 #Number of simulated trials for each true scenario

#Parallely loop over every true scenario
library(foreach); library(doParallel); library(doRNG)
RNGkind("L'Ecuyer"); registerDoParallel(cores=12)
#Loop over targets
for(target in targets){

  #Comparator Design parameters
  crm.skeleton = get.skeleton(ndose,target,delta=0.25*target) #CRM skeleton
  bmacrm.skeletons = rbind(crm.skeleton,get.skeleton(ndose,target,prior.best=2,delta=0.4*target),get.skeleton(ndose,target,prior.best=ndose,delta=0.15*target))
  boin.int = as.numeric(BOIN::get.boundary(target,1,3)[1:2]) #BOIN dose escalation interval
  
  #Generate <nscens> sets of true DLT probabilties via pseudo-uniform algorithm
  set.seed(1985) #This allows us to recover true DLT probabilities in each scenario
  scen.type = rep(0:(ndose+1),each=nrep) #Scenario Type
  scenarios = lapply(scen.type,function(scen) get.dlt.probs(ndose=ndose,target=target,delta=0.25*target,scen=scen)) #True Randomly Generated Scenarios

  #Loop over ncohorts.vec  
  for(ncohorts in ncohorts.vec){
    
    #Loop over scenarios
    for(iter in 1:length(scen.type)){

      ###True DLT Probabilities
      scen = scenarios[[iter]]$scen
      dlt.probs = scenarios[[iter]]$dlt.probs
      
      ###Simulate <ntrials> trials for each design and get summary statistics
      #RIDE (with default working value OR=2)
      ride.dat = foreach(m=1:ntrials) %dorng% { simulate.trial(dlt.probs,cohort.size=cohort.size,ncohorts=ncohorts,target=target,or=2,min.esc.size=3) }
      ride.stats = get.stats(ride.dat,scen,dlt.probs,nmax=ncohorts*cohort.size,target=target,delta=delta)
      write.table(cbind(rbind(iter),round(ride.stats,4)),sep=",",row.names=FALSE,
                  file=paste("Results/Contemporary/RIDE-Stats-N",ncohorts*cohort.size,"-P",100*target,".txt",sep=""),
                  col.names=!file.exists(paste("Results/Contemporary/RIDE-Stats-N",ncohorts*cohort.size,"-P",100*target,".txt",sep="")),
                  append=file.exists(paste("Results/Contemporary/RIDE-Stats-N",ncohorts*cohort.size,"-P",100*target,".txt",sep="")))
      
      #mTPI
      mtpi.dat = foreach(m=1:ntrials) %dorng% { simulate.mtpi.trial(dlt.probs,cohort.size=cohort.size,ncohorts=ncohorts,target=target,delta=0.25*target,min.esc.size=3) }
      mtpi.stats = get.stats(mtpi.dat,scen,dlt.probs,nmax=ncohorts*cohort.size,target=target,delta=delta)
      write.table(cbind(rbind(iter),round(mtpi.stats,4)),sep=",",row.names=FALSE,
                  file=paste("Results/Contemporary/mTPI-Stats-N",ncohorts*cohort.size,"-P",100*target,".txt",sep=""),
                  col.names=!file.exists(paste("Results/Contemporary/mTPI-Stats-N",ncohorts*cohort.size,"-P",100*target,".txt",sep="")),
                  append=file.exists(paste("Results/Contemporary/mTPI-Stats-N",ncohorts*cohort.size,"-P",100*target,".txt",sep="")))
      
      #Keyboard
      key.dat = foreach(m=1:ntrials) %dorng% { simulate.keyboard.trial(dlt.probs,cohort.size=cohort.size,ncohorts=ncohorts,target=target,delta=0.25*target,min.esc.size=3) }
      key.stats = get.stats(key.dat,scen,dlt.probs,nmax=ncohorts*cohort.size,target=target,delta=delta)
      write.table(cbind(rbind(iter),round(key.stats,4)),sep=",",row.names=FALSE,
                  file=paste("Results/Contemporary/Key-Stats-N",ncohorts*cohort.size,"-P",100*target,".txt",sep=""),
                  col.names=!file.exists(paste("Results/Contemporary/Key-Stats-N",ncohorts*cohort.size,"-P",100*target,".txt",sep="")),
                  append=file.exists(paste("Results/Contemporary/Key-Stats-N",ncohorts*cohort.size,"-P",100*target,".txt",sep="")))
      
      #CRM (with default skeleton)
      crm.dat = foreach(m=1:ntrials) %dorng% { simulate.crm.trial(dlt.probs,skeleton=crm.skeleton,cohort.size=cohort.size,ncohorts=ncohorts,target=target,min.esc.size=3) }
      crm.stats = get.stats(crm.dat,scen,dlt.probs,nmax=ncohorts*cohort.size,target=target,delta=delta,crm=TRUE,skeleton=crm.skeleton)
      write.table(cbind(rbind(iter),round(crm.stats,4)),sep=",",row.names=FALSE,
                  file=paste("Results/Contemporary/CRM-Stats-N",ncohorts*cohort.size,"-P",100*target,".txt",sep=""),
                  col.names=!file.exists(paste("Results/Contemporary/CRM-Stats-N",ncohorts*cohort.size,"-P",100*target,".txt",sep="")),
                  append=file.exists(paste("Results/Contemporary/CRM-Stats-N",ncohorts*cohort.size,"-P",100*target,".txt",sep="")))
      
      #BMA-CRM (with default skeletons)
      bmacrm.dat = foreach(m=1:ntrials) %dorng% { simulate.bmacrm.trial(dlt.probs,skeletons=bmacrm.skeletons,cohort.size=cohort.size,ncohorts=ncohorts,target=target,min.esc.size=3) }
      bmacrm.stats = get.stats(bmacrm.dat,scen,dlt.probs,nmax=ncohorts*cohort.size,target=target,delta=delta,bmacrm=TRUE,skeleton=bmacrm.skeletons)
      write.table(cbind(rbind(iter),round(bmacrm.stats,4)),sep=",",row.names=FALSE,
                  file=paste("Results/Contemporary/BMACRM-Stats-N",ncohorts*cohort.size,"-P",100*target,".txt",sep=""),
                  col.names=!file.exists(paste("Results/Contemporary/BMACRM-Stats-N",ncohorts*cohort.size,"-P",100*target,".txt",sep="")),
                  append=file.exists(paste("Results/Contemporary/BMACRM-Stats-N",ncohorts*cohort.size,"-P",100*target,".txt",sep="")))
      
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
library(extrafont) #For eps figures
#Plot summarizing scenarios
set.seed(1985) #This allows us to recover true DLT probabilities in each scenario
scen.type = rep(0:(ndose+1),each=100) #Scenario Type
scenarios = lapply(scen.type,function(scen) get.dlt.probs(ndose=ndose,target=0.30,delta=0.07,scen=scen)) #True Randomly Generated Scenarios
scens = as.vector(sapply(0:7*100,function(x) sample(x:(x+99)+1,2)))

postscript("Scenario-Boxplot.eps", width = 4, height = 4, horizontal = FALSE, family="serif")
#pdf(file="Scenario-Boxplot.pdf")
par(lab=c(6,10,7),family="serif")
boxplot(t(sapply(1:800,function(i) scenarios[[i]]$dlt.probs)),las=1,ylab=expression(paste("DLT Probability ( ",pi[j]," )"),sep=""),xlab="Dose Level (j)",main="Scenario Coverage",pch=20)
for(i in 1:16) lines(1:6,scenarios[[scens[i]]]$dlt.probs,type='b',pch=20,col="grey")
abline(h=0.30,col="red",lty=3,lwd=2)
dev.off()

##############
#Figure (2 x 3) of overall sB and/or AI by target and ncohorts
targets = c(0.1,0.2,0.3); ncohorts.vec = c(8,12); cohort.size = 3 
stat.vec = c("sBD","sAD","sOD","sUD","sAI","mDLT","tBD","tAD","tOD","tUD","tAI","pOD40","pOD70","pUD40","pUD70")
ylab.vec = c("Proportion Selecting Best Dose",
            "Proportion Selecting Adequate Dose",
            "Proportion Selecting Overdose",
            "Proportion Selecting Underdose",
            "Selection Accuracy Index",
            "Average # of DLTs",
            "Average Percent Given Best Dose",
            "Average Percent Given Adequate Dose",
            "Average Percent Given Overdose",
            "Average Percent Given Underdose",
            "Treatment Accuracy Index",
            "Proportion Giving >40% Overdose",
            "Proportion Giving >70% Overdose",
            "Proportion Giving >40% Underdose",
            "Proportion Giving >70% Underdose")
for(s in 1:length(stat.vec)){
  stat = stat.vec[s]; yl = ylab.vec[s]
  postscript(paste("Contemporary-",stat,".eps",sep=""),family="serif",height=4,width=6,horizontal=FALSE)
  #pdf(paste("Contemporary-",stat,".pdf",sep=""),family="serif",height=4,width=6)
  par(mfcol=c(2,3),mar=c(5,4,3,0)+0.1,lab=c(6,10,7),family="serif")
  for(target in targets){
    for(ncohorts in ncohorts.vec){
      #Load Data
      ride.stats = read.table(paste("Results/Contemporary/RIDE-Stats-N",ncohorts*cohort.size,"-P",100*target,".txt",sep=""),sep=",",header=TRUE)
      mtpi.stats = read.table(paste("Results/Contemporary/mTPI-Stats-N",ncohorts*cohort.size,"-P",100*target,".txt",sep=""),sep=",",header=TRUE)
      key.stats = read.table(paste("Results/Contemporary/Key-Stats-N",ncohorts*cohort.size,"-P",100*target,".txt",sep=""),sep=",",header=TRUE)
      boin.stats = read.table(paste("Results/Contemporary/BOIN-Stats-N",ncohorts*cohort.size,"-P",100*target,".txt",sep=""),sep=",",header=TRUE)
      crm.stats = read.table(paste("Results/Contemporary/CRM-Stats-N",ncohorts*cohort.size,"-P",100*target,".txt",sep=""),sep=",",header=TRUE)
      bmacrm.stats = read.table(paste("Results/Contemporary/BMACRM-Stats-N",ncohorts*cohort.size,"-P",100*target,".txt",sep=""),sep=",",header=TRUE)
      
      #Compile relevant statistic
      res = cbind(loRIDE = ride.stats[,stat], mTPI = mtpi.stats[,stat], Keyboard = key.stats[,stat], BOIN = boin.stats[,stat], CRM = crm.stats[,stat], BMACRM = bmacrm.stats[,stat])
      if(target==0.1) boxplot(res,las=2,pch=20,range=0,ylab=yl,main=bquote(expression(paste(pi,"* = ",.(target),", n =",.(ncohorts),"x",.(cohort.size),sep=""))))
      if(target!=0.1) boxplot(res,las=2,pch=20,range=0,ylab=" ",main=bquote(expression(paste(pi,"* = ",.(target),", n =",.(ncohorts),"x",.(cohort.size),sep=""))))
    }
  }
  dev.off()
}

#Supplementary Comparative Tables
for(target in targets){
  for(ncohorts in ncohorts.vec){
    #Load Data
    ride.stats = read.table(paste("Results/Contemporary/RIDE-Stats-N",ncohorts*cohort.size,"-P",100*target,".txt",sep=""),sep=",",header=TRUE)
    mtpi.stats = read.table(paste("Results/Contemporary/mTPI-Stats-N",ncohorts*cohort.size,"-P",100*target,".txt",sep=""),sep=",",header=TRUE)
    key.stats = read.table(paste("Results/Contemporary/Key-Stats-N",ncohorts*cohort.size,"-P",100*target,".txt",sep=""),sep=",",header=TRUE)
    boin.stats = read.table(paste("Results/Contemporary/BOIN-Stats-N",ncohorts*cohort.size,"-P",100*target,".txt",sep=""),sep=",",header=TRUE)
    crm.stats = read.table(paste("Results/Contemporary/CRM-Stats-N",ncohorts*cohort.size,"-P",100*target,".txt",sep=""),sep=",",header=TRUE)
    bmacrm.stats = read.table(paste("Results/Contemporary/BMACRM-Stats-N",ncohorts*cohort.size,"-P",100*target,".txt",sep=""),sep=",",header=TRUE)

      #Averaging across all scenarios
    tab <- rbind(colMeans(ride.stats[,-1],na.rm = TRUE),
                 colMeans(mtpi.stats[,-1],na.rm = TRUE),
                 colMeans(key.stats[,-1],na.rm = TRUE),
                 colMeans(boin.stats[,-1],na.rm = TRUE),
                 colMeans(crm.stats[,-1],na.rm = TRUE),
                 colMeans(bmacrm.stats[,-1],na.rm = TRUE))
    rownames(tab) = c("loRIDE","mTPI","Keyboard","BOIN","CRM","BMA-CRM")
    colnames(tab) = c("Best Dose (sBD)","Adequate Dose (sAD)","Overdose (sOD)","Underdose (sUD)","Accuracy Index (sAI)",
                      "Average DLTs (mDLTs)","Best Dose (tBD)","Adequate Dose (tAD)","Overdose (tOD)","Underdose (tUD)","Accuracy Index (sAI)",
                      "40% Overdose","70% Overdose","40% Underdose","70% Underdose")
    write.table(round(t(tab),3),sep=" & ",quote=F)
    print(" \n ")
  }
}




###################
#Figure (3 x 3) overall & by type sB and/or AI, barchart?
#Break-out by scenario type for a particular target and sample size
targets = c(0.1,0.2,0.3); ncohorts.vec = c(8,12); cohort.size = 3
stat = "sBD"; yl = "Proportion Selecting Best Dose"
#stat = "sAI"; yl = "Selection Accuracy Index"

for(target in targets){
  for(ncohorts in ncohorts.vec){
    #Load Data
    ride.stats = read.table(paste("Results/Contemporary/RIDE-Stats-N",ncohorts*cohort.size,"-P",100*target,".txt",sep=""),sep=",",header=TRUE)
    mtpi.stats = read.table(paste("Results/Contemporary/mTPI-Stats-N",ncohorts*cohort.size,"-P",100*target,".txt",sep=""),sep=",",header=TRUE)
    key.stats = read.table(paste("Results/Contemporary/Key-Stats-N",ncohorts*cohort.size,"-P",100*target,".txt",sep=""),sep=",",header=TRUE)
    boin.stats = read.table(paste("Results/Contemporary/BOIN-Stats-N",ncohorts*cohort.size,"-P",100*target,".txt",sep=""),sep=",",header=TRUE)
    crm.stats = read.table(paste("Results/Contemporary/CRM-Stats-N",ncohorts*cohort.size,"-P",100*target,".txt",sep=""),sep=",",header=TRUE)
    bmacrm.stats = read.table(paste("Results/Contemporary/BMACRM-Stats-N",ncohorts*cohort.size,"-P",100*target,".txt",sep=""),sep=",",header=TRUE)
    
    postscript(paste(stat,"-by-Type-N",ncohorts*cohort.size,"-p",round(100*target),".eps",sep=""),family="serif",height=6,width=6,horizontal=FALSE)
    #pdf(paste(stat,"-by-Type-N",ncohorts*cohort.size,"-p",round(100*target),".pdf",sep=""),family="serif",height=6,width=6)
    par(mfrow=c(3,3),mar=c(5,4,2,0)+0.1,lab=c(6,10,7),family="serif")
    
    #All scenarios
    scens = 1:800
    
    #Compile relevant statistic
    res = cbind(loRIDE = ride.stats[scens,stat], mTPI = mtpi.stats[scens,stat], Keyboard = key.stats[scens,stat], BOIN = boin.stats[scens,stat], CRM = crm.stats[scens,stat], BMACRM = bmacrm.stats[scens,stat])
    boxplot(res,las=2,pch=20,range=0,ylab=yl,main="Overall")
    
    #write out table
    tab <- rbind(colMeans(ride.stats[scens,-1],na.rm = TRUE),
                 colMeans(mtpi.stats[scens,-1],na.rm = TRUE),
                 colMeans(key.stats[scens,-1],na.rm = TRUE),
                 colMeans(boin.stats[scens,-1],na.rm = TRUE),
                 colMeans(crm.stats[scens,-1],na.rm = TRUE),
                 colMeans(bmacrm.stats[scens,-1],na.rm = TRUE))
    rownames(tab) = c("loRIDE","mTPI","Keyboard","BOIN","CRM","BMA-CRM")
    write.table(round(t(tab),3),sep=" & ",quote=F)
    print("\n")
    
    #By scenario type
    for(type in 0:7){
      scens = 1:100 + 100*type
      
      #Compile relevant statistic and Plot
      res = cbind(loRIDE = ride.stats[scens,stat], mTPI = mtpi.stats[scens,stat], Keyboard = key.stats[scens,stat], BOIN = boin.stats[scens,stat], CRM = crm.stats[scens,stat], BMACRM = bmacrm.stats[scens,stat])
      boxplot(res,las=2,ylab=yl,pch=20,main=paste("Scenario Type ",type,sep=""))
      
      #write out table
      tab <- rbind(colMeans(ride.stats[scens,-1],na.rm = TRUE),
                   colMeans(mtpi.stats[scens,-1],na.rm = TRUE),
                   colMeans(key.stats[scens,-1],na.rm = TRUE),
                   colMeans(boin.stats[scens,-1],na.rm = TRUE),
                   colMeans(crm.stats[scens,-1],na.rm = TRUE),
                   colMeans(bmacrm.stats[scens,-1],na.rm = TRUE))
      rownames(tab) = c("loRIDE","mTPI","Keyboard","BOIN","CRM","BMA-CRM")
      write.table(round(t(tab),3),sep=" & ",quote=F)
      print("\n")  
    }
    dev.off()
  }
}





##############
#Comparing versions of loRIDE
#Figure (2 x 3) of overall sB and/or AI by target and ncohorts
targets = c(0.1,0.2,0.3); ncohorts.vec = c(8,12); cohort.size = 3; ors = c(1.2,1.6,2,2.5,3.4) 
stat.vec = c("sBD","sAD","sOD","sUD","sAI","mDLT","tBD","tAD","tOD","tUD","tAI","pOD40","pOD70","pUD40","pUD70")
ylab.vec = c("Proportion Selecting Best Dose",
             "Proportion Selecting Adequate Dose",
             "Proportion Selecting Overdose",
             "Proportion Selecting Underdose",
             "Selection Accuracy Index",
             "Average # of DLTs",
             "Average Percent Given Best Dose",
             "Average Percent Given Adequate Dose",
             "Average Percent Given Overdose",
             "Average Percent Given Underdose",
             "Treatment Accuracy Index",
             "Proportion Giving >40% Overdose",
             "Proportion Giving >70% Overdose",
             "Proportion Giving >40% Underdose",
             "Proportion Giving >70% Underdose")
for(s in 1:length(stat.vec)){
  stat = stat.vec[s]; yl = ylab.vec[s]
  postscript(paste("loRIDE-",stat,".eps",sep=""),family="serif",height=4,width=6,horizontal=FALSE)
  #pdf(paste("loRIDE-",stat,".pdf",sep=""),family="serif",height=4,width=6)
  par(mfcol=c(2,3),mar=c(5,4,3,0)+0.1,lab=c(6,10,7),family="serif")
  for(target in targets){
    for(ncohorts in ncohorts.vec){
      #Load Data
      ride1.stats = read.table(paste("Results/Contemporary/RIDE-Stats-OR1.2-N",ncohorts*cohort.size,"-P",100*target,".txt",sep=""),sep=",",header=TRUE)
      ride2.stats = read.table(paste("Results/Contemporary/RIDE-Stats-OR1.6-N",ncohorts*cohort.size,"-P",100*target,".txt",sep=""),sep=",",header=TRUE)
      ride3.stats = read.table(paste("Results/Contemporary/RIDE-Stats-OR2-N",ncohorts*cohort.size,"-P",100*target,".txt",sep=""),sep=",",header=TRUE)
      ride4.stats = read.table(paste("Results/Contemporary/RIDE-Stats-OR2.5-N",ncohorts*cohort.size,"-P",100*target,".txt",sep=""),sep=",",header=TRUE)
      ride5.stats = read.table(paste("Results/Contemporary/RIDE-Stats-OR3.4-N",ncohorts*cohort.size,"-P",100*target,".txt",sep=""),sep=",",header=TRUE)
      
      #Compile relevant statistic
      res = cbind(loRIDE12 = ride1.stats[,stat], loRIDE16 = ride2.stats[,stat], loRIDE20 = ride3.stats[,stat], loRIDE25 = ride4.stats[,stat], loRIDE35 = ride5.stats[,stat])
      if(target==0.1) boxplot(res,las=2,pch=20,range=0,ylab=yl,main=bquote(expression(paste(pi,"* = ",.(target),", n =",.(ncohorts),"x",.(cohort.size),sep=""))))
      if(target!=0.1) boxplot(res,las=2,pch=20,range=0,ylab=" ",main=bquote(expression(paste(pi,"* = ",.(target),", n =",.(ncohorts),"x",.(cohort.size),sep=""))))
    }
  }
  dev.off()
}


#Supplementary Comparative Tables
for(target in targets){
  for(ncohorts in ncohorts.vec){
    #Load Data
    ride1.stats = read.table(paste("Results/Contemporary/RIDE-Stats-OR1.2-N",ncohorts*cohort.size,"-P",100*target,".txt",sep=""),sep=",",header=TRUE)
    ride2.stats = read.table(paste("Results/Contemporary/RIDE-Stats-OR1.6-N",ncohorts*cohort.size,"-P",100*target,".txt",sep=""),sep=",",header=TRUE)
    ride3.stats = read.table(paste("Results/Contemporary/RIDE-Stats-OR2-N",ncohorts*cohort.size,"-P",100*target,".txt",sep=""),sep=",",header=TRUE)
    ride4.stats = read.table(paste("Results/Contemporary/RIDE-Stats-OR2.5-N",ncohorts*cohort.size,"-P",100*target,".txt",sep=""),sep=",",header=TRUE)
    ride5.stats = read.table(paste("Results/Contemporary/RIDE-Stats-OR3.4-N",ncohorts*cohort.size,"-P",100*target,".txt",sep=""),sep=",",header=TRUE)
    
    #Averaging across all scenarios
    tab <- rbind(colMeans(ride1.stats[,-1],na.rm = TRUE),
                 colMeans(ride2.stats[,-1],na.rm = TRUE),
                 colMeans(ride3.stats[,-1],na.rm = TRUE),
                 colMeans(ride4.stats[,-1],na.rm = TRUE),
                 colMeans(ride5.stats[,-1],na.rm = TRUE))
    rownames(tab) = c("loRIDE12","loRIDE16","loRIDE20","loRIDE25","loRIDE34")
    write.table(round(t(tab),3),sep=" & ",quote=F)
    print(" \n ")
  }
}



###################
#Figure (3 x 3) overall & by type sB and/or AI, barchart?
#Break-out by scenario type for a particular target and sample size
target = 0.3; ncohorts = 12; cohort.size = 3

#Load Data
ride1.stats = read.table(paste("Results/Contemporary/RIDE-Stats-OR1.2-N",ncohorts*cohort.size,"-P",100*target,".txt",sep=""),sep=",",header=TRUE)
ride2.stats = read.table(paste("Results/Contemporary/RIDE-Stats-OR1.6-N",ncohorts*cohort.size,"-P",100*target,".txt",sep=""),sep=",",header=TRUE)
ride3.stats = read.table(paste("Results/Contemporary/RIDE-Stats-OR2-N",ncohorts*cohort.size,"-P",100*target,".txt",sep=""),sep=",",header=TRUE)
ride4.stats = read.table(paste("Results/Contemporary/RIDE-Stats-OR2.5-N",ncohorts*cohort.size,"-P",100*target,".txt",sep=""),sep=",",header=TRUE)
ride5.stats = read.table(paste("Results/Contemporary/RIDE-Stats-OR3.4-N",ncohorts*cohort.size,"-P",100*target,".txt",sep=""),sep=",",header=TRUE)

stat = "sBD"; yl = "Proportion Selecting Best Dose"
#stat = "sAI"; yl = "Selection Accuracy Index"
postscript(paste("loRIDE-",stat,"-by-Type.eps",sep=""),family="serif",height=6,width=6,horizontal=FALSE)
#pdf(paste("loRIDE-",stat,"-by-Type.eps",sep=""),family="serif",height=6,width=6)
par(mfrow=c(3,3),mar=c(5,4,2,0)+0.1,lab=c(6,10,7),family="serif")

#All scenarios
scens = 1:800

#Compile relevant statistic
res = cbind(loRIDE12 = ride1.stats[scens,stat], loRIDE16 = ride2.stats[scens,stat], loRIDE20 = ride3.stats[scens,stat], loRIDE25 = ride4.stats[scens,stat], loRIDE35 = ride5.stats[scens,stat])
boxplot(res,las=2,pch=20,range=0,ylab=yl,main="Overall")

#write out table
tab <- rbind(colMeans(ride1.stats[scens,-1],na.rm = TRUE),
             colMeans(ride2.stats[scens,-1],na.rm = TRUE),
             colMeans(ride3.stats[scens,-1],na.rm = TRUE),
             colMeans(ride4.stats[scens,-1],na.rm = TRUE),
             colMeans(ride5.stats[scens,-1],na.rm = TRUE))
rownames(tab) = rownames(tab) = c("loRIDE12","loRIDE16","loRIDE20","loRIDE25","loRIDE34")
write.table(round(t(tab),3),sep=" & ",quote=F)
print("\n")

#By scenario type
for(type in 0:7){
  scens = 1:100 + 100*type
  
  #Compile relevant statistic and Plot
  res = cbind(loRIDE12 = ride1.stats[scens,stat], loRIDE16 = ride2.stats[scens,stat], loRIDE20 = ride3.stats[scens,stat], loRIDE25 = ride4.stats[scens,stat], loRIDE35 = ride5.stats[scens,stat])
  boxplot(res,las=2,ylab=yl,pch=20,main=paste("Scenario Type ",type,sep=""))
  
  #write out table
  tab <- rbind(colMeans(ride1.stats[scens,-1],na.rm = TRUE),
               colMeans(ride2.stats[scens,-1],na.rm = TRUE),
               colMeans(ride3.stats[scens,-1],na.rm = TRUE),
               colMeans(ride4.stats[scens,-1],na.rm = TRUE),
               colMeans(ride5.stats[scens,-1],na.rm = TRUE))
  rownames(tab) = rownames(tab) = c("loRIDE12","loRIDE16","loRIDE20","loRIDE25","loRIDE34")
  write.table(round(t(tab),3),sep=" & ",quote=F)
  print("\n")  
}
dev.off()
