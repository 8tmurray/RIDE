### This code validates the implementations of BOIN and CRM
### Written by Thomas A Murray on 11/05/2018
setwd("C:/Users/Thomas/Dropbox/JADE/Software") #Set working directory
source('jade-functions.R') #Load JADE functions
source('comparator-functions.R') #Load comparator functions

#Simulation Parameters
target = 0.3 #Target DLT Probability
delta = 0.07 #Deviation from target that defines overdose
ndose = 6 #Number of investigational dose levels
cohort.size = 3 #Cohort Size
ncohorts = 12 #Maximum number of cohorts
nrep = 100 #Number of times each type of scenario is repeated
ntrials = 2000 #Number of simulated trials for each true scenario

#Comparator Design parameters
crm.skeleton = get.skeleton(ndose,target,delta) #CRM skeleton
boin.int = as.numeric(get.boundary(target,ncohorts,cohort.size)[1:2]) #BOIN dose escalation interval

#Generate <nscens> sets of true DLT probabilties via pseudo-uniform algorithm
set.seed(1985) #This allows us to recover true DLT probabilities in each scenario
scen.type = rep(0:(ndose+1),each=nrep) #Scenario Type
scenarios = lapply(scen.type,function(scen) get.dlt.probs(ndose=ndose,scen=scen,target=target,delta=delta)) #True Randomly Generated Scenarios

#Load results
crm.stats = read.table("Results3/CRM-Stats.txt",sep=",",header=TRUE) %>% as_tibble() %>% arrange(X)
boin.stats = read.table("Results3/BOIN-Stats.txt",sep=",",header=TRUE) %>% as_tibble() %>% arrange(X)

#Validate BOIN
library(BOIN)
iter = 450
#BOIN::get.boundary(target=target,ncohort=ncohorts,cohortsize=cohort.size,p.saf=boin.int[1],p.tox=boin.int[2])
BOIN::get.oc(target=target,p.true=scenarios[[iter]]$dlt.probs,ncohort=ncohorts,cohortsize=cohort.size,p.saf=boin.int[1],p.tox=boin.int[2],ntrial=ntrials)
boin.stats[which(boin.stats$X==iter),]
#Difference is very likely attributable to the change in the early stopping/excessive toxicity rule...makes my BOIN implementation perform better in this case

#Validate CRM
library(dfcrm)
iter = 450
dfcrm::crmsim(target=target,prior=crm.skeleton,PI=scenarios[[iter]]$dlt.probs,n=ncohorts*cohort.size,mcohort=cohort.size,x0=1,scale=sqrt(2),nsim=ntrials)
crm.stats[which(boin.stats$X==iter),]
#Correct selection percentage is within the simulation error bounds (+/- 0.02)
