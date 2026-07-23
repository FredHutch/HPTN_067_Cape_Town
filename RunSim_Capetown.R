#####Things to Change
filename = paste0("data/Capetown.067.",Sys.Date(),".rda")
#Adherence regime switching
adherence.target= 0    #0 default, >0 to trigger adherence targeted regime switching
#067 participants only or 10k simulated participants
participants = "10k" #"067" or "10k"
#211 vs daily PrEP or all regimens
regimens = "211" #"211" or "all"

#Packages
library(plyr)
library(lme4)
library(scales)
library(zoo)
library(R.matlab)
library(expm)

#Load functions 
data_location = "~/PrEPSIM/data/" 
setwd("~/PrEPSIM/code")
source('Load067_ciswomen.R')
source('HPTN067_GLMM_ciswomen.R')
source('Efficacy_Calculation_ciswomen.R')
source('Crossover_Simulation_Code.R')

#Efficacy values  
#Moore efficacy 0.593 0.838 0.959 Moore M, et al. Dosing forgiveness of oral PrEP for cisgender women remains uncertain. J Int AIDS Soc. 2025 May;28(5):e26496. doi: 10.1002/jia2.26496. PMID: 40384398; PMCID: PMC12086362.
#Moore minimum efficacy 0.299 0.517 0.726 (same as above)
#Marazzo efficacy 61.61943  89.75357 100 Marrazzo J, et al. HIV Preexposure Prophylaxis With Emtricitabine and Tenofovir Disoproxil Fumarate Among Cisgender Women. Jama. 2024;331(11):930-7.
maxefficacy=0.959     #7 pills/week 
fullefficacy=0.838    #4 pills/week 
partialefficacy=0.593 #2 pills/week

condomefficacy=.90 #Foss, A.M., et al., A systematic review of published evidence on intervention impact on condom use in sub-Saharan Africa and Asia. Sexually Transmitted Infections, 2007. 83(7): p. 510-516.

eventefficacy  = 0.838 #fullefficacy, 0.96 with perfect 211 adherence (sensitivity analysis)
r211efficacy   = 0.838 #1 #not calculated without special case, 0.96 with perfect 211 adherence (sensitivity analysis)

#Crossover sim values
pillratefactor.eff=0.01  #0.01 to include all individuals with increased efficacy
efficacy.diff.eff =0.1   #0.1 default, 0 to include all individuals with increased efficacy

#Subset
#load("~/PrEPSIM/PK_validation_analysis/ciswomen_conc_mismatch_ids.rda")
mismatch_subset = c(0,0)       #conc_pill_mismatch_ids: remove all ids with any mismatch in pill report vs PK; conc0_ids: remove all ids with reported 2+ pills & 0 drug concentration; c(0,0): no subsetting

pillfile.capetown=readMat(paste0(data_location, "data_v1/people_data_C.mat")) 

#Participant Data
participant.dat=read.csv(paste0(data_location, "data_v1/participants_B.csv"))
participant.dat.C=subset(participant.dat,site=="Cape Town|South Africa")
conc_mismatch_location = which(participant.dat.C$pub_id %in% mismatch_subset)
#conc_mismatch_location = which(participant.dat.C$pub_id %in% conc_pill_mismatch_ids)
participant.dat.C1=subset(participant.dat.C,pub_id %notin% mismatch_subset)

pd.C=processdata2(pillfile.capetown,participant.dat.C1)
dat.capetown=pd.C[[1]]
sexactdata.capetown=pd.C[[2]]
participant.df.capetown=pd.C[[3]]
participantstats.capetown=pd.C[[4]]
dailysex.capetown=pd.C[[5]]

ARM.capetown=participantstats.capetown$ARM

for(i in seq_along(ARM.capetown)){
  dailysex.capetown[[i]]$id=i
  if(nrow(sexactdata.capetown[[i]])>0){
    sexactdata.capetown[[i]]$id=i
  }
  sex.day.of=(dailysex.capetown[[i]]$sexonly+dailysex.capetown[[i]]$combopre+dailysex.capetown[[i]]$combopost+dailysex.capetown[[i]]$comboprepost)>0
  nopre.day.of=(dailysex.capetown[[i]]$sexonly+dailysex.capetown[[i]]$combopost)>0 
  nopost.day.of=(dailysex.capetown[[i]]$sexonly+dailysex.capetown[[i]]$combopre)>0
  N.days=length(sex.day.of)
  analsex.day.of=dailysex.capetown[[i]]$sexonly.a+dailysex.capetown[[i]]$combopre.a+dailysex.capetown[[i]]$combopost.a+dailysex.capetown[[i]]$comboprepost.a
  dailysex.capetown[[i]]$sexnextday=NA
  dailysex.capetown[[i]]$analsexnextday=NA
  dailysex.capetown[[i]]$sexnextday=c(sex.day.of[-1],NA)
  dailysex.capetown[[i]]$analsexnextday=c(analsex.day.of[-1],NA)
  dailysex.capetown[[i]]$sexdayof=sex.day.of
  dailysex.capetown[[i]]$sexdaybefore=c(NA,sex.day.of[-N.days])
  dailysex.capetown[[i]]$nopredayof=nopre.day.of
  dailysex.capetown[[i]]$noprenextday=c(nopre.day.of[-1],NA)
  dailysex.capetown[[i]]$nopostdayof=nopost.day.of
  dailysex.capetown[[i]]$nopostdaybefore=c(NA,nopost.day.of[-N.days])
  dailysex.capetown[[i]]$analsexdayof=analsex.day.of
  dailysex.capetown[[i]]$eventpredict=as.numeric(sex.day.of|dailysex.capetown[[i]]$sexnextday)
  dailysex.capetown[[i]]$nosex=1-as.numeric(sex.day.of|dailysex.capetown[[i]]$sexnextday)
  dailysex.capetown[[i]]$eventpredictanal=as.numeric(analsex.day.of>0|dailysex.capetown[[i]]$analsexnextday>0)
  dailysex.capetown[[i]]$pilldayof=dailysex.capetown[[i]]$comboprepost+dailysex.capetown[[i]]$pillonly+dailysex.capetown[[i]]$combopost+dailysex.capetown[[i]]$combopre>0
  dailysex.capetown[[i]]$pilldaybefore=c(NA,dailysex.capetown[[i]]$pilldayof[-N.days])
  dailysex.capetown[[i]]$pillweekbefore=NA
  dailysex.capetown[[i]]$pilllasttwo=NA
  dailysex.capetown[[i]]$pilllastthree=NA
  if(N.days>7){
    dailysex.capetown[[i]]$pillweekbefore=c(rep(NA,7),dailysex.capetown[[i]]$pillonly[-seq(N.days-6,N.days)]) #not sure what this is capturing? 
    y=cumsum(dailysex.capetown[[i]]$pillonly+dailysex.capetown[[i]]$combopre+dailysex.capetown[[i]]$combopost+dailysex.capetown[[i]]$comboprepost)
    dailysex.capetown[[i]]$pilllasttwo=c(0,y[1:(N.days-1)])-c(0,0,0,y[1:(N.days-3)])
    dailysex.capetown[[i]]$pilllastthree=c(0,y[1:(N.days-1)])-c(0,0,0,0,y[1:(N.days-4)])
  }
}

ad.daily=sapply(dailysex.capetown[ARM.capetown==1],getadherence.daily)
ad.time=sapply(dailysex.capetown[ARM.capetown==2],getadherence.time)
ad.event=sapply(dailysex.capetown[ARM.capetown==3],getadherence.event)
ad.211=sapply(dailysex.capetown[ARM.capetown==3],getadherence.211)

spec.event=sapply(dailysex.capetown[ARM.capetown==3], evaluate.specificity)
sens.event=sapply(dailysex.capetown[ARM.capetown==3], evaluate.sensitivity)

spec.211=sapply(dailysex.capetown[ARM.capetown==3], evaluate.specificity)
sens.211=sapply(dailysex.capetown[ARM.capetown==3], evaluate.sensitivity)

participantstats.capetown$adherence.regimen=c(ad.daily,ad.event,ad.time)#,ad.211)

eff.daily = get.efficacy.trial(subset(participantstats.capetown, ARM==1))
eff.time = get.efficacy.trial(subset(participantstats.capetown, ARM==2))
eff.event = get.efficacy.trial(subset(participantstats.capetown, ARM==3))
eff.211 = get.efficacy.trial(subset(participantstats.capetown, ARM==3))

################### Construct data for analysis #line 488 
data=construct.allsex.capetown() #in HPTN067GLMMR
model.list=fit.models.all(data)  #in HPTN067GLMMR

################### Simulations
if(participants=="10k"){ 
  ##Generate individual parameters (both for the true 067 individuals and 10000 “fake” individuals
  params.rand=get.random.params(10000,model.list) #in HPTN067GLMMR

  ##Calculate the efficacy associated with those parameters (via get.efficacy). 
  #All trial arms are simulated separately.
  efficacy.rand=get.efficacy(params.rand) #in Crossoversimulationcode

  #Perform cross over simulations of both the 067 and random individuals. This is done by linking up subject behavior 
  #in the different arms and then assuming they are the same person. This is done according to various different rules 
  #(using reassign.efficacy.by.class).
  efficacy.assign.by.quintile.list=reassign.efficacy.byclass(efficacy.rand, params.rand) #in Crossoversimulationcode

  #efficacy.randomize.list=reassign.efficacy.byclass(efficacy.rand, params.rand, N.bins = 1)
  #efficacy.assign.by.quintile.separate.list=reassign.efficacy.byclass(efficacy.rand, params.rand, attribute.time = "adherence.regular", attribute.event = "adherence.regular")
  #efficacy.assign.by.quintile.twostage.list=reassign.efficacy.byclass(efficacy.rand, params.rand, attribute.time = "adherence.regular", attribute.event = "adherence.sex", linkevent.to.time = TRUE)

} else if(participants =="067"){
  ##Generate individual parameters (both for the true 067 individuals and 10000 “fake” individuals
  params.067=get.HPTN067.params(model.list)      #in HPTN067GLMMR
  ##Calculate the efficacy associated with those parameters (via get.efficacy). 
  #All trial arms are simulated separately.
  efficacy.067=get.efficacy(params.067)  #in Crossoversimulationcode
  #Perform cross over simulations of both the 067 and random individuals. This is done by linking up subject behavior 
  #in the different arms and then assuming they are the same person. This is done according to various different rules 
  #(using reassign.efficacy.by.class).
  efficacy.assign.by.quintile.list=reassign.efficacy.byclass(efficacy.067, params.067)
}

################## 
if(regimens == "211"){ 
  efficacy.twoway.211=process.efficacy(efficacy.assign.by.quintile.list[[1]], efficacy.diff1 = 0.1, event=0, r211=1, twoway_comparison=1)
  setwd("../")
  save(model.list,efficacy.twoway.211,efficacy.assign.by.quintile.list,
     file=filename)

} else if (regimens == "all"){
  efficacy.combine.by.quintile=process.efficacy(efficacy.assign.by.quintile.list[[1]], efficacy.diff1 = 0.1, event=1, r211=1) #in Crossoversimulationcode
  #efficacy.combine.by.quintile.211=process.efficacy(efficacy.assign.by.quintile.list[[1]], efficacy.diff = 0.1, event=0, r211=1) #in Crossoversimulationcode
  #efficacy.combine.by.quintile.event=process.efficacy(efficacy.assign.by.quintile.list[[1]], efficacy.diff = 0.1, event=1, r211=0) #in Crossoversimulationcode
  #efficacy.combine.randomize=process.efficacy(efficacy.randomize.list[[1]], efficacy.diff = 0.1)
  #efficacy.combine.separate=process.efficacy(efficacy.assign.by.quintile.separate.list[[1]], efficacy.diff = 0.1)
  #efficacy.combine.twostage=process.efficacy(efficacy.assign.by.quintile.twostage.list[[1]], efficacy.diff = 0.1)
  #efficacy.twoway.event=process.efficacy(efficacy.assign.by.quintile.twostage.list[[1]], efficacy.diff1 = 0.1, event=1, r211=0, twoway_comparison=1) #in Crossoversimulationcode
  #efficacy.twoway.211=process.efficacy(efficacy.assign.by.quintile.twostage.list[[1]], efficacy.diff1 = 0.1, event=0, r211=1, twoway_comparison=1) #in Crossoversimulationcode
  #efficacy.twoway=process.efficacy(efficacy.assign.by.quintile.twostage.list[[1]], efficacy.diff = 0.1, event = 1, r211=1, twoway_comparison=0) #in Crossoversimulationcode
  setwd("../")
  save(model.list,efficacy.combine.by.quintile,efficacy.assign.by.quintile.list,#efficacy.assign.by.quintile.twostage.list,
     file=filename)
}
