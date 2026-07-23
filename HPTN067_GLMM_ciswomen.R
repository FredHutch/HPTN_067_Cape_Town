logit = function(x){log(x/(1-x))}
expit=function(x){exp(x)/(1+exp(x))}

assignclass=function(i1,i2,i3){ 
  N.bins=5
  
  class.daily=numeric(length(i1))
  rank.daily=rank(ad.daily, ties.method='random')/length(ad.daily) # rank scaled 0-1
  
  class.time=numeric(length(i2))
  rank.time=rank(ad.time, ties.method='random')/length(ad.time)
  
  class.event=numeric(length(i3))
  rank.event=rank(ad.event, ties.method='random')/length(ad.time) # rank scaled 0-1
  
  class.211=numeric(length(i3)) #using ad.event bc no ad.211
  rank.211=rank(ad.event, ties.method='random')/length(ad.time) # rank scaled 0-1
  
  for(i in seq(N.bins)){
    class.daily[rank.daily<=(i/N.bins)&rank.daily>((i-1)/N.bins)]=i
    class.time[rank.time<=(i/N.bins)&rank.time>((i-1)/N.bins)]=i
    class.event[rank.event<=(i/N.bins)&rank.event>((i-1)/N.bins)]=i
    class.211[rank.211<=(i/N.bins)&rank.211>((i-1)/N.bins)]=i
  }
  
  class.adherence=c(class.daily,class.event,class.time,class.211)
  
} 
#assigns ranks 1-5 based on ad.daily/event/time
#output: class.adherence (class.daily,class.event,class.time) in construct.allsex.capetown (below)
#input: i1=which(ARM.capetown==1),i2=which(ARM.capetown==2),i3=which(ARM.capetown==3)

construct.allsex.capetown=function(i1=which(ARM.capetown==1),i2=which(ARM.capetown==2),i3=which(ARM.capetown==3)){
  
  class.adherence=assignclass(i1,i2,i3) #assigns ranks 1-5 based on ad.daily/event/time
  allsex.capetown.daily=Reduce('rbind',dailysex.capetown[i1]) #combine by rows, each successive individual's observations
  allsex.capetown.time=Reduce('rbind',dailysex.capetown[i2])  #each person's recorded days are the observations, each person sequentially
  allsex.capetown.event=Reduce('rbind',dailysex.capetown[i3])
  
  asad.c.daily=Reduce('rbind',sexactdata.capetown[i1]) #combine by rows, each successive component
  asad.c.time=Reduce('rbind',sexactdata.capetown[i2])  #each observation is a person's sex act, each person sequentially
  asad.c.event=Reduce('rbind',sexactdata.capetown[i3])
  
  allsex.capetown.daily$ARM=1 #adds ARM variable
  allsex.capetown.time$ARM=2
  allsex.capetown.event$ARM=3
  
  asad.c.daily$ARM=1 #adds ARM variable
  asad.c.time$ARM=2
  asad.c.event$ARM=3
  
  as.c=rbind(allsex.capetown.daily,allsex.capetown.time,allsex.capetown.event) #combine by rows: all recorded day observations
  asad.c=rbind(asad.c.daily,asad.c.time,asad.c.event) #combine by rows: all sex acts
  
  as.c$notdaily=(as.c$ARM!=1)
  as.c$eventnosex=as.numeric((as.c$ARM==3)&(as.c$eventpredict==0))
  as.c$eventsex=as.numeric((as.c$ARM==3)&(as.c$eventpredict==1))
  as.c$timenosex=as.numeric((as.c$ARM==2)&(as.c$eventpredict==0))
  as.c$timesex=as.numeric((as.c$ARM==2)&(as.c$eventpredict==1))
  as.c$nondailynosex=as.numeric((as.c$ARM!=1)&(as.c$eventpredict==0))
  as.c$nondailysex=as.numeric((as.c$ARM!=1)&(as.c$eventpredict==1))
  
  as.c$eventsexdayof=as.numeric((as.c$ARM==3)&(as.c$sexdayof>0))
  as.c$timesexdayof=as.numeric((as.c$ARM==2)&(as.c$sexdayof>0))
  
  as.c$eventnosexdayof=as.numeric((as.c$ARM==3)&(as.c$sexdayof==0))
  as.c$timenosexdayof=as.numeric((as.c$ARM==2)&(as.c$sexdayof==0))
  
  as.c$eventsexnextday=as.numeric((as.c$ARM==3)&(as.c$sexnextday>0))
  as.c$timesexnextday=as.numeric((as.c$ARM==2)&(as.c$sexnextday>0))
  as.c$nondailysexnextday=as.numeric((as.c$ARM!=1)&(as.c$sexnextday>0))
  
  as.c$eventnosexnextday=as.numeric((as.c$ARM==3)&(as.c$sexnextday==0))
  as.c$timenosexnextday=as.numeric((as.c$ARM==2)&(as.c$sexnextday==0))
  
  as.c$eventnoanalsex=as.numeric((as.c$ARM==3)&(as.c$eventpredictanal==0))
  as.c$eventanalsex=as.numeric((as.c$ARM==3)&(as.c$eventpredictanal==1))
  as.c$timenoanalsex=as.numeric((as.c$ARM==2)&(as.c$eventpredictanal==0))
  as.c$timeanalsex=as.numeric((as.c$ARM==2)&(as.c$eventpredictanal==1))
  as.c$nondailynoanalsex=as.numeric((as.c$ARM!=1)&(as.c$eventpredictanal==0))
  as.c$nondailyanalsex=as.numeric((as.c$ARM!=1)&(as.c$eventpredictanal==1))
  
  as.c$eventanalsexdayof=as.numeric((as.c$ARM==3)&(as.c$analsexdayof>0))
  as.c$timeanalsexdayof=as.numeric((as.c$ARM==2)&(as.c$analsexdayof>0))
  
  as.c$eventnoanalsexdayof=as.numeric((as.c$ARM==3)&(as.c$analsexdayof==0))
  as.c$timenoanalsexdayof=as.numeric((as.c$ARM==2)&(as.c$analsexdayof==0))
  
  as.c$eventanalsexnextday=as.numeric((as.c$ARM==3)&(as.c$analsexnextday>0))
  as.c$timeanalsexnextday=as.numeric((as.c$ARM==2)&(as.c$analsexnextday>0))
  as.c$nondailyanalsexnextday=as.numeric((as.c$ARM!=1)&(as.c$analsexnextday>0))
  
  as.c$eventnoanalsexnextday=as.numeric((as.c$ARM==3)&(as.c$analsexnextday==0))
  as.c$timenoanalsexnextday=as.numeric((as.c$ARM==2)&(as.c$analsexnextday==0))
  
  asad.c.pre=asad.c
  asad.c.post=asad.c
  
  asad.c.pre$combopill=asad.c.pre$pretime>=(-1/12) #combo pill if within +-2 hours of sex
  asad.c.post$combopill=asad.c.pre$posttime<=(1/12)
  asad.c.pre$ispost=FALSE
  asad.c.post$ispost=TRUE
  
  asad.c.prepost=rbind(asad.c.pre,asad.c.post) #acts vector twice, pre first then post
  
  asad.c.prepost$prescribed=(asad.c.prepost$ARM!=1&asad.c.prepost$ispost)|asad.c.prepost$ARM==3
  
  asad.c.prepost$takepill=asad.c.prepost$combopill
  
  as.c1=as.c #obs= all recorded days for an individual
  as.c2=merge(as.c,asad.c.prepost,by=c('day','id','ARM'),all.y=TRUE) #obs=sex days twice pre&post
  
  as.c1[,setdiff(names(asad.c.prepost),names(as.c))]=NA #variables only in asad.c.prepost now all NA
  as.c1$type=ifelse(as.c1$eventpredict,'predict',ifelse(as.c1$nopostdaybefore,'makeup','notpredict'))
  as.c1$type2=as.c1$type
  as.c2$type=ifelse(as.c2$ispost,'post',ifelse(as.c2$prepill,'extrapre','neededpre')) #? confused about needed/extra pre
  as.c2$type2="combo"
  
  as.c1$timetype=ifelse(as.c1$pilllasttwo==0,"needed",ifelse(as.c1$nopostdaybefore,"makeup","notneeded"))
  as.c2$timetype=ifelse(as.c2$ispost,"post","pre")
  
  as.c1$takepill=as.c1$pillonly>0
  as.c2$takepill=as.c2$combopill
  
  ap=rbind(as.c1,as.c2)
  list(allsex.capetown=as.c,allsexactdata.capetown=asad.c, 
       allsexactdata.capetown.prepost=asad.c.prepost,
       allpill=ap, 
       class.adherence=class.adherence)
} 
#creates as.c, asad.c, as.c1, as.c2 vectors of acts and recorded days from 067 data
#output: allsex.capetown (all recorded days), allsexactdata.capetown (all sex acts),
#   allsexactdata.capetown.prepost, allpill, class.adherence in fit.models.all (below)
#calls: assignclass (above)
#input: i1=which(ARM.capetown==1),i2=which(ARM.capetown==2),i3=which(ARM.capetown==3)

fit.pill.daily=function(allsex){
  glm.pill.daily=glmer(I(pillonly+combopre+combopost+comboprepost)>0~0+(0+pilldaybefore|id)+pilldaybefore,allsex,
                       subset=ARM==1,
                       family='binomial')
  ran.daily=as.data.frame(coefficients(glm.pill.daily)$id) #pill day before T/F
  ran.daily=ran.daily
  ran.daily$id=as.numeric(row.names(ran.daily))
  list(glm.pill.daily=glm.pill.daily,ran.daily=ran.daily)
} 
#fits generalized linear mixed-effects model
#data=data$allsex.capetown
#response:(pillonly+combopre+combopost+comboprepost)>0
#terms:0+(0+pilldaybefore|id)+pilldaybefore conditional on individual

fit.pill.time=function(allpill){
  glm.pill.time=glmer(takepill~0 + (0+timetype|id)+timetype,allpill,
                      subset=ARM==2,
                      family='binomial',
                      control=glmerControl(optimizer='bobyqa'))
  ran.time=as.data.frame(coefficients(glm.pill.time)$id)
  ran.time=ran.time
  ran.time$id=as.numeric(row.names(ran.time))
  list(glm.pill.time=glm.pill.time,ran.time=ran.time)
}
#fits generalized linear mixed-effects model
#data=data$allpill
#response:takepill
#terms:0 + (0+timetype|id)+timetype conditional on individual

fit.pill.event=function(allpill){
  glm.pill.event=glmer(takepill~0+(0+type|id)+type,allpill,
                       subset=ARM==3,
                       family='binomial',
                       control=glmerControl(optimizer='bobyqa'))
  ran.event=as.data.frame(coefficients(glm.pill.event)$id)
  
  ran.event=ran.event
  
  p0=expit(ran.event$typepredict) #predict = sex day of or day after
  prepill=1-(1-p0)^2
  combopre=expit(ran.event$typeneededpre)
  precoverage=prepill+(1-prepill)*combopre
  
  ran.event$sensitivity=prepill/precoverage
  ran.event$pre.adherence.logit=logit(precoverage)
  p3=expit(ran.event$typenotpredict) #not predict = no sex day before, day of, day after
  ran.event$specificity= 1 - p3
  ran.event$id=as.numeric(row.names(ran.event))
  
  list(glm.pill.event=glm.pill.event,ran.event=ran.event)
}
#fits generalized linear mixed-effects model
#data=data$allpill
#response:takepill
#terms:0+(0+type|id)+type conditional on individual
#ran.event with expit transforms?

#allpill=data$allpill
fit.pill.211=function(allpill){
  glm.pill.211=glmer(takepill~0+(0+type|id)+type,allpill, #since using type, output = type(each variable under type)
                       subset=ARM==3,
                       family='binomial',
                       control=glmerControl(optimizer='bobyqa'))
  ran.211=as.data.frame(coefficients(glm.pill.211)$id)
  
  ran.211=ran.211
  
  p0=expit(ran.211$typepredict) #predict = sex day of or day after
  prepill=1-(1-p0)^2
  combopre=expit(ran.211$typeneededpre)
  precoverage=prepill+(1-prepill)*combopre
  
  ran.211$sensitivity=prepill/precoverage
  ran.211$pre.adherence.logit=logit(precoverage)
  p3=expit(ran.211$typenotpredict) #not predict = no sex day before, day of, day after
  ran.211$specificity= 1 - p3
  ran.211$id=as.numeric(row.names(ran.211))
  
  list(glm.pill.211=glm.pill.211,ran.211=ran.211)
}
#fits generalized linear mixed-effects model
#data=data$allpill
#response:takepill
#terms:0+(0+type|id)+type conditional on individual
#ran.event with expit transforms?

fit.sex.all=function(data){
  with(data,{
    idlist=as.character(seq(179))
    glm.sexfrequency=glmer(I(combopre+combopost+comboprepost+sexonly)~(1|id),allsex.capetown, family='poisson')
    
    sexfrequency=coefficients(glm.sexfrequency)$id[idlist,1]
    
    glm.analfrequency=glmer(anal~(1|id),allsexactdata.capetown, family='binomial')
    
    
    glm.condom=glmer(condom~0+(0+anal|id)+anal,allsexactdata.capetown, family='binomial')
    
    condomanal=coefficients(glm.condom)$id[idlist,'analTRUE']
    condomvaginal=coefficients(glm.condom)$id[idlist,'analFALSE']
    
    
    ranef.af=ranef(glm.analfrequency)$id
    
    
    hifreq=as.numeric(rownames(ranef.af)[which(ranef.af[,1]>-2)])
    
    #glm.analfrequency=glmer(anal~(1|id),allsexactdata.capetown, family='binomial', subset=!id%in%hifreq)#replaces above?
    #in capetown data, hifreq contains all ids
    analfrequency=expit(coefficients(glm.analfrequency)$id[idlist,1])
    
    analfrequency[hifreq]=1 #gives anal freq a max of 1
    
    
    sexparams.ind=data.frame(id=idlist,
                             analfrequency=analfrequency,
                             sexfrequency=sexfrequency,
                             condomanal=condomanal,
                             condomvaginal=condomvaginal)
    
    
    list(mer.mod.sex=glm.sexfrequency,
         mer.mod.analfrequency=glm.analfrequency,
         analonly=mean(ranef.af[,1]>-2),
         mer.mod.condom=glm.condom,
         ran.sex=sexparams.ind)
  })
}
#fits generalized linear mixed-effects model
#for sexfrequency, analfrequency, anal
#data=data$allsex.capetown
#response:combopre+combopost+comboprepost+sexonly
#terms:(1|id) conditional on individual

fit.models.all=function(data){
  with(as.list(data),{
    list(glm.pill.daily=fit.pill.daily(allsex.capetown),
         glm.pill.time=fit.pill.time(allpill),
         glm.pill.event=fit.pill.event(allpill),
         glm.pill.211=fit.pill.211(allpill),
         glm.sex=fit.sex.all(data))
  })
} 
#output: model.list in get.random.params, get.HPTN067.params (below)
#calls: fit.pill.daily .time .event and fit.sex.all
#input: data (from construct.allsex.capetown() above) 

random.matrix=function(N,mer.mod){
  SIG=VarCorr(mer.mod)[[1]] #calculates estimated variance, sd, correlations between random effect terms
  MU=fixef(mer.mod) #extracts fixed-effect estimates
  m=nrow(SIG)
  X=matrix(rnorm(N*m),nrow=N,ncol=m) #creates matrix with normal dist, mean=0, sd=1
  Y=X%*%sqrtm(SIG)
  colnames(Y)=row.names(SIG)
  t(t(Y)+MU)
} 
#combines random and fixed effect estimates into matrix of specified size
#output: ran.daily/time/event/... (in get.random.params below)
#input: model.list$glm.pill.daily/time/event, model.list$glm.sex (from fit.models.all above)

get.random.params=function(N,model.list){
  with(model.list,{
    ran.daily=random.matrix(N,glm.pill.daily$glm.pill.daily)
    ran.time=random.matrix(N,glm.pill.time$glm.pill.time)
    ran.event=random.matrix(N,glm.pill.event$glm.pill.event)
    ran.211=random.matrix(N,glm.pill.211$glm.pill.211)
    ran.sex.freq=random.matrix(N,glm.sex$mer.mod.sex)
    ran.anal=random.matrix(N,glm.sex$mer.mod.analfrequency)
    ran.condom=random.matrix(N,glm.sex$mer.mod.condom)
    ran.analonly=rbinom(N,size=1,prob=mean(glm.sex$analonly))
    ran.anal[ran.analonly==1]=1 #change ran.anal to 1 if binomial dist=1
    
    ran.sex=Reduce('cbind',list(ran.sex.freq,ran.anal,ran.condom)) #combine sex,anal,condom
    
    colnames(ran.sex)=c('sexfrequency','analfrequency','condomvaginal','condomanal')
    
    list(params.daily=as.data.frame(cbind(ran.daily,ran.sex)),
         params.time=as.data.frame(cbind(ran.time,ran.sex)),
         params.event=as.data.frame(cbind(ran.event,ran.sex)),
         params.211=as.data.frame(cbind(ran.211,ran.sex)))
  })
} 
#creates matrices for pills, sex freq, anal freq, condom vaginal, condom anal 
#output: params in get.efficacy in CrossoverSimulationcode.R and in reassign.efficacy.by.class in run sim 
#calls: random.matrix 
#input: data (from glms): model.list$glm.pill.daily/time/event, model.list$glm.sex (from fit.models.all above)

get.HPTN067.params=function(model.list){
  with(model.list,{
    list(params.daily=merge(glm.pill.daily$ran.daily,glm.sex$ran.sex),
         params.time=merge(glm.pill.time$ran.time,glm.sex$ran.sex),
         params.event=merge(glm.pill.event$ran.event,glm.sex$ran.sex),
         params.211=merge(glm.pill.211$ran.211,glm.sex$ran.sex))
  })
}
#combines glm.pill.daily/time/event$ran.daily/time/event and glm.sex$ran.sex
#output: params.067 in get.efficacy (in CrossoverSimulationcode.R) and in reassign.efficacy.by.class (in run sim) 
#input: data (from glms): model.list$glm.pill.daily/time/event, model.list$glm.sex (from fit.models.all above)

random.params=function(){
  risk.anal=rnorm(1,DATAanalmed,sd=DATAanalse)
  risk.vaginal=rnorm(1,DATAvaginalmed,sd=DATAvaginalse)
  condom.efficacy=rnorm(1,CONDOMeffmed,sd=CONDOMeffse)
  list(risk.anal=expit(risk.anal),
       risk.vaginal=expit(risk.vaginal),
       condom.efficacy=expit(condom.efficacy))
}#not called
#create normal distributions with means DATAanalmed, DATAvaginalmed, CONDOMeffmed