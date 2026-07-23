get.efficacy=function(params.list, N.days=180){
  with(params.list,{
    list(efficacy.daily=adply(params.daily,1,calc.efficacy.daily.rand, N.days=N.days),
         efficacy.time=adply(params.time,1,calc.efficacy.time, N.days=N.days),
         efficacy.event=adply(params.event,1,calc.efficacy.event, N.days=N.days),
         efficacy.211=adply(params.211,1,calc.efficacy.211, N.days=N.days))
  })
}
#calls efficacy calculators
#output: efficacy.067 and efficacy.rand in reassign.efficacy.by.class in run sim
#calls: calc.efficacy.daily.rand, calc.efficacy.time, calc.efficacy.event (in EfficacyCalculation4.R)
#input: params.067/params.rand: params.daily, params.time. params.event (from get.HPTN067/random.params in HPTN067GLMM.R)

assign.adherence.class=function(N.bins, N.classes, efficacy.df, attribute){
  rank.adherence=rank(efficacy.df[ ,attribute], ties.method='random')
  class=ceiling(N.classes*rank.adherence/length(rank.adherence))
  bin=ceiling(N.bins*rank.adherence/length(rank.adherence))
  efficacy.df$class.adherence=class
  efficacy.df$bin.adherence=bin
  efficacy.df
}
#creates efficacy.df with bin & class based on rank.adherence
#output: efficacy.daily/event/time in reassign.efficacy.byclass (below)
#input: efficacy.df=efficacy.daily/event/time, attribute=attribute.daily/event/time

reorder.subjects=function(class.original,class.new){
  bins=unique(class.original)
  
  id=class.original
  for(i in bins){
    class.i=which(class.original==i)
    new.class.i=which(class.new==i)
    id[class.i]=sample(new.class.i, length(class.i), replace=TRUE)
  }
  id
}
#new id from a sample n=length(class.i) 
#output: id.time/event in reassign.efficacy.byclass (below)
#input: efficacy.daily$bin.adherence,efficacy.time/event$bin.adherence

reassign.efficacy.byclass=function(efficacy.list,params.list, N.bins = 5, N.classes = 5, N.days=180,
                                   attribute.daily = "adherence", attribute.time = "adherence", attribute.event = "adherence",
                                   linkevent.to.time = FALSE, attribute.time2 = "adherence.sex",attribute.211 = "adherence"){
  with(as.list(c(params.list,efficacy.list)),{
    efficacy.daily=assign.adherence.class(N.bins, N.classes, efficacy.daily, attribute.daily)
    efficacy.event=assign.adherence.class(N.bins, N.classes, efficacy.event, attribute.event)
    efficacy.time=assign.adherence.class(N.bins, N.classes, efficacy.time, attribute.time)
    efficacy.211=assign.adherence.class(N.bins, N.classes, efficacy.211, attribute.211)
    
    id.time=reorder.subjects(efficacy.daily$bin.adherence,efficacy.time$bin.adherence)
    
    id.event=reorder.subjects(efficacy.daily$bin.adherence,efficacy.event$bin.adherence)
    
    id.211=reorder.subjects(efficacy.daily$bin.adherence,efficacy.211$bin.adherence)
    
    if(linkevent.to.time){
      efficacy.time2=assign.adherence.class(N.bins, N.classes, efficacy.time, attribute.time2)
      id.event=reorder.subjects(efficacy.time2$bin.adherence[id.time], efficacy.event$bin.adherence)
      id.211=reorder.subjects(efficacy.time2$bin.adherence[id.time], efficacy.211$bin.adherence)
    }
    
    n.time=setdiff(colnames(params.time),colnames(params.daily))
    n.event=setdiff(colnames(params.event),colnames(params.daily))
    n.211=setdiff(colnames(params.211),colnames(params.daily))
    params.combine=cbind(params.daily, params.time[id.time,n.time], params.event[id.event,n.event],params.211[id.211,n.211])

    efficacy.time=adply(params.combine,1,calc.efficacy.time, N.days=N.days)
    efficacy.event=adply(params.combine,1,calc.efficacy.event, N.days=N.days)
    efficacy.211=adply(params.combine,1,calc.efficacy.211, N.days=N.days)
    
    efficacy.daily.followup=adply(params.combine,1,calc.efficacy.daily.rand, N.days=N.days)
    efficacy.time.followup=adply(params.combine,1,calc.efficacy.time, N.days=N.days)
    efficacy.event.followup=adply(params.combine,1,calc.efficacy.event, N.days=N.days)
    efficacy.211.followup=adply(params.combine,1,calc.efficacy.211, N.days=N.days)
    
    names(efficacy.daily)=sapply(names(efficacy.daily),function(x){paste(x,'daily',sep='.')})
    names(efficacy.time)=sapply(names(efficacy.time),function(x){paste(x,'time',sep='.')})
    names(efficacy.event)=sapply(names(efficacy.event),function(x){paste(x,'event',sep='.')})
    names(efficacy.211)=sapply(names(efficacy.211),function(x){paste(x,'211',sep='.')})
    names(efficacy.daily.followup)=sapply(names(efficacy.daily.followup),function(x){paste(x,'daily.followup',sep='.')})
    names(efficacy.time.followup)=sapply(names(efficacy.time.followup),function(x){paste(x,'time.followup',sep='.')})
    names(efficacy.event.followup)=sapply(names(efficacy.event.followup),function(x){paste(x,'event.followup',sep='.')})
    names(efficacy.211.followup)=sapply(names(efficacy.211.followup),function(x){paste(x,'211.followup',sep='.')})
    
    list(cbind(efficacy.daily,efficacy.time,efficacy.event,efficacy.211,efficacy.daily.followup,efficacy.time.followup,efficacy.event.followup,efficacy.211.followup),params.combine)
  })
}
#output: crossover simulations in run sim
#calls: assign.adherence.class, reorder.subjects (above), calc.efficacy.time/event/daily.rand (in EfficacyCalculation4.R)
#input: efficacy (from get.efficacy above), params (from get.HPTN067/random.params in HPTN067GLMM.R)

combine.efficacy=function(efficacy.list){
  with(efficacy.list,{
    efficacy.daily$ARM=1
    efficacy.time$ARM=2
    efficacy.event$ARM=3
    efficacy.211$ARM=3 #4?
    neededcols=c('analfrequency','sexfrequency','pillrate','efficacy','condomanal','condomvaginal','id','ARM','sensitivity','specificity')
    efficacy.all=rbind(efficacy.daily[,neededcols],efficacy.time[,neededcols],efficacy.event[,neededcols],efficacy.211[,neededcols])
    efficacy.all$class=class.adherence[efficacy.all$id]
    efficacy.all
  })
}
#combines efficacy.daily/time/event into efficacy.all with columns 
# 'analfrequency','sexfrequency','pillrate','efficacy','condomanal','condomvaginal','id','ARM','sensitivity','specificity'
#output: efficacy.all
#input: efficacy.list

compare.regimes=function(ef1,ef2,p1,p2,prf1,prf2,ed1,ed2){
  cond1=ef1<(ef2+ed1)&(p2<p1*prf1) #regime 2= efficacy higher or within 0.1 & pills less than half regime 1's: [efficacy1 < (efficacy2+0.1)] & [pillrate2 < (pillrate1*0.5)]
  cond2=ef1<ef2&p2<p1              #regime 2= higher efficacy & fewer pills: (efficacy1 < efficacy2) & (pillrate2 < pillrate1)
  cond3=ef1<(ef2-ed2)&(p2<p1/prf2) #regime 2= efficacy 0.1+ higher & pills less than double regime 1's: [efficacy1 < (efficacy2-0.1)] & [pillrate2 < (pillrate1/0.5)]
  
  cond1|cond2|cond3 # regime 2: efficacy close/better & much fewer pills|higher efficacy & fewer pills|much higher efficacy & pills ok
}  
#compares regimes based on 3 conditions
#output: timeoverdaily, eventoverdaily, timeoverevent in process.efficacy (below)
#input: ef1=efficacy.df$efficacy.daily/daily/event,ef2=efficacy.df$efficacy.time/event/time,
#   p1=efficacy.df$pillrate.daily/daily/event,p2=efficacy.df$pillrate.time/event/time,prf=pillratefactor,ed=efficacy.diff

compare.regimes2=function(ef1,ef2,p1,p2,prf1,prf2,ed1,ed2){
  cond1=ef1<(ef2+ed1)&(p2<p1*prf1) #regime 2= efficacy higher or within 0.1 & pills less than half regime 1's: [efficacy1 < (efficacy2+0.1)] & [pillrate2 < (pillrate1*0.5)]
  cond2=ef1<ef2&p2<p1            #regime 2= higher efficacy & fewer pills: (efficacy1 < efficacy2) & (pillrate2 < pillrate1)
  cond3=ef1<(ef2-ed2)&(p2<p1/prf2) #regime 2= efficacy 0.1+ higher & pills less than double regime 1's: [efficacy1 < (efficacy2-0.1)] & [pillrate2 < (pillrate1/0.5)]
  #reason 1= regime 2: efficacy close/better & much fewer pills
  #reason 2= higher efficacy & fewer pills
  #reason 3= much higher efficacy & pills ok
  reason = ifelse(cond2,2,ifelse(cond1,1,ifelse(cond3,3,0)))
  reason
}  
#compares regimes based on 3 conditions, outputs reason that one is better
#output: timeoverdaily, eventoverdaily, timeoverevent in process.efficacy (below)
#input: ef1=efficacy.df$efficacy.daily/daily/event,ef2=efficacy.df$efficacy.time/event/time,
#   p1=efficacy.df$pillrate.daily/daily/event,p2=efficacy.df$pillrate.time/event/time,prf=pillratefactor,ed=efficacy.diff

compare.regimes.adherence=function(ad.rate,ad.target){
  cond1 = ad.rate < ad.target
  cond1 
}  

process.efficacy=function(efficacy.df, pillratefactor1=0.5, efficacy.diff1=0.1, pillratefactor2=pillratefactor.eff, efficacy.diff2=efficacy.diff.eff, event=1, r211=1, twoway_comparison=0,ad.target=adherence.target){
  if (adherence.target > 0) {
      nondailyoverdaily=compare.regimes.adherence(efficacy.df$adherence.daily,ad.target)
    
      if (r211==0 & twoway_comparison==1){ 
        efficacy.df$bestarm=ifelse(nondailyoverdaily,3,1)
      } else {
        efficacy.df$bestarm=ifelse(nondailyoverdaily,4,1)
      }
    
    } else {
      timeoverdaily=compare.regimes(efficacy.df$efficacy.daily,efficacy.df$efficacy.time,efficacy.df$pillrate.daily,efficacy.df$pillrate.time,pillratefactor1,pillratefactor2,efficacy.diff1,efficacy.diff2)
      eventoverdaily=compare.regimes(efficacy.df$efficacy.daily,efficacy.df$efficacy.event,efficacy.df$pillrate.daily,efficacy.df$pillrate.event,pillratefactor1,pillratefactor2,efficacy.diff1,efficacy.diff2)
      r211overdaily=compare.regimes(efficacy.df$efficacy.daily,efficacy.df$efficacy.211,efficacy.df$pillrate.daily,efficacy.df$pillrate.211,pillratefactor1,pillratefactor2,efficacy.diff1,efficacy.diff2)
      timeoverevent=compare.regimes(efficacy.df$efficacy.event,efficacy.df$efficacy.time,efficacy.df$pillrate.event,efficacy.df$pillrate.time,pillratefactor1,pillratefactor2,efficacy.diff1,efficacy.diff2)
      timeover211=compare.regimes(efficacy.df$efficacy.211,efficacy.df$efficacy.time,efficacy.df$pillrate.211,efficacy.df$pillrate.time,pillratefactor1,pillratefactor2,efficacy.diff1,efficacy.diff2)
      eventover211=compare.regimes(efficacy.df$efficacy.211,efficacy.df$efficacy.event,efficacy.df$pillrate.211,efficacy.df$pillrate.event,pillratefactor1,pillratefactor2,efficacy.diff1,efficacy.diff2)
  
      if (event == 1 & r211==1 & twoway_comparison==0){ 
        efficacy.df$bestarm=ifelse(eventoverdaily,ifelse(timeoverevent,ifelse(timeover211,2,4),ifelse(eventover211,3,4)),ifelse(timeoverdaily,ifelse(timeover211,2,4),ifelse(r211overdaily,4,1)))
      } else if (r211==0 & twoway_comparison==0){ 
        efficacy.df$bestarm=ifelse(eventoverdaily,ifelse(timeoverevent,2,3),ifelse(timeoverdaily,2,1))
      } else if (event==0 & twoway_comparison==0){
        efficacy.df$bestarm=ifelse(r211overdaily,ifelse(timeover211,2,4),ifelse(timeoverdaily,2,1))
      } else if (r211==0 & twoway_comparison==1){ 
      efficacy.df$bestarm=ifelse(eventoverdaily,3,1)
      } else {
      efficacy.df$bestarm=ifelse(r211overdaily,4,1)
      }
  }
  newefficacy=efficacy.df$efficacy.daily #assigns new efficacy to daily
  newefficacy.followup=efficacy.df$efficacy.daily.followup
  pillrate = efficacy.df$pillrate.daily
  pillrate.followup = efficacy.df$pillrate.daily.followup
  reason = numeric(length(pillrate))
  
  newefficacy[efficacy.df$bestarm==2] = efficacy.df$efficacy.time[efficacy.df$bestarm==2] #reassigns new efficacy if another regime is better
  newefficacy[efficacy.df$bestarm==3] = efficacy.df$efficacy.event[efficacy.df$bestarm==3]
  newefficacy[efficacy.df$bestarm==4] = efficacy.df$efficacy.211[efficacy.df$bestarm==4]
  
  newefficacy.followup[efficacy.df$bestarm==2] = efficacy.df$efficacy.time.followup[efficacy.df$bestarm==2] #reassigns if another regime is better
  newefficacy.followup[efficacy.df$bestarm==3] = efficacy.df$efficacy.event.followup[efficacy.df$bestarm==3]
  newefficacy.followup[efficacy.df$bestarm==4] = efficacy.df$efficacy.211.followup[efficacy.df$bestarm==4]
  
  pillrate[efficacy.df$bestarm==2] = efficacy.df$pillrate.time[efficacy.df$bestarm==2] #reassigns if another regime is better
  pillrate[efficacy.df$bestarm==3] = efficacy.df$pillrate.event[efficacy.df$bestarm==3]
  pillrate[efficacy.df$bestarm==4] = efficacy.df$pillrate.211[efficacy.df$bestarm==4]
  
  pillrate.followup[efficacy.df$bestarm==2] = efficacy.df$pillrate.time.followup[efficacy.df$bestarm==2] #reassigns if another regime is better
  pillrate.followup[efficacy.df$bestarm==3] = efficacy.df$pillrate.event.followup[efficacy.df$bestarm==3]
  pillrate.followup[efficacy.df$bestarm==4] = efficacy.df$pillrate.211.followup[efficacy.df$bestarm==4]
  
  #reason other regimes beat daily
  timeoverdailyreason =compare.regimes2(efficacy.df$efficacy.daily,efficacy.df$efficacy.time,efficacy.df$pillrate.daily,efficacy.df$pillrate.time,pillratefactor1,pillratefactor2,efficacy.diff1,efficacy.diff2)
  eventoverdailyreason=compare.regimes2(efficacy.df$efficacy.daily,efficacy.df$efficacy.event,efficacy.df$pillrate.daily,efficacy.df$pillrate.event,pillratefactor1,pillratefactor2,efficacy.diff1,efficacy.diff2)
  r211overdailyreason =compare.regimes2(efficacy.df$efficacy.daily,efficacy.df$efficacy.211,efficacy.df$pillrate.daily,efficacy.df$pillrate.211,pillratefactor1,pillratefactor2,efficacy.diff1,efficacy.diff2)
  timeovereventreason =compare.regimes2(efficacy.df$efficacy.event,efficacy.df$efficacy.time,efficacy.df$pillrate.event,efficacy.df$pillrate.time,pillratefactor1,pillratefactor2,efficacy.diff1,efficacy.diff2)
  timeover211reason   =compare.regimes2(efficacy.df$efficacy.211,efficacy.df$efficacy.time,efficacy.df$pillrate.211,efficacy.df$pillrate.time,pillratefactor1,pillratefactor2,efficacy.diff1,efficacy.diff2)
  eventover211reason  =compare.regimes2(efficacy.df$efficacy.211,efficacy.df$efficacy.event,efficacy.df$pillrate.211,efficacy.df$pillrate.event,pillratefactor1,pillratefactor2,efficacy.diff1,efficacy.diff2)
  eventovertimereason =compare.regimes2(efficacy.df$efficacy.time,efficacy.df$efficacy.event,efficacy.df$pillrate.time,efficacy.df$pillrate.event,pillratefactor1,pillratefactor2,efficacy.diff1,efficacy.diff2)
  r211overeventreason =compare.regimes2(efficacy.df$efficacy.event,efficacy.df$efficacy.211,efficacy.df$pillrate.event,efficacy.df$pillrate.211,pillratefactor1,pillratefactor2,efficacy.diff1,efficacy.diff2)
  r211overtimereason  =compare.regimes2(efficacy.df$efficacy.time,efficacy.df$efficacy.211,efficacy.df$pillrate.time,efficacy.df$pillrate.211,pillratefactor1,pillratefactor2,efficacy.diff1,efficacy.diff2)
  
  reason[efficacy.df$bestarm==2] = timeoverdailyreason[efficacy.df$bestarm==2]+((timeoverdailyreason[efficacy.df$bestarm==2]==0)*timeovereventreason[efficacy.df$bestarm==2]+(timeoverdailyreason[efficacy.df$bestarm==2]==0)*(timeovereventreason[efficacy.df$bestarm==2]==0)*timeover211reason[efficacy.df$bestarm==2])*(twoway_comparison==0)  
  reason[efficacy.df$bestarm==3] = eventoverdailyreason[efficacy.df$bestarm==3]+((eventoverdailyreason[efficacy.df$bestarm==3]==0)*eventovertimereason[efficacy.df$bestarm==3]+(eventoverdailyreason[efficacy.df$bestarm==3]==0)*(eventovertimereason[efficacy.df$bestarm==3]==0)*eventover211reason[efficacy.df$bestarm==3])*(twoway_comparison==0) 
  reason[efficacy.df$bestarm==4] = r211overdailyreason[efficacy.df$bestarm==4]+((r211overdailyreason[efficacy.df$bestarm==4]==0)*r211overeventreason[efficacy.df$bestarm==4]+(r211overdailyreason[efficacy.df$bestarm==4]==0)*(r211overeventreason[efficacy.df$bestarm==4]==0)*r211overtimereason[efficacy.df$bestarm==4])*(twoway_comparison==0)  
  
  efficacy.df$newefficacy = newefficacy
  efficacy.df$newefficacy.followup = newefficacy.followup
  efficacy.df$pillrate = pillrate
  efficacy.df$pillrate.followup = pillrate.followup
  efficacy.df$reason = reason
  
  efficacy.df
}
#assigns the best efficacy regime (daily/event/time) 
#output: efficacy.combine.by.quintile/efficacy.combine.randomize/... (in runsim)
#calls: compare.regimes (above)
#input: efficacy.assign.by.quintile.list[[1]]/efficacy.randomize.list[[1]]/..., efficacy.diff = 0.1 (in runsim)

getclass=function(x,Nbins){
  ceiling(Nbins*rank(x, ties.method='random')/length(x))
} 
#returns the bin for x
#output: efficacy.combine.twostage$sex.class/specificty.class/sensitivity.class in runsim
#input: efficacy.combine.twostage$sexacts.daily/efficacy.combine.twostage$specificity.event,efficacy.combine.twostage$sensitivity.event, 3 from process.efficacy (above)

process.subset=function(efficacy.df){
  samples=nrow(efficacy.df)
  choose.daily=sum(efficacy.df$bestarm==1)/samples #% in each regime
  choose.time=sum(efficacy.df$bestarm==2)/samples
  choose.event=sum(efficacy.df$bestarm==3)/samples
  choose.211=sum(efficacy.df$bestarm==4)/samples
  
  #assign efficacy/sexacts/... based on best regime
  efficacy=ifelse(efficacy.df$bestarm==1,efficacy.df$efficacy.daily,ifelse(efficacy.df$bestarm==2,efficacy.df$efficacy.time,ifelse(efficacy.df$bestarm==3,efficacy.df$efficacy.event,efficacy.df$efficacy.211)))
  
  sexacts=ifelse(efficacy.df$bestarm==1,efficacy.df$sexacts.daily,ifelse(efficacy.df$bestarm==2,efficacy.df$sexacts.time,ifelse(efficacy.df$bestarm==3,efficacy.df$sexacts.event,efficacy.df$sexacts.211)))
  
  pillrate=ifelse(efficacy.df$bestarm==1,efficacy.df$pillrate.daily,ifelse(efficacy.df$bestarm==2,efficacy.df$pillrate.time,ifelse(efficacy.df$bestarm==3,efficacy.df$pillrate.event,efficacy.df$pillrate.211)))
  
  efficacy.followup=ifelse(efficacy.df$bestarm==1,efficacy.df$efficacy.daily.followup,ifelse(efficacy.df$bestarm==2,efficacy.df$efficacy.time.followup,ifelse(efficacy.df$bestarm==3,efficacy.df$efficacy.event.followup,efficacy.df$efficacy.211.followup)))
  
  sexacts.followup=ifelse(efficacy.df$bestarm==1,efficacy.df$sexacts.daily.followup,ifelse(efficacy.df$bestarm==2,efficacy.df$sexacts.time.followup,ifelse(efficacy.df$bestarm==3,efficacy.df$sexacts.event.followup,efficacy.df$sexacts.211.followup)))
  
  pillrate.followup=ifelse(efficacy.df$bestarm==1,efficacy.df$pillrate.daily.followup,ifelse(efficacy.df$bestarm==2,efficacy.df$pillrate.time.followup,ifelse(efficacy.df$bestarm==3,efficacy.df$pillrate.event.followup,efficacy.df$pillrate.211.followup)))
  
  #calculate original and new efficacy and pillrate
  efficacy.orig=sum(efficacy.df$efficacy.daily*efficacy.df$sexacts.daily)/sum(efficacy.df$sexacts.daily)
  
  pillrate.orig=mean(efficacy.df$pillrate.daily)
  
  efficacy.new=sum(efficacy*sexacts)/sum(sexacts)
  
  pillrate.new=mean(pillrate)
  
  efficacy.orig.followup=sum(efficacy.df$efficacy.daily.followup*efficacy.df$sexacts.daily.followup)/sum(efficacy.df$sexacts.daily.followup)
  
  pillrate.orig.followup=mean(efficacy.df$pillrate.daily.followup)
  
  efficacy.new.followup=sum(efficacy.followup*sexacts.followup)/sum(sexacts.followup)
  
  pillrate.new.followup=mean(pillrate.followup)
  #calculate efficacy and pillrate of switchers
  switch = efficacy.df$bestarm!=1
  pillrate.orig.switch = mean(efficacy.df$pillrate.daily[switch])
  pillrate.new.switch = mean(pillrate[switch])
  
  efficacy.orig.switch = sum(efficacy.df$efficacy.daily[switch]*efficacy.df$sexacts.daily[switch])/sum(efficacy.df$sexacts.daily[switch])
  efficacy.new.switch = sum(efficacy[switch]*sexacts[switch])/sum(sexacts[switch])
  
  efficacy.orig.personal = mean(efficacy.df$efficacy.daily)
  efficacy.orig.personal.sd = sd(efficacy.df$efficacy.daily)
  
  efficacy.followup.personal = mean(efficacy.followup)
  efficacy.followup.personal.sd = sd(efficacy.followup)
  
  pillrate.orig.sd = sd(efficacy.df$pillrate.daily)
  pillrate.new.followup.sd = sd(pillrate.followup)
  
  efficacy.orig.noswitch = sum(efficacy.df$efficacy.daily[!switch]*efficacy.df$sexacts.daily[!switch])/sum(efficacy.df$sexacts.daily.followup[!switch])
  efficacy.new.noswitch = sum(efficacy[!switch]*sexacts[!switch])/sum(sexacts[!switch])
  
  daily.q=as.numeric(quantile(efficacy.df$efficacy.daily,prob=c(.025,.5,.975), na.rm=T))
  time.q=as.numeric(quantile(efficacy.df$efficacy.time,prob=c(.025,.5,.975), na.rm=T))
  event.q=as.numeric(quantile(efficacy.df$efficacy.event,prob=c(.025,.5,.975), na.rm=T))
  r211.q=as.numeric(quantile(efficacy.df$efficacy.211,prob=c(.025,.5,.975), na.rm=T))
  
  c(choose.daily=choose.daily,
    choose.time=choose.time,
    choose.event=choose.event,
    choose.211=choose.211,
    samples=samples,
    efficacy.orig=efficacy.orig*100,
    efficacy.new=efficacy.new*100,
    pillrate.orig=pillrate.orig,
    pillrate.orig.sd=pillrate.orig.sd,
    pillrate.new=pillrate.new,
    pillrate.new.followup=pillrate.new.followup,
    pillrate.new.followup.sd=pillrate.new.followup.sd,
    efficacy.orig.followup=efficacy.orig.followup*100,
    efficacy.new.followup=efficacy.new.followup*100,
    pillrate.orig.followup=pillrate.orig.followup,
    pillrate.new.followup=pillrate.new.followup,
    pillrate.new.switch = pillrate.new.switch,
    pillrate.orig.switch = pillrate.orig.switch,
    efficacy.new.switch = efficacy.new.switch*100,
    efficacy.orig.switch = efficacy.orig.switch*100,
    efficacy.new.noswitch = efficacy.new.noswitch*100,
    efficacy.orig.noswitch = efficacy.orig.noswitch*100,
    efficacy.orig.personal = efficacy.orig.personal*100,
    efficacy.orig.personal.sd = efficacy.orig.personal.sd*100,
    efficacy.followup.personal = efficacy.followup.personal*100,
    efficacy.followup.personal.sd = efficacy.followup.personal.sd*100,
    daily.lo=daily.q[1]*100,
    daily.med=daily.q[2]*100,
    daily.hi=daily.q[3]*100,
    event.lo=event.q[1]*100,
    event.med=event.q[2]*100,
    event.hi=event.q[3]*100,
    time.lo=time.q[1]*100,
    time.med=time.q[2]*100,
    time.hi=time.q[3]*100,
    r211.lo=r211.q[1]*100,
    r211.med=r211.q[2]*100,
    r211.hi=r211.q[3]*100)
}
#assigns efficacy, sexacts, pillrate, ... based on best regime arm and calculates original and final efficacy & pillrate
#output: % in each regime, efficacy and pillrate originally and with best regime in plot.pills.effectiveness.change (in runsim)
#input: efficacy.df
