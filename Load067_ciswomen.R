
processdata=function(pillfile, participantdata){
  pilldata=pillfile$patient
  
  #Get arm of patient (1-daily, 2-time, 3-event)
  patient.arm=as.numeric(pilldata[2,,])
  
  #Time stamps of each pill
  patient.pilldays=pilldata[3,,]
  
  #Time stamps of each sexact
  patient.sexdaysA=pilldata[5,,] #AI
  
  #Whether or not condom is used AI
  patient.condomA=pilldata[6,,]
  
  #Time stamps of each sexact
  patient.sexdaysV=pilldata[8,,] #vaginal act
  
  #Whether or not condom is used VI
  patient.condomV=pilldata[9,,]
  
  N.patients=length(patient.arm)
  
  #Get the indices of the patients in each arm
  index.daily=which(patient.arm==1) #1-59
  index.time=which(patient.arm==2)  #120-179
  index.event=which(patient.arm==3) #60-119
  
  #make empty df/vector for data
  dat=as.data.frame(matrix(ncol=10,nrow=0)) # day by day data for each patient: all recorded days
  sexactdata=vector(N.patients,mode='list') # data for each sex act for each patient
  
  daily.sexpills=vector(N.patients,mode='list')
  ndays=numeric(N.patients)
  npills=ndays
  for(i in seq(N.patients)){
    pill = patient.pilldays[[i]][1,]
    npills[i] = length(pill)
    sexA = patient.sexdaysA[[i]][1,]
    condomA=patient.condomA[[i]][1,]
    sexV = patient.sexdaysV[[i]][1,]
    condomV=patient.condomV[[i]][1,]
    pilltime=round(patient.pilldays[[i]][1,])
    sextimeA=round(patient.sexdaysA[[i]][1,])
    sextimeV=round(patient.sexdaysV[[i]][1,])
    
    #In Dobromir's files, there is a quirk where a sextime of zero is given if 
    #there are no sex acts. A sex time of zero corresponds to some time
    #in 1970, so we have to just eliminate that act entirely
    if(sextimeA[1]==0){
      sextimeA=numeric(0)
      sexA=numeric(0)
      condomA=numeric(0)
    }
    
    if(sextimeV[1]==0){
      sextimeV=numeric(0)
      sexV=numeric(0)
      condomV=numeric(0)
    }
    
    sex=c(sexA,sexV)
    condom=c(condomA,condomV)
    daily.sexpills[[i]]=merge(makedaily(pill,sex),makedaily(pill,sexA),by='day',suffixes=c("",".a"))
    
    #Get the first and last times referenced in the pill and sex data
    first=floor(min(c(pill,sexA,sexV)))
    last=floor(max(c(pill,sexA, sexV)))
    
    #This gives the trial day on which pills and sex occur (relative to recorded start rather than absolute)
    adjpilltime=floor(pill)-first+1
    adjsextimeA=floor(sexA)-first+1
    adjsextimeV=floor(sexV)-first+1
    adjsextime=c(adjsextimeA,adjsextimeV)
    
    #Vectors of sex acts
    NsexA=length(sexA)
    NsexV=length(sexV)
    pre=numeric(NsexA+NsexV) #vector indicating precoverage #vector of 0s
    post=pre  #post coverage
    total=pre #total pills taken during the week
    pretime=pre
    posttime=pre
    efficacy=pre #efficacy
    coverage=pre #coverage level
    daybefore=pre
    prepill=pre #Was there a pill 2-48hrs prior to sex?
    anal=c(rep(T,NsexA),rep(F,NsexV))
    
    
    #Vectors of days #vectors of 0s
    havesexA=numeric(last-first+1) #number of times a person has anal sex on a given day
    havesexV=numeric(last-first+1) #number of times a person has anal sex on a given day
    takepill=numeric(last-first+1) #whether they take a pill
    precoverage =numeric(last-first+1)   #whether they have precoverage
    postcoverage = numeric(last-first+1) #whether they have postcoverage
    totalpill = numeric(last-first+1)    #how many pills in the 5 days before and 2 days after
    
    if(length(adjpilltime)>0){
      for(j in adjpilltime){
        takepill[j]=takepill[j]+1
      }
    }
    
    if(NsexA+NsexV>0){
      for(j in seq(NsexA+NsexV)){
        if(anal[j]==1){
          havesexA[adjsextime[j]]=havesexA[adjsextime[j]]+1
        }
        else{
          havesexV[adjsextime[j]]=havesexV[adjsextime[j]]+1
        }
        s=sex[j]
        rel=pill-s
        pre[j]=length(which(pill>=s-4&pill<s))>0    #pills 4-0 days before act
        post[j]=length(which(pill<=s+1&pill>=s))>0  #pills 0-1 days after act
        total[j]=length(which(pill<=s+2&pill>=s-5)) #pills 5 days before - 2 days after act 
        daybefore[j]=ifelse(adjsextime[j]>1,takepill[adjsextime[j]-1],0) #pill day before act
        posttime[j]=min(rel[rel>=0]) 
        pretime[j]=max(rel[rel<0])
        prepill[j]=length(which(pill>=s-2&pill<s-1/12))>0 #pill 2 days-1 hour before act
      }
    }
    
    havesex=havesexA+havesexV
    
    # for ts in length of time have data, calculate pills taken before/after each day
    for(j in seq(last-first+1)){
      t=first+j-1
      precoverage[j]=length(which(pill>t-4&pill<t))>0 
      postcoverage[j]=length(which(pill<t+1&pill>t))>0
      totalpill[j]=length(which(pill<t+2&pill>t-5))
    }
    #efficacy for daily use?
    coverage=pre+post 
    efficacy[which(total>3)]=fullefficacy #4 or more pills get the high efficacy value
    efficacy[which(total>1&total<4)]=partialefficacy #2-3 pills the middle efficacy value
    
    efficacy[which(daybefore&posttime<=1/12)]=eventefficacy #for event driven
    efficacy[which(total>6)]=maxefficacy
    
    days=seq(last-first+1)
    ndays[i]=last-first+1
    
    id=rep(i, length(days))
    arm=rep(patient.arm[i],length(days))
    dat.i=data.frame(day=days,pill=takepill,sex=havesex,sexA=havesexA,sexV=havesexV,arm=arm,id=id,precoverage=precoverage,postcoverage=postcoverage,totalpill=totalpill)
    
    sex.i=data.frame(day=adjsextime,
                     precoverage=pre,
                     postcoverage=post,
                     totalpills=total,
                     efficacy=efficacy,
                     coverage=coverage,
                     condom=condom,
                     anal=anal,
                     posttime=posttime,
                     pretime=pretime,
                     prepill=prepill)
    sexactdata[[i]]=sex.i
    dat=rbind(dat,dat.i)
  }
  
  participant.daily=subset(participantdata, RANarm=="Daily usage")
  ord.id.daily=order(participant.daily$pub_id)
  participant.time=subset(participantdata, RANarm=="Time-driven usage")
  ord.id.time=order(participant.time$pub_id)
  participant.event=subset(participantdata, RANarm=="Event-driven usage")
  ord.id.event=order(participant.event$pub_id)
  
  participant.df=rbind(participant.daily[ord.id.daily,],participant.time[ord.id.time,],participant.event[ord.id.event,])
  
  dat$presexacts=participant.df[dat$id,"SEXACT"] 
  dat$presexparts=participant.df[dat$id,"ACPARTN"]
  
  participantstats=as.data.frame(t(sapply(sexactdata,GetSexStats)))
  
  participantstats$numpills=npills
  participantstats$SEXACTSPRIOR=participant.df$SEXACT
  participantstats$PARTNERSPRIOR=participant.df$ACPARTN
  participantstats$ARM=patient.arm
  participantstats$DAYS=ndays
  participantstats$efficacy.new=sapply(sexactdata, getefficacy.fromdata)
  list(dat,sexactdata,participant.df,participantstats,daily.sexpills)
}
#pulls in individual data and calculates precoverage, postcoverage, sex acts, ...
#output: dat,sexactdata,participant.df,participantstats,daily.sexpills 
  #pd.H/B/C (below)
#input: pillfile, participantdata (pillfile.harlem,participant.dat.H,... from input data)

processdata2=function(pillfile, participantdata){
  pilldata=pillfile$patient
  
  #Get arm of patient (1-daily, 2-time, 3-event)
  patient.arm=as.numeric(pilldata[2,,])
  #patient.arm=patient.arm[-conc_mismatch_location]
  #Time stamps of each pill
  patient.pilldays=pilldata[3,,]
  #patient.pilldays=patient.pilldays[-conc_mismatch_location]
  #Time stamps of each sexact
  patient.sexdaysA=pilldata[5,,] #AI
  #patient.sexdaysA=patient.sexdaysA[-conc_mismatch_location]
  #Whether or not condom is used AI
  patient.condomA=pilldata[6,,]
  #patient.condomA=patient.condomA[-conc_mismatch_location]
  #Time stamps of each sexact
  patient.sexdaysV=pilldata[8,,] #vaginal act
  #patient.sexdaysV=patient.sexdaysV[-conc_mismatch_location]
  #Whether or not condom is used VI
  patient.condomV=pilldata[9,,]
  #patient.condomV=patient.condomV[-conc_mismatch_location]
  
  N.patients=length(patient.arm)
  
  #Get the indices of the patients in each arm
  index.daily=which(patient.arm==1) #1-59
  index.time=which(patient.arm==2)  #120-179
  index.event=which(patient.arm==3) #60-119
  
  #make empty df/vector for data
  dat=as.data.frame(matrix(ncol=10,nrow=0)) # day by day data for each patient: all recorded days
  sexactdata=vector(N.patients,mode='list') # data for each sex act for each patient
  
  daily.sexpills=vector(N.patients,mode='list')
  ndays=numeric(N.patients)
  npills=ndays
  for(i in seq(N.patients)){
    pill = patient.pilldays[[i]][1,]
    npills[i] = length(pill)
    sexA = patient.sexdaysA[[i]][1,]
    condomA=patient.condomA[[i]][1,]
    sexV = patient.sexdaysV[[i]][1,]
    condomV=patient.condomV[[i]][1,]
    pilltime=round(patient.pilldays[[i]][1,])
    sextimeA=round(patient.sexdaysA[[i]][1,])
    sextimeV=round(patient.sexdaysV[[i]][1,])
    
    #In Dobromir's files, there is a quirk where a sextime of zero is given if 
    #there are no sex acts. A sex time of zero corresponds to some time
    #in 1970, so we have to just eliminate that act entirely
    if(sextimeA[1]==0){
      sextimeA=numeric(0)
      sexA=numeric(0)
      condomA=numeric(0)
    }
    
    if(sextimeV[1]==0){
      sextimeV=numeric(0)
      sexV=numeric(0)
      condomV=numeric(0)
    }
    
    sex=c(sexA,sexV)
    condom=c(condomA,condomV)
    daily.sexpills[[i]]=merge(makedaily(pill,sex),makedaily(pill,sexA),by='day',suffixes=c("",".a"))
    
    #Get the first and last times referenced in the pill and sex data
    first=floor(min(c(pill,sexA,sexV)))
    last=floor(max(c(pill,sexA, sexV)))
    
    #This gives the trial day on which pills and sex occur (relative to recorded start rather than absolute)
    adjpilltime=floor(pill)-first+1
    adjsextimeA=floor(sexA)-first+1
    adjsextimeV=floor(sexV)-first+1
    adjsextime=c(adjsextimeA,adjsextimeV)
    
    #Vectors of sex acts
    NsexA=length(sexA)
    NsexV=length(sexV)
    pre=numeric(NsexA+NsexV) #vector indicating precoverage #vector of 0s
    post=pre  #post coverage
    total=pre #total pills taken during the week
    pretime=pre
    posttime=pre
    efficacy=pre #efficacy
    coverage=pre #coverage level
    daybefore=pre
    prepill=pre #Was there a pill 2-48hrs prior to sex?
    anal=c(rep(T,NsexA),rep(F,NsexV))
    
    
    #Vectors of days #vectors of 0s
    havesexA=numeric(last-first+1) #number of times a person has anal sex on a given day
    havesexV=numeric(last-first+1) #number of times a person has anal sex on a given day
    takepill=numeric(last-first+1) #whether they take a pill
    precoverage =numeric(last-first+1)   #whether they have precoverage
    postcoverage = numeric(last-first+1) #whether they have postcoverage
    totalpill = numeric(last-first+1)    #how many pills in the 5 days before and 2 days after
    
    if(length(adjpilltime)>0){
      for(j in adjpilltime){
        takepill[j]=takepill[j]+1
      }
    }
    
    if(NsexA+NsexV>0){
      for(j in seq(NsexA+NsexV)){
        if(anal[j]==1){
          havesexA[adjsextime[j]]=havesexA[adjsextime[j]]+1
        }
        else{
          havesexV[adjsextime[j]]=havesexV[adjsextime[j]]+1
        }
        s=sex[j]
        rel=pill-s
        pre[j]=length(which(pill>=s-4&pill<s))>0    #pills 4-0 days before act
        post[j]=length(which(pill<=s+1&pill>=s))>0  #pills 0-1 days after act
        total[j]=length(which(pill<=s+2&pill>=s-5)) #pills 5 days before - 2 days after act 
        daybefore[j]=ifelse(adjsextime[j]>1,takepill[adjsextime[j]-1],0) #pill day before act
        posttime[j]=min(rel[rel>=0]) 
        pretime[j]=max(rel[rel<0])
        prepill[j]=length(which(pill>=s-2&pill<s-1/12))>0 #pill 2 days-1 hour before act
      }
    }
    
    havesex=havesexA+havesexV
    
    # for ts in length of time have data, calculate pills taken before/after each day
    for(j in seq(last-first+1)){
      t=first+j-1
      precoverage[j]=length(which(pill>t-4&pill<t))>0 
      postcoverage[j]=length(which(pill<t+1&pill>t))>0
      totalpill[j]=length(which(pill<t+2&pill>t-5))
    }
    #efficacy for daily use?
    coverage=pre+post 
    efficacy[which(total>3)]=fullefficacy #4 or more pills get the high efficacy value
    efficacy[which(total>1&total<4)]=partialefficacy #2-3 pills the middle efficacy value
    
    efficacy[which(daybefore&posttime<=1/12)]=fullefficacy #for event driven
    efficacy[which(total>6)]=maxefficacy
    
    days=seq(last-first+1)
    ndays[i]=last-first+1
    
    id=rep(i, length(days))
    arm=rep(patient.arm[i],length(days))
    dat.i=data.frame(day=days,pill=takepill,sex=havesex,sexA=havesexA,sexV=havesexV,arm=arm,id=id,precoverage=precoverage,postcoverage=postcoverage,totalpill=totalpill)
    
    sex.i=data.frame(day=adjsextime,
                     precoverage=pre,
                     postcoverage=post,
                     totalpills=total,
                     efficacy=efficacy,
                     coverage=coverage,
                     condom=condom,
                     anal=anal,
                     posttime=posttime,
                     pretime=pretime,
                     prepill=prepill)
    sexactdata[[i]]=sex.i
    dat=rbind(dat,dat.i)
  }
  
  participant.daily=subset(participantdata, RANarm=="Daily usage")
  ord.id.daily=order(participant.daily$pub_id)
  participant.time=subset(participantdata, RANarm=="Time-driven usage")
  ord.id.time=order(participant.time$pub_id)
  participant.event=subset(participantdata, RANarm=="Event-driven usage")
  ord.id.event=order(participant.event$pub_id)
  
  participant.df=rbind(participant.daily[ord.id.daily,],participant.time[ord.id.time,],participant.event[ord.id.event,])
  
  dat$presexacts=participant.df[dat$id,"SEXACT"] 
  dat$presexparts=participant.df[dat$id,"ACPARTN"]
  
  participantstats=as.data.frame(t(sapply(sexactdata,GetSexStats)))
  
  participantstats$numpills=npills
  participantstats$SEXACTSPRIOR=participant.df$SEXACT
  participantstats$PARTNERSPRIOR=participant.df$ACPARTN
  participantstats$ARM=patient.arm
  participantstats$DAYS=ndays
  participantstats$efficacy.new=sapply(sexactdata, getefficacy.fromdata)
  list(dat,sexactdata,participant.df,participantstats,daily.sexpills)
}
#removes elements from lists to eliminate individuals with drug concentrations =/= reported pills taken
#pulls in individual data and calculates precoverage, postcoverage, sex acts, ...
#output: dat,sexactdata,participant.df,participantstats,daily.sexpills 
#pd.H/B/C (below)
#input: pillfile, participantdata (pillfile.harlem,participant.dat.H,... from input data)

GetSexStats=function(X){
  nocoverage=sum(X$coverage==0)
  partialcoverage=sum(X$coverage==1)
  fullcoverage=sum(X$coverage==2)
  partialefficacy=sum(X$efficacy[X$coverage==1])
  fullefficacy=sum(X$efficacy[X$coverage==2])
  condomuse=sum(X$condom)
  condomusepart=sum(X$condom[X$coverage==1])
  condomusefull=sum(X$condom[X$coverage==2])
  pillcondomefficacypart=sum(X$efficacy[X$coverage==1])
  pillcondomefficacyfull=sum(X$efficacy[X$coverage==2])
  vaginal=sum(X$anal==0)
  efficacy=(partialefficacy+fullefficacy)/(nocoverage+fullcoverage+partialcoverage)
  c(none=nocoverage,partiallycoveredacts=partialcoverage,fullycoveredacts=fullcoverage,
    partiallycoveredefficacy=partialefficacy,fullycoveredefficacy=fullefficacy,condomuse=condomuse,
    condompart=condomusepart,condomfull=condomusefull,pillcondomefficacypart=pillcondomefficacypart,
    pillcondomefficacyfull=pillcondomefficacyfull, numacts=nocoverage+partialcoverage+fullcoverage, 
    efficacy=efficacy, numvaginal=vaginal)
} 
#calculates pill coverage levels 
#output: participantstats in processdata (above)
#input: sexactdata from processdata (above)

makedaily=function(pills,sex){
  comboppre=numeric(0)
  combospre=numeric(0)
  comboppost=numeric(0)
  combospost=numeric(0)
  combopprepost=numeric(0)
  combosprepost=numeric(0)
  pill24=numeric(0)
  pill48=numeric(0)
  pill72=numeric(0)
  sex24=numeric(0)
  sex48=numeric(0)
  sex72=numeric(0)
  pillaftersex=numeric(0)
  for(i in seq_along(sex)){
    s=sex[i]
    post=pills>=s&pills<=(s+1/12) #within 2 hours after sex
    pre=pills<s&pills>=(s-1/12)   #within 2 hours before sex
    if(sum(post)>0&sum(pre)==0){
      p=which(post)[1]
      comboppost=c(comboppost,p)
      combospost=c(combospost,i)
    }
    if(sum(post)==0&sum(pre)>0){
      p=which(pre)[1]
      comboppre=c(comboppre,p)
      combospre=c(combospre,i)
    }
    if(sum(post)>0&sum(pre)>0){
      p1=which(pre)[1]
      p2=which(post)[1]
      combopprepost=c(combopprepost,p1,p2)
      combosprepost=c(combosprepost,i)
    }
    pre2=pills<=s
    if(sum(pre2)>0){ #pill 24/48/72 hrs before sex
      p=max(pills[pre2]) #last pill before sex
      diff=s-p           #how long between pill & sex
      if(diff<1){
        sex24=c(sex24,s)
      }
      else if(diff<2){
        sex48=c(sex48,s)
      }
      else if(diff<3){
        sex72=c(sex72,s)
      }
    }
  }
  
  combopre=sex[combospre]
  combopost=sex[combospost]
  comboprepost=sex[combosprepost]
  combo=c(combopre,combopost,comboprepost) #pills within 2 hours of sex
  
  sexonly=setdiff(sex,combo)
  pillonly=setdiff(pills,pills[c(comboppre,comboppost,combopprepost)])
  
  for(i in seq_along(pillonly)){ 
    p=pillonly[i]
    post=sex>=p 
    pre=sex<p
    if(sum(post)>0){
      s=sex[which(post)[1]]
      diff=s-p 
      if(diff<1){
        pill24=c(pill24,p)
      }
      else if(diff<2){
        pill48=c(pill48,p)
      }
      else if(diff<3){
        pill72=c(pill72,p)
      }
      if(sum(pre)>0&s>=ceiling(p)){
        s2=sex[which(pre)[1]]
        if(s2>=p-1){
          pillaftersex=c(pillaftersex,p)
        }
      }
    }
  }
  
  pillother=setdiff(pillonly,c(pill24,pill48,pill72))
  first=floor(min(c(pills,sex)))
  last=floor(max(c(pills,sex)))
  
  adj.f=function(x){floor(x)-first+1}
  
  Ndays=last-first+1
  combopost=adj.f(combopost)
  combopre=adj.f(combopre)
  comboprepost=adj.f(comboprepost)
  pillonly=adj.f(pillonly)
  sexonly=adj.f(sexonly)
  pillother=adj.f(pillother)
  pill24=adj.f(pill24)
  pill48=adj.f(pill48)
  pill72=adj.f(pill72)
  pillaftersex=adj.f(pillaftersex)
  sex24=adj.f(sex24)
  sex48=adj.f(sex48)
  sex72=adj.f(sex72)
  
  sexonly.vec=numeric(Ndays)
  pillonly.vec=numeric(Ndays)
  pill24.vec=numeric(Ndays)
  pill48.vec=numeric(Ndays)
  pill72.vec=numeric(Ndays)
  pillaftersex.vec=numeric(Ndays)
  sex24.vec=numeric(Ndays)
  sex48.vec=numeric(Ndays)
  sex72.vec=numeric(Ndays)
  pillother.vec=numeric(Ndays)
  combopre.vec=numeric(Ndays)
  combopost.vec=numeric(Ndays)
  comboprepost.vec=numeric(Ndays)
  day=seq(Ndays)
  
  for(i in day){
    sexonly.vec[i]=sum(sexonly==i)
    pillonly.vec[i]=sum(pillonly==i)
    combopre.vec[i]=sum(combopre==i)
    combopost.vec[i]=sum(combopost==i)
    comboprepost.vec[i]=sum(comboprepost==i)
    pill24.vec[i]=sum(pill24==i)
    pill48.vec[i]=sum(pill48==i)
    pill72.vec[i]=sum(pill72==i)
    sex24.vec[i]=sum(sex24==i) #pills 24 hours before sex
    sex48.vec[i]=sum(sex48==i) #pills 48 hours before sex
    sex72.vec[i]=sum(sex72==i) #pills 72 hours before sex
    pillother.vec[i]=sum(pillother==i)
    pillaftersex.vec[i]=sum(pillaftersex==i)
  }
  
  data.frame(day=day,sexonly=sexonly.vec,pillonly=pillonly.vec,combopre=combopre.vec,
             combopost=combopost.vec, comboprepost=comboprepost.vec,
             pill24=pill24.vec,pill48=pill48.vec,pill72=pill72.vec,
             sex24=sex24.vec,sex48=sex48.vec,sex72=sex72.vec,
             pillother=pillother.vec,pillaftersex=pillaftersex.vec)
}
#calculates when pills are taken within 2/24/48/72 hours of sex or days with sex only or pill only
#output: daily.sexpills in processdata (above)
#input: pill,sex/sexA where pill = patient.pilldays[[i]][1,], sexA = patient.sexdaysA[[i]][1,] sexV = patient.sexdaysV[[i]][1,] ,sex=sexA+sexV

getefficacy.fromdata=function(sexactdata,precutoff=3,postcutoff=1/12){
  event.post=ifelse(sexactdata$posttime<=precutoff,TRUE,FALSE) #why reversed?
  event.pre=ifelse(sexactdata$pretime>=-postcutoff,TRUE,FALSE)
  
  event.covered=event.pre&event.post
  
  efficacy=sexactdata$efficacy
  
  efficacy[event.covered]=pmax(efficacy[event.covered],eventefficacy)
  mean(efficacy)
}
#if there were pills <=3 days after and >= 2 hours before sex, makes efficacy=eventefficacy
#  outputs mean efficacy with the new efficacy numbers
#output: participantstats$efficacy.new in processdata (above)
#input: sexactdata from processdata (above)

getadherence.daily=function(dailysex){
  success=sum(dailysex$pillonly+dailysex$combopre+dailysex$combopost+dailysex$comboprepost>0)
  failure=sum(dailysex$pillonly+dailysex$combopre+dailysex$combopost+dailysex$comboprepost==0)
  
  success/(success+failure)
} 
#determines proportion of time there were any pills taken around a sex act
#output: ad.daily, participantstats.harlem$adherence.regimen (below)
#input: dailysex.harlem[ARM.harlem==1] (dailysex.harlem=pd.H[[5]])

getadherence.time=function(dailysex){
  N.days=length(dailysex$pillonly)
  N.weeks=floor(N.days/7)
  N.atleasttwo=0
  for(i in seq(N.weeks)){
    j=seq(7)+(i-1)*7
    n.pills=sum(dailysex$pillonly[j]+dailysex$combopre[j]+dailysex$combopost[j]+dailysex$comboprepost[j])
    N.atleasttwo=N.atleasttwo+1*(n.pills>1)
  }
  adherence.regular=N.atleasttwo/N.weeks
  adherence.sex = sum(dailysex$combopost+dailysex$comboprepost)/sum(dailysex$combopost+dailysex$comboprepost
                                                                    + dailysex$sexonly+dailysex$combopre)

  
  (adherence.regular + adherence.sex)/2
}
#proportion weeks had at least 2 pills taken and sex acts with post pills taken
#output: ad.time, participantstats.harlem$adherence.regimen (below)
#input: dailysex.harlem[ARM.harlem==2] (dailysex.harlem=pd.H[[5]])

getadherence.event=function(dailysex){
  pills=dailysex$combopost+dailysex$combopre+dailysex$comboprepost+dailysex$pillonly
  sex=dailysex$combopost+dailysex$combopre+dailysex$comboprepost+dailysex$sexonly
  N.days=length(pills)
  success=0
  failure=0
  for(i in seq(1,N.days)){
    if(i>1){
    if(pills[i]+pills[i-1]>0){
      success=success+sex[i] #took pill around sex day of or day before AND had sex = success
    }
    else{
      failure=failure+sex[i] #no pill around sex day of or day before AND had sex = failure
    }}
    
    success=success+dailysex$combopost[i]+dailysex$comboprepost[i]
    failure=failure+dailysex$sexonly[i]+dailysex$combopre[i]
  }
  success/(success+failure)
}
#calculates proportion days pills were taken before/after sex
#output: ad.event, participantstats.harlem$adherence.regimen (below)
#input: dailysex.harlem[ARM.harlem==3] (dailysex.harlem=pd.H[[5]])

getadherence.211=function(dailysex){
  pills=dailysex$combopost+dailysex$combopre+dailysex$comboprepost+dailysex$pillonly#
  sex=dailysex$combopost+dailysex$combopre+dailysex$comboprepost+dailysex$sexonly
  N.days=length(pills)
  success=0
  failure=0
  for(i in seq(1,N.days)){
    if(i>1){
      if(pills[i]+pills[i-1]>0){ 
        success=success+sex[i] #took pill around sex day of or day before AND had sex = success
      }
      else{
        failure=failure+sex[i] #no pill around sex day of or day before AND had sex = failure
      }}
    
    success=success+dailysex$combopost[i]+dailysex$comboprepost[i]
    failure=failure+dailysex$sexonly[i]+dailysex$combopre[i]
  }
  success/(success+failure)
}
#calculates proportion days pills were taken before/after sex
#output: ad.211, participantstats.harlem$adherence.regimen (below)
#input: dailysex.harlem[ARM.harlem==3] (dailysex.harlem=pd.H[[5]])

evaluate.specificity = function(dailysex){
  sex = dailysex$sexonly+dailysex$combopre+dailysex$combopost+dailysex$comboprepost
  days.no.sex = which(sex==0)
  predays.no.sex = pmax(1,days.no.sex-1)
  nopill.dayof = dailysex$pillonly[days.no.sex]==0
  nopill.daybefore = dailysex$pillonly[predays.no.sex]==0
  sex.daybefore = sex[predays.no.sex] > 0
  mean(nopill.dayof & (nopill.daybefore|sex.daybefore))
}
#specificity of when pills are taken relative to sex acts in event-based regimes
#output: spec.event (below)
#input: dailysex.capetown[ARM.capetown==3]

evaluate.sensitivity = function(dailysex){
  sex = dailysex$sexonly+dailysex$combopre+dailysex$combopost+dailysex$comboprepost
  pills = dailysex$pillonly
  combopills = dailysex$combopre+dailysex$combopost+dailysex$comboprepost
  days.sex = which(sex>0)
  predays.sex = pmax(1,days.sex-1)
  pill.needed = days.sex[pills[days.sex]==0&pills[predays.sex]==0&combopills[predays.sex]==0]
  pill.not.needed = setdiff(days.sex,pill.needed)
  
  needed.pills.taken = (dailysex$combopre[pill.needed] + dailysex$comboprepost[pill.needed])>0
  sum(sex[pill.not.needed])/(sum(needed.pills.taken) + sum(sex[pill.not.needed]))
}
#sensitivity of when pills are taken relative to sex acts in event-based regimes
#output: sens.event (below)
#input: dailysex.capetown[ARM.capetown==3]

get.efficacy.trial=function(df){
  eff=sum(df$efficacy*df$numacts, na.rm = T)/sum(df$numacts)
  round(100*eff)
}
#sum of efficacy for each act/ number acts
#output: eff.daily, eff.time, eff.event (below)
#input: participantstats.capetown, ARM==1/2/3

#`%notin%` <- function(x,y) !(x %in% y) 
`%notin%` <- Negate(`%in%`)

