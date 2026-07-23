#Autocorrelated binomial

simulate.sex=function(N.days,sex.params){
  with(as.list(sex.params),{
    s=exp(sexfrequency)
    a=analfrequency
    ca=expit(condomanal)
    cv=expit(condomvaginal)
    
    S=rpois(N.days,lambda=s)
    boundaryacts=sum(S[c(seq(5),N.days-1,N.days)])
    if(sum(S)==boundaryacts){
      mainindex=round(N.days/2)
      S[mainindex]=1
    }
    S
  })
}
#simulates sex acts based on a poisson distribution with labda=exp(sexfrequency)
#output: S in calc.efficacy.daily.rand, calc.efficacy.time, calc.efficacy.event (below)
#input: sex.params from adherence.params

calc.efficacy=function(weeklytotal,S,S.post,Needed, p.event=0){
  if(sum(S)==0){
    return(NA)
  }
  PrePill = !Needed
  
  eventcoverage=sum(S.post[(weeklytotal<4&PrePill)])
  coverage.7=sum(S[weeklytotal>=7])/sum(S)
  coverage.4=sum(S[weeklytotal>=4&(weeklytotal<7)])/sum(S)
  coverage.2=sum(S[weeklytotal>=2&(weeklytotal<4)])/sum(S)-sum(eventcoverage)/sum(S)
  coverage.event=sum(eventcoverage)/sum(S)
  efficacy = maxefficacy*coverage.7+coverage.4*fullefficacy+coverage.2*partialefficacy+coverage.event*eventefficacy#+coverage.211*r211efficacy
  
  if(efficacy == 0){
    p.7 = mean(weeklytotal>7)
    p.4 = mean(weeklytotal%in%c(4,5,6))
    p.2 = mean(weeklytotal%in%c(2,3))
    efficacy = maxefficacy*p.7+p.4*fullefficacy+pmax(0, p.2 - p.event)*partialefficacy+p.event*eventefficacy
  }
  efficacy
}
#calculates efficacy based on number of pills taken in a week (99% for 7 pills, ...)
#output: efficacy in calc.efficacy.daily.rand, calc.efficacy.time, calc.efficacy.event (below)
#input: weeklytotal,S,S.post,Needed, p.event=0 (within fn called in below)

calc.efficacy.daily.rand=function(adherence.params,N.days){
  with(as.list(c(adherence.params)),{
    #Probability of pill if no pill day before
    p0=expit(pilldaybeforeFALSE)
    
    #Probability of pill if pill day before
    p1=expit(pilldaybeforeTRUE)
    
    S=simulate.sex(N.days,adherence.params)
    P=numeric(N.days+1)
    P[1]=1 #Assume pill taken on day one
    for(i in seq(N.days)){
      p=ifelse(P[i]>0,p1,p0)
      P[i+1]=P[i+1]+rbinom(1,size=1,prob=p)
    }
    
    #Fraction of pills taken within 2 hour window
    g.post=2/24
    
    P=P[-1]
    
    S.post=rbinom(N.days,size=1,prob=(1-(1-g.post)^S)*P)
    
    #E=S>0|c(S[-1],0)>0
    #M=c(0,S.nopost[-N.days])>0
    #N=!E&!M
    
    #pillprob=ifelse(E,p0,ifelse(M,p1,p3))
    #P=rbinom(N.days,size=1,prob=pillprob)
    
    totalpills=P
    
    runningtotal=cumsum(c(totalpills,0,0))
    
    weeklytotal=runningtotal[3:(N.days+2)]-c(rep(0,5),runningtotal[1:(N.days-5)])
    Needed=c(F,P[-N.days]==0)&(P==0)
    
    #Don't count acts near boundary
    S.post[c(seq(5),N.days-1,N.days)]=0
    S[c(seq(5),N.days-1,N.days)]=0
    
    #Probability of pill overall
    p=as.numeric(p0/(1+p0-p1))
    
    adherence=mean(P>0)
    
    adherence=ifelse(adherence>0&adherence<1, adherence, p)
    
    efficacy=calc.efficacy(weeklytotal,S,S.post,Needed, p.event = p*p1/12)
    c(efficacy=efficacy, logitefficacy = logit(efficacy),pillrate=mean(P),sensitivity=NA,specificity=NA,adherence=adherence,adherence.regular = adherence, adherence.sex = NA, sexacts=sum(S), logitadherence = logit(adherence), logitadherence.regular = logit(adherence), logitadherence.sex = NA)
  })
}
#simulates sex acts to find efficacy based on adherence params
#output: efficacy.daily in get.efficacy (in Crossoversimulationcode.R)
#calls: expit, simulate.sex, calc.efficacy
#input: params.daily from params.067/params.rand

calc.efficacy.time=function(adherence.params, N.days=10000){
  
  with(as.list(adherence.params),{
    #Probability of needed regular pill
    p0=expit(timetypeneeded)
    
    #sexfrequency
    s=exp(sexfrequency)
    ps=1-exp(-s)
    
    #Probability of pillday after sex
    p1=expit(timetypemakeup)
    
    #Probability of pill two hours after sex
    p2=expit(timetypepost)
    
    #Probability of pill when not needed
    p3=expit(timetypenotneeded)
    
    p4=expit(timetypepre)
    
    S=simulate.sex(N.days,adherence.params)
    S.post=pmin(rbinom(N.days,size=S,prob=p2),1)
    S.both=rbinom(N.days,size=S.post,prob=p4)
    S.postonly=S.post-S.both
    S.nopost=S-S.post
    S.preonly=pmin(rbinom(N.days,size=S.nopost,prob=p4),1)
    
    
    P=c(0,0,S.postonly+S.preonly+2*S.both)
    for(i in seq(N.days)){
      p=ifelse(P[i]+P[i+1]==0,p0,ifelse(S.nopost[i]>0,p1,p3))
      P[i+2]=P[i+2]+rbinom(1,size=1,prob=p)
    }
    
    P=P[-c(1,2)]
    
    Needed=c(T,P[-N.days]==0)&(P==0)
    
    totalpills=P
    
    runningtotal=cumsum(c(totalpills,0,0))
    
    weeklytotal=runningtotal[3:(N.days+2)]-c(rep(0,5),runningtotal[1:(N.days-5)])
    
    #Don't count acts near boundary
    S.post[c(seq(5),N.days-1,N.days)]=0
    S[c(seq(5),N.days-1,N.days)]=0
    
    pills.per.week=weeklytotal[seq(7,N.days,by=7)]
    
    adherence.1 = sum(pills.per.week>1)/(length(pills.per.week))
    adherence.2 = p2
    
    adherence.full = (adherence.1 + adherence.2)/2
    efficacy=calc.efficacy(weeklytotal,S,S.post, Needed, p.event = p2*p4)
    c(efficacy = efficacy, logitefficacy = logit(efficacy), pillrate=mean(totalpills),sensitivity=NA,specificity=NA,adherence=adherence.full, adherence.regular = adherence.1, adherence.sex = adherence.2, sexacts=sum(S), logitadherence = logit(adherence.full), logitadherence.regular = logit(adherence.1), logitadherence.sex = logit(adherence.2))
  })
}
#simulates sex acts to find efficacy based on adherence params
#output: efficacy.time in get.efficacy (in Crossoversimulationcode.R)
#calls: expit, simulate.sex, calc.efficacy
#input: params.time from params.067/params.rand

calc.efficacy.event=function(adherence.params,N.days=100000){
  with(as.list(adherence.params),{
    #Probability of pill in advance of sex (sex day of or day after)
    p0=expit(typepredict)
    
    #sexfrequency
    s=exp(sexfrequency)
    ps=1-exp(-s)
    
    #Probability of pillday after sex (sex and no post-sex pill day before, no sex day of or day after)
    p1=expit(typemakeup)
    
    #Probability of pill two hours after sex
    p2=expit(typepost)
    
    #Probability of pill not in advance of sex (no sex day before, day of, day after)
    p3=expit(typenotpredict)
    
    #Probability of pill two hours before sex (when needed)
    p4=expit(typeneededpre)
    
    #Probability of pill two hours before sex (when not needed)
    p5=expit(typeextrapre)
    
    S=simulate.sex(N.days,adherence.params)
    
    S.post=pmin(rbinom(N.days,size=S,prob=p2),1)
    S.nopost=S-S.post
    
    E=S>0|c(S[-1],0)>0
    M=c(0,S.nopost[-N.days])>0
    
    pillprob=ifelse(E,p0,ifelse(M,p1,p3))
    P=rbinom(N.days,size=1,prob=pillprob)
    Needed=c(T,P[-N.days]==0)&(P==0)
    
    pill.pre=numeric(N.days)
    pill.pre[!Needed]=pmin(rbinom(sum(!Needed),size=S[!Needed],prob=p5),1)
    
    pill.pre[Needed]=pmin(rbinom(sum(Needed),size=S[Needed],prob=p4),1)
    
    totalpills=S.post+P+pill.pre
    
    runningtotal=cumsum(c(totalpills,0,0))
    
    weeklytotal=runningtotal[3:(N.days+2)]-c(rep(0,5),runningtotal[1:(N.days-5)])
    
    #Don't count acts near boundary
    S.post[c(seq(5),N.days-1,N.days)]=0
    pill.pre[c(seq(5),N.days-1,N.days)]=0
    S[c(seq(5),N.days-1,N.days)]=0
    
    adherence.1 = (sum(S[!Needed])+sum(pill.pre[Needed]))/sum(S)
    adherence.2 = p2
    
    adherence.1 = ifelse(sum(S)>=10&(sum(S[!Needed])+sum(pill.pre[Needed]))>0&((sum(S[!Needed])+sum(pill.pre[Needed]))<sum(S)),adherence.1, 1-(1-p0)^2*(1-p4))
    #if: sum(S)>=10 & [0 < (sum(S[!Needed])+sum(pill.pre[Needed])) < sum(S)] then: adherence.1
    # else: 1-(1-[pill in advance of sex])^2*(1-[pill 2 hours before sex])
    #to avoid very small numbers of sex acts--> use theoretical number if number sex acts too small
    adherence = adherence.1/2 + adherence.2/2
    
    prepillstaken=(sum(pill.pre[Needed])+sum(S[!Needed]))
    theoreticalsensitivity = (1 - (1 - p0)^2)/((1 - (1 - p0)^2) + (1 - p0)^2*p4)
    
    theoreticalspecificity = (1 - p3)^2
    
    sensitivity=ifelse(prepillstaken>0&sum(S[!Needed])>0&sum(S[!Needed])<prepillstaken,(sum(S[!Needed]))/prepillstaken,theoreticalsensitivity)
    specificity=1-mean(P[!E&!M])
    specificity=ifelse(sum(!E&!M)>0&mean(P[!E&!M])>0&mean(P[!E&!M])<1,specificity, theoreticalspecificity)
    efficacy=calc.efficacy(weeklytotal,S,S.post, Needed, p.event = adherence.1*adherence.2)
    c(efficacy = efficacy, logitefficacy = logit(efficacy),pillrate=mean(totalpills),sensitivity=sensitivity,specificity=specificity,adherence=adherence,adherence.regular = adherence.1, adherence.sex = adherence.2, sexacts=sum(S), logitadherence = logit(adherence), logitadherence.regular = logit(adherence.1), logitadherence.sex = logit(adherence.2), logitsensitivity = logit(sensitivity))
  })
}
#simulates sex acts to find efficacy based on adherence params
#output: efficacy.event in get.efficacy (in Crossoversimulationcode.R)
#calls: expit, simulate.sex, calc.efficacy
#input: params.event from params.067/params.rand

calc.efficacy.211=function(adherence.params,N.days=100000,r211pill2=0.9){
  with(as.list(adherence.params),{
    #Probability of pill in advance of sex (sex day of or day after)
    p0=expit(typepredict)
    
    #sexfrequency
    s=exp(sexfrequency)
    ps=1-exp(-s)
    
    #Probability of pill day after sex (sex and no post-sex pill day before, no sex day of or day after)
    p1=expit(typemakeup)
    
    #Probability of pill two (24) hours after sex 
    p2=expit(typepost)
    
    #Probability of pill not in advance of sex (no sex day before, day of, day after)
    p3=expit(typenotpredict)
    
    #Probability of pill two (24) hours before sex (when needed)
    p4=expit(typeneededpre)
    
    #Probability of pill two (24) hours before sex (when not needed)
    p5=expit(typeextrapre)
    
    S=simulate.sex(N.days,adherence.params)
    
    S.post=pmin(rbinom(N.days,size=S,prob=p2),1) #post sex pill, only possible if simulated sex that day
    S.nopost=S-S.post
    
    E=S>0|c(S[-1],0)>0 #sex today or tomorrow
    M=c(0,S.nopost[-N.days])>0 #no post pill & sex day before 
    
    pillprob=ifelse(E,p0,ifelse(M,p1,p3)) #sex & pill, sex & makeup pill, pill not needed 
    P=rbinom(N.days,size=1,prob=pillprob) #pill taken
    
    Needed=c(T,P[-N.days]==0)&(P==0)      #pill 0-2 hours prior to sex, no other pre pill 2-48 hours before
    
    pill.pre=numeric(N.days)
    pill.pre[!Needed]=pmin(rbinom(sum(!Needed),size=S[!Needed],prob=p5),1) #not needed pill 2 hours before sex
    pill.pre[Needed]=pmin(rbinom(sum(Needed),size=S[Needed],prob=p4),1) #needed pill 2 hours before sex 
    
    #extra pills for 211
    nopill=(P+pill.pre+S.post)==0
    #second after pill
    S.post2a=pmin(rbinom(N.days,size=S.post,prob=r211pill2),1) #+ (S.post==0)*(background pill rate) #2nd post pill in 211 (if postpill 1 taken, postpill2 has r211pill2=.9 prob)
    S.post2b=c(0,S.post2a[1:(N.days-1)]) #post pills next day
    S.post2=S.post2b*nopill #only take pill if have not already taken one that day
    #double first pill
    totalpills2=P+pill.pre+S.post+S.post2
    totalpills.lastweek=c(0,totalpills2[1],sum(totalpills2[1:2]),sum(totalpills2[1:3]),sum(totalpills2[1:4]),
                           sum(totalpills2[1:5]),sum(totalpills2[1:6]),rollapply(totalpills2[1:(N.days-1)], 7, sum))
    nopill.lastweek=(totalpills.lastweek==0)*(totalpills2<=1) #no pills last week, 1 or 0 pills today (to prevent double pill when already taking 2+ on a given day)
  
    pillnum.P =ifelse(E&nopill.lastweek,2,ifelse(M&nopill.lastweek,2,1))    #211: double pill in advance of sex & as makeup, single pill when not needed
    pillnum.pre =ifelse(Needed&nopill.lastweek,2,1)
    
    totalpills=(P*pillnum.P)+(pill.pre*pillnum.pre)+S.post+S.post2 #double prepills  
    runningtotal=cumsum(c(totalpills,0,0))
    weeklytotal=runningtotal[3:(N.days+2)]-c(rep(0,5),runningtotal[1:(N.days-5)]) #pills 4 days before, 2 days after today
    
    #r211 total pills = 7 days before to 2 days after (don't double first pill unless no pills 7 days before, so have 10 days to get to 4 pills total)
    #r211total=c(sum(totalpills[1:3]),sum(totalpills[1:4]),sum(totalpills[1:5]),sum(totalpills[1:6]),sum(totalpills[1:7]),
    #            sum(totalpills[1:8]),sum(totalpills[1:9]),rollapply(totalpills, 10, sum),sum(totalpills[(N.days-8):N.days]),sum(totalpills[(N.days-7):N.days]))
    #r211adherent=r211total>3
    
    #r211 pre- and post- sex pills (need 2 pre/day of pills & 2 post pills)
    r211pre=c(sum(totalpills[1]),sum(totalpills[1:2]),sum(totalpills[1:3]),sum(totalpills[1:4]),sum(totalpills[1:5]),
                sum(totalpills[1:6]),sum(totalpills[1:7]),rollapply(totalpills, 8, sum))
    r211pre.adherent=r211pre>1
    r211post=c(rollapply(totalpills[2:N.days], 2, sum),totalpills[N.days],0)
    r211post.adherent=r211post>1
    
    #Don't count acts near boundary
    S.post[c(seq(5),N.days-1,N.days)]=0
    S.post2[c(seq(5),N.days-1,N.days)]=0
    pill.pre[c(seq(5),N.days-1,N.days)]=0
    S[c(seq(5),N.days-1,N.days)]=0
    
    #adherence = (sum(S[r211adherent]))/sum(S) #took at least 4 pills -7/+2 days around sex  
    adherence.1 = (sum(S[r211pre.adherent]))/sum(S) #took at least 2 prepills before/day of sex
    adherence.2 = (sum(S[r211post.adherent]))/sum(S)#took at least 2 postpills after sex
    adherence = adherence.1/2 + adherence.2/2
    
    prepillstaken=sum(S[r211pre.adherent]) #covered acts not covered pills
    
    theoreticalsensitivity = (1 - (1 - p0)^2)/((1 - (1 - p0)^2) + (1 - p0)^2*p4)
    theoreticalspecificity = (1 - p3)^2
    
    sensitivity=ifelse(prepillstaken>0&sum(S[!Needed])>0&sum(S[!Needed])<prepillstaken,(sum(S[!Needed]))/prepillstaken,theoreticalsensitivity)
    specificity=1-mean(P[!E&!M])
    specificity=ifelse(sum(!E&!M)>0&mean(P[!E&!M])>0&mean(P[!E&!M])<1,specificity, theoreticalspecificity)
    efficacy=calc.efficacy(weeklytotal,S,S.post, Needed, p.event = adherence.1*adherence.2)
    c(efficacy = efficacy, logitefficacy = logit(efficacy),pillrate=mean(totalpills),sensitivity=sensitivity,
      specificity=specificity,adherence=adherence,adherence.regular = adherence.1, adherence.sex = adherence.2, 
      sexacts=sum(S), logitadherence = logit(adherence), logitadherence.regular = logit(adherence.1), 
      logitadherence.sex = logit(adherence.2), logitsensitivity = logit(sensitivity))
  })
}  
#simulates sex acts to find efficacy based on adherence params
#output: efficacy.211 in get.efficacy, efficacy.211 & efficacy.211.followup in reassign.efficacy.byclass (in Crossoversimulationcode.R)
#calls: expit, simulate.sex, calc.efficacy
#input: params.211 from params.067/params.rand (from get.HPTN067.params(model.list)/get.random.params in  HPTN067GLMMR.R) 

calc.efficacy.daily=function(adherence.params){
  with(as.list(adherence.params),{
    #Probability of pill if no pill day before
    p0=expit(pilldaybeforeFALSE)
    
    #Probability of pill if pill day before
    p1=expit(pilldaybeforeTRUE)
    
    #Probability of pill overall
    p=as.numeric(p0/(1+p0-p1))
    
    #Days when pre pill must be taken
    predays=c(4)
    
    #Day when post pill must be taken
    dayof=5
    
    #Fraction of pills taken within 2 hour window
    g.post=2/24
    
    #Probability vectors
    P0=numeric(8) #No pill
    P1=numeric(8) #pill
    Q0=P0 #No pill and precoverage
    Q1=P0 #Pill and precoverage
    C0=P0 #Covered by event
    C1=P1 #Covered by event
    
    newP0=P0
    newP1=P1
    
    newQ0=Q0
    newQ1=Q1
    
    newC0=C0
    newC1=C1
    
    
    P0[1]=1-p
    P0[2]=p
    
    for(i in seq(2,7)){
      if(i %in% predays){ #=4
        for(j in seq(7)){
          newP0[j]=(1-p0)*P0[j]+(1-p1)*P1[j]
          newP1[j+1]=0
          newQ0[j]=(1-p0)*Q0[j]+(1-p1)*Q1[j]
          newQ1[j+1]=(p0)*Q0[j]+(p1)*Q1[j]+(p0)*P0[j]+(p1)*P1[j]
          newC0[j]=(1-p0)*C0[j]+(1-p1)*C1[j]
          newC1[j+1]=(p0)*C0[j]+(p1)*C1[j]
        }
      }
      else if(i==dayof){ #=5
        for(j in seq(7)){
          newP0[j]=(1-p0)*P0[j]+(1-p1)*P1[j]
          newP1[j+1]=(p0)*P0[j]+(p1)*P1[j]
          newQ0[j]=(1-p0)*Q0[j]+(1-p1)*Q1[j]
          newQ1[j+1]=(1-g.post)*(p0*Q0[j]+p1*Q1[j])
          newC0[j]=(1-p0)*C0[j]+(1-p1)*C1[j]
          newC1[j+1]=(p0)*C0[j]+(p1)*C1[j]+g.post*(p0*Q0[j]+p1*Q1[j])
        }
      }
      else{
        for(j in seq(7)){
          newP0[j]=(1-p0)*P0[j]+(1-p1)*P1[j]
          newP1[j+1]=(p0)*P0[j]+(p1)*P1[j]
          newQ0[j]=(1-p0)*Q0[j]+(1-p1)*Q1[j]
          newQ1[j+1]=(p0)*Q0[j]+(p1)*Q1[j]
          newC0[j]=(1-p0)*C0[j]+(1-p1)*C1[j]
          newC1[j+1]=(p0)*C0[j]+(p1)*C1[j]
        }
      }
      P0=newP0
      P1=newP1
      Q0=newQ0
      Q1=newQ1
      C0=newC0
      C1=newC1
    }
    
    efficacy.by.pill.1=c(0,0,.76,.76,.96,.96,.96,.99)
    efficacy.by.pill.2=c(0,0,.96,.96,.96,.96,.96,.99)
    
    c(efficacy=sum(efficacy.by.pill.1*(P0+P1+Q0+Q1)+efficacy.by.pill.2*(C0+C1)),pillrate=p, sensitivity=NA,specificity=NA,adherence=p)
  })
} #not called
#output: efficacy
#calls: expit
#input: adherence.params

calc.risk=function(riskparams,efficacy,sexparams){
  with(as.list(c(sexparams,riskparams)),{
    s=exp(sexfrequency)
    a=analfrequency
    ca=expit(condomanal)
    cv=expit(condomvaginal)
    risk.v=risk.vaginal*(1-condom.efficacy*cv)*(1-efficacy)
    risk.a=risk.anal*(1-condom.efficacy*ca)*(1-efficacy)
    c(risk.v=risk.v,risk.a=risk.a,risk=risk.a*a+risk.v*(1-v))
  })  
} #not called
#output: risk.v,risk.a,risk=risk.a*a+risk.v*(1-v)
#input: sexparams,riskparams,efficacy