#    TREATMENT CHOICE COX MODELS
#    AUTHOR - JAMES F. TROENDLE, NOVEMBER, 2017
#    PACKAGE TCCOX



#---------------------------------------------
# Define functions of Package tccox
#---------------------------------------------


treatinit=function(dataset,nmaxint) {
#---------------------------------------------
#------------------------------------------------------------
# TREATINIT
# R VERSION
#Count number of subjects in input dataset
#Count number of subjects who initiate treatment and those who then have an event
#Create sequential ID
#Define start_times[], stop_times[], and followup[] vectors
#Define tti[], tts[] vectors
#---------------------------------------------

obsnum=length(dataset$id)
idnum=0
idprev=0
i=1
while (i <= obsnum) {
  idcurr=dataset$id[i]
  if (idcurr != idprev) {idnum=idnum+1}
  dataset$id2[i]=idnum
  idprev=idcurr
  i=i+1
}
dataset$id=dataset$id2
drop=c("id2")
dataset=dataset[,!(names(dataset) %in% drop)]


#---------------------------------------------
#cat('Number of observations=',obsnum,'\n')
#cat('Dataset with sequential IDs=','\n')
#print(dataset[1:50,])
#cat('OBSNUM=',obsnum,'\n')
#---------------------------------------------
#---------------------------------------------

ntreatedevents=0
nperson=dataset$id[obsnum]
followup=rep(0,times=nperson)
treated=rep(0,times=nperson)
lastobs=rep(0,times=obsnum)
tti=rep(0,times=obsnum)
tts=rep(0,times=obsnum)
start_times=rep(0,times=nperson)
stop_times=rep(0,times=nperson)
idprev=0
i=1
while (i <=obsnum) {
  idcurr=dataset$id[i]
  followup[idcurr]=dataset$stop[i]
  tti[i]=start_times[idcurr]
  tts[i]=stop_times[idcurr]
  if (idcurr != idprev) {
    iobs=1
    itreated=0
    istarted=0
    istopped=0
    lastobs[i]=1
  } else {
    iobs=iobs+1
    lastobs[i]=1
    if (i != 1) {lastobs[i-1]=0}
  }
  if (dataset$treatment[i]==1 & itreated==0 & istarted==0) {
    start_times[idcurr]=dataset$start[i]
    j=1
    while (j <= iobs) {
      tti[i-j+1]=start_times[idcurr]
      j=j+1
    }
    treated[idcurr]=1
    itreated=1
    istarted=1
  }
  if (dataset$treatment[i]==0 & itreated==1 & istopped==0) {
    stop_times[idcurr]=dataset$start[i]
    j=1
    while (j <= iobs) {
      tts[i-j+1]=stop_times[idcurr]
      j=j+1
    }
    itreated=0
    istopped=1
  }
  if (dataset$status[i]==1 & treated[idcurr]==1) {ntreatedevents=ntreatedevents+1}

  idprev=idcurr
  i=i+1
}
ntreated=sum(treated)
numevents=sum(dataset$status)
nuntreatedevents=numevents-ntreatedevents
medianfollowup=stats::median(followup)
dataset=transform(dataset,lastobs=lastobs,tti=tti,tts=tts)
dataset=dataset[1:obsnum,]

treatmissing=0
numtreattrial=0
treatmissing=sum(treated==0)
numtreattrial=sum(treated==1)
numtreat=nperson-treatmissing
maxobs=nmaxint*nperson


#------------------------------------------------------------
#cat('Dataset after R VERSION of treatinit=','\n')
#print(dataset[1:50,])

#cat('OBSNUM=',obsnum,'\n')
#cat('NUM SUBJECTS=',nperson,'\n')
#cat('MAXOBS=',maxobs,'\n')
#cat('NUM EVENTS=',numevents,'\n')
#cat('NUM TREATED SUBJECTS=',ntreated,'\n')
#cat('NUM TREATED EVENTS=',ntreatedevents,'\n')

#------------------------------------------------------------
#------------------------------------------------------------
# END TREATINIT FUNCTION
#------------------------------------------------------------
list(dataset,obsnum,nperson,maxobs,numevents,start_times,stop_times,followup,tti,tts,medianfollowup)
}








itcadd=function(dataset,nmaxint,interval_width,min_exp_events,nitc_fixed,n_start_fixed,n_stop_fixed,interval_stop_beginning) {
#-----------------------------------------------------
# ITCADD FUNCTION
#-----------------------------------------------------
#-----------------------------------------------------
#DETERMINE NUMBER AND LOCATION OF ITC STARTING INTERVALS
#-----------------------------------------------------

if (nitc_fixed==0) {
  itc_start_endpoint=list()
  no_start=0
  previous_interval=0.0
  lastobs=subset(dataset,lastobs==1)
  event_rate=sum(lastobs$status)/length(lastobs$status)
#  cat('event_rate=',event_rate,'\n')

  nitc_start=0
  j=1
  while (j <= nmaxint) {
# Adjust itc_start_endpoint based on available data
#
    candidate=j*interval_width
    num_cell_s=sum((lastobs$tti > candidate-interval_width & lastobs$tti <= candidate & lastobs$stop > candidate))
    num_future_s=sum((lastobs$tti > candidate-interval_width & lastobs$stop > candidate))
    num_cell_u=sum((lastobs$stop > candidate & (lastobs$tti==0 | lastobs$tti > candidate)))
    num_events=sum((dataset$stop > candidate & dataset$status==1))
    if (no_start==0) {
      if (min(num_cell_s,num_cell_u)*event_rate >= min_exp_events) {
        itc_start_endpoint=c(itc_start_endpoint,j*interval_width)
        previous_interval=j*interval_width
        nitc_start=nitc_start+1
      }
    }
    if (no_start==0 & num_future_s*event_rate < min_exp_events) {
      no_start=1
      j=floor(previous_interval/interval_width)
    }
    if ((num_cell_s+num_cell_u)*event_rate < 2.*min_exp_events) {break}
    if (no_start==1) {break}
    j=j+1
  }
#-----------------------------------------------------
#DETERMINE NUMBER AND LOCATION OF ITC STOPPING INTERVALS
#-----------------------------------------------------
  itc_stop_endpoint=list()
  previous_interval=0
  no_stop=0

  nitc_stop=0
  j=2
  while (j <= nmaxint) {
# Adjust itc_stop_endpoint based on available data
#
    candidate=j*interval_width
    num_cell_stop=sum((lastobs$tti != 0 & lastobs$tti <= candidate-interval_width & lastobs$tts > candidate-interval_width & lastobs$tts <= candidate & lastobs$stop > candidate))
    num_cell_unstop=sum((lastobs$tti != 0 & lastobs$tti <= candidate-interval_width & (lastobs$tts != 0 | lastobs$tts > candidate) & lastobs$stop > candidate))
    num_future_stop=sum((lastobs$tts > candidate-interval_width & lastobs$stop > candidate))
    num_events=sum((dataset$stop > candidate & dataset$status==1))
    if (no_stop==0) {
      if (min(num_cell_stop,num_cell_unstop)*event_rate >= min_exp_events) {
        itc_stop_endpoint=c(itc_stop_endpoint,j*interval_width)
        previous_interval=j*interval_width
        nitc_stop=nitc_stop+1
      }
    }
    if (no_stop==0 & num_future_stop*event_rate < min_exp_events) {
      no_stop=1
      j=floor(previous_interval/interval_width)
    }
    if (no_stop==1) {break}
    j=j+1
  }
} else {

  itc_start_endpoint=list()
  nitc_start=n_start_fixed
  j=1
  while (j <= nitc_start) {
    itc_start_endpoint=c(itc_start_endpoint,j*interval_width)
    j=j+1
  }
  itc_stop_endpoint=list()
  nitc_stop=n_stop_fixed
  j=1
  while (j <= nitc_stop) {
    itc_stop_endpoint=c(itc_stop_endpoint,interval_stop_beginning+(j-1)*interval_width)
    j=j+1
  }
}

list(nitc_start,itc_start_endpoint,nitc_stop,itc_stop_endpoint)
}
# END OF ITCADD FUNCTION
#-------------------------------------------------










addtc=function(dataset,ncov,maxfollow,start_times,stop_times,min_future_events,numevents,nperson,nmaxint,maxobs,interval_width,
               nitc_start,itc_start_endpoint,nitc_stop,itc_stop_endpoint,tti,tts,followup) {
#-----------------------------------------------------
# ADDTC FUNCTION
#-----------------------------------------------------

order_start_times=sort(start_times[which(start_times!=0)])
start_times=unique(start_times)
start_times=sort(start_times[which(start_times!=0)])
nstart_times=length(start_times)
#cat('Number of unique start times=',nstart_times,'\n')
nstarts=length(order_start_times)-ceiling(min_future_events*nperson/numevents)
if (nstarts < 0) {nstarts=0}
if (nstarts > 0) {
  last_start_time=order_start_times[nstarts]
#  cat('Last start time utilized=',last_start_time,'\n')
  index=c(1:length(start_times))
  nstart_fit=max(index[start_times <= last_start_time])
} else {
  last_start_time=0
#  cat('Last start time utilized=',last_start_time,'\n')
  nstart_fit=0
}
if (nstart_fit < 0) {nstart_fit=0}
#cat('Number of unique start times used to guarantee expected events from future starters>=',min_future_events,'is ',nstart_fit,'\n')
if (nstart_fit > 0) {start_times=start_times[1:nstart_fit]}


order_stop_times=sort(stop_times[which(stop_times!=0)])
stop_times=unique(stop_times)
stop_times=sort(stop_times[which(stop_times!=0)])
nstop_times=length(stop_times)
#cat('Number of unique stop times=',nstop_times,'\n')
nstops=length(order_stop_times)-ceiling(min_future_events*nperson/numevents)
if (nstops < 0) {nstops=0}
if (nstops > 0) {
  last_stop_time=order_stop_times[nstops]
#  cat('Last stop time utilized=',last_stop_time,'\n')
  index=c(1:length(stop_times))
  nstop_fit=max(index[stop_times <= last_stop_time])
} else {
  last_stop_time=0
#  cat('Last stop time utilized=',last_stop_time,'\n')
  nstop_fit=0
}
if (nstop_fit < 0) {nstop_fit=0}
#cat('Number of unique stop times used to guarantee expected events from future stopers>=',min_future_events,'is ',nstop_fit,'\n')
if (nstop_fit > 0) {stop_times=stop_times[1:nstop_fit]}


n=length(dataset$id)
id2=dataset$id
id=c(dataset$id,rep(0,times=maxobs-n))
start2=dataset$start
start=c(dataset$start,rep(0,times=maxobs-n))
stop2=dataset$stop
stop=c(dataset$stop,rep(0,times=maxobs-n))
status2=dataset$status
status=c(dataset$status,rep(0,times=maxobs-n))
tstartp=rep(0,times=maxobs*nmaxint)
dim(tstartp)=c(maxobs,nmaxint)
tstopp=rep(0,times=maxobs*nmaxint)
dim(tstopp)=c(maxobs,nmaxint)
tti2=tti
tti=c(dataset$tti,rep(0,times=maxobs-n))
tts2=tts
tts=c(dataset$tts,rep(0,times=maxobs-n))
cov=rep(0,times=ncov*maxobs)
dim(cov)=c(maxobs,ncov)
cov2=rep(0,times=ncov*n)
dim(cov2)=c(n,ncov)
j=1
while (j <= ncov) {
  cov2[,j]=dataset[,j]
  cov[,j]=c(dataset[,j],rep(0,times=maxobs-n))
  j=j+1
}
tstart=rep(0,times=maxobs)
tstart1=rep(0,times=maxobs)
tstart0=rep(0,times=maxobs)
tstop=rep(0,times=maxobs)
tstop1=rep(0,times=maxobs)
tstop0=rep(0,times=maxobs)
itc_start_endpoint=c(itc_start_endpoint,rep(0,times=nmaxint-nitc_start))
itc_stop_endpoint=c(itc_stop_endpoint,rep(0,times=nmaxint-nitc_stop))

#cat('Before R Version of ADDWHI2.F','\n')
#cat('maxobs=',maxobs,'\n')
#cat('nmaxint=',nmaxint,'\n')
#cat('nitc_start=',nitc_start,'\n')
#cat('nitc_stop=',nitc_stop,'\n')
#cat('n=',n,'\n')
#cat('hormone=',hormone,'\n')
#cat('hormone1=','\n')
#cat('length(start_times)=',length(start_times),'\n')
#cat('length(stop_times)=',length(stop_times),'\n')
#cat('tti=','\n')
#print(tti[1:50])
#cat('length(tti)=',length(tti),'\n')
#cat('tts=','\n')
#print(tts[1:50])
#cat('length(tts)=',length(tts),'\n')
#cat('followup=','\n')
#print(followup[1:50])
#cat('length(followup)=',length(followup),'\n')


#-----------------------------------------
# Reduce list of starttime intervals to those with at least one start
#-----------------------------------------
startint=rep(0,times=nmaxint)
nstartint=0
startcurr=0
istartcurr=0
j=1
while (j <= nmaxint) {
  startcurr=startcurr+interval_width
  startcurr=round(startcurr,digits=2)
  if (istartcurr < nstart_fit) {
    if (start_times[istartcurr+1] < startcurr+.000001) {
      nstartint=nstartint+1
      istartcurr=istartcurr+1
      startint[nstartint]=startcurr
      if (istartcurr < nstart_fit) {
        while (istartcurr < nstart_fit & start_times[istartcurr+1] < startcurr+.000001) {istartcurr=istartcurr+1}
      }
    }
  }
  j=j+1
}

nriskstart=rep(0,times=nstartint)
n1start=rep(0,times=nstartint)
fstart=rep(0,times=nstartint)
npind=rep(0,times=nstartint)

newnitc_start=0
j=1
while (j <= nstartint) {
#  cat('j=',j,'\n')
  jj=newnitc_start+1
  while (jj <= nitc_start) {
#    cat('jj=',jj,'\n')
#    cat('startint[j]=',startint[j],'\n')
#    cat('itc_start_endpoint[jj]=','\n')
#    print(itc_start_endpoint[jj])
    if (abs(startint[j]-itc_start_endpoint[jj]) < .000001) {
      npind[j]=1
      newnitc_start=newnitc_start+1
      itc_start_endpoint[newnitc_start]=itc_start_endpoint[jj]
    }
    jj=jj+1
  }
  j=j+1
}
nitc_start=newnitc_start

idprev=0
idobs=0
i=1
while (i <= n) {
  idcurr=id[i]
  idobs=idobs+1
  if (idcurr != idprev) {idobs=1}
  if ((idobs==1 & cov2[i,1] <= .000001) | (idobs==2 & cov2[i,1] <= .000001)) {
    j=1
    while (j <= nstartint) {
      if (stop[i]-startint[j]+interval_width > .000001 & followup[idcurr]-startint[j] > .000001) {
        nriskstart[j]=nriskstart[j]+1
        if (tti[i] != 0 & tti[i]-startint[j] < .000001) {n1start[j]=n1start[j]+1}
        j=j+1
      } else {
        j=nstartint+1
      }
    }
  }
  idprev=idcurr
  i=i+1
}

j=1
while (j <= nstartint) {
  fstart[j]=n1start[j]/nriskstart[j]
#  cat('j=',j,'\n')
#  cat('nriskstart=',nriskstart[j],'\n')
#  cat('fstart=',fstart[j],'\n')
  j=j+1
}

#cat('Num of intervals with >= 1 start used in fitting=',nstartint,'\n')
#cat('starttime interval endpoints=',startint[1:nstartint],'\n')
#cat('----------------------------------------------------','\n')
#cat('nitc_start=',nitc_start,'\n')
#cat('itc_start_endpoint=','\n')
#print(itc_start_endpoint[1:nitc_start])
#cat('----------------------------------------------------','\n')

#-----------------------------------------
# Section on Stop times
#-----------------------------------------
stopint=rep(0,times=nstop_fit)
j=1
while (j <= nstop_fit) {
  stopint[j]=(floor(stop_times[j]/interval_width)+1)*interval_width
  j=j+1
}

nstopint=0
current=0
j=1
while (j <= nstop_fit) {
  if (abs(stopint[j]-current) > .000001) {
    nstopint=nstopint+1
    stopint[nstopint]=stopint[j]
    current=stopint[nstopint]
  }
  j=j+1
}
#cat('Num of intervals with >= 1 stop=',nstopint,'\n')
#cat('stoptime interval endpoints=',stopint[1:nstopint],'\n')

nriskstop=rep(0,times=nstopint)
n1stop=rep(0,times=nstopint)
fstop=rep(0,times=nstopint)
npindstop=rep(0,times=nstopint)

idobs=0
idprev=0
i=1
while (i <= n) {
  idcurr=id[i]
  idobs=idobs+1
  if (idcurr != idprev) {idobs=1}
  if ((idobs==1 & cov2[i,1] >= .999999) | (idobs==2 & cov2[i,1] >= .999999)) {
    j=1
    while (j <= nstopint) {
      if (stop[i]-stopint[j]+interval_width > .000001 & followup[idcurr]-stopint[j] > .000001 & tti[i] != 0 & tti[i]-stopint[j]+interval_width < .000001) {
        nriskstop[j]=nriskstop[j]+1
        if (tts[i] != 0 & tts[i]-stopint[j] < .000001) {n1stop[j]=n1stop[j]+1}
        j=j+1
      } else {
        j=j+1
      }
    }
  }
  idprev=idcurr
  i=i+1
}

#-----------------------------------------
# Reduce list of stoptime intervals to those with at least one at risk stop and at risk no-stop
#-----------------------------------------
nstopintnew=0
istopcurr=0
j=1
while (j <= nstopint) {
  if (n1stop[j] > 0 & nriskstop[j]-n1stop[j] > 0) {
    nstopintnew=nstopintnew+1
    stopint[nstopintnew]=stopint[j]
    n1stop[nstopintnew]=n1stop[j]
    nriskstop[nstopintnew]=nriskstop[j]
    fstop[nstopintnew]=fstop[j]
  }
  j=j+1
}
nstopint=nstopintnew


newnitc_stop=0
j=1
while (j <= nstopint) {
  jj=newnitc_stop+1
  while (jj <= nitc_stop) {
    if (abs(stopint[j]-itc_stop_endpoint[jj]) < .000001) {
      npindstop[j]=1
      newnitc_stop=newnitc_stop+1
      itc_stop_endpoint[newnitc_stop]=itc_stop_endpoint[jj]
    }
    jj=jj+1
  }
  j=j+1
}
nitc_stop=newnitc_stop

j=1
while (j <= nstopint) {
  fstop[j]=n1stop[j]/nriskstop[j]
#  cat('j=',j,'\n')
#  cat('nriskstop=',nriskstop[j],'\n')
#  cat('fstop=',fstop[j],'\n')
  j=j+1
}

#cat('Num of intervals with >= 1 stop used in fitting=',nstopint,'\n')
#cat('stoptime interval endpoints used=',stopint[1:nstopint],'\n')
#cat('----------------------------------------------------','\n')
#cat('nitc_stop=',nitc_stop,'\n')
#cat('itc_stop_endpoint=','\n')
#print(itc_stop_endpoint[1:nitc_stop])
#cat('----------------------------------------------------','\n')

#-----------------------------------------
# Create combined list of interval endpoints
#-----------------------------------------
combint=rep(0,times=2*nmaxint)
icombstart=rep(0,times=2*nmaxint)
icombstop=rep(0,times=2*nmaxint)
npcombindstart=rep(0,times=2*nmaxint)
npcombindstop=rep(0,times=2*nmaxint)
fcombstart=rep(0,times=2*nmaxint)
fcombstop=rep(0,times=2*nmaxint)
nfit=0
jj=1
j=1
while (j <= nstartint) {
  while (jj <= nstopint & stopint[jj]-startint[j] < -.000001) {
    nfit=nfit+1
    combint[nfit]=stopint[jj]
    icombstop[nfit]=1
    npcombindstop[nfit]=0
    if (npindstop[jj]==1) {npcombindstop[nfit]=1}
    fcombstop[nfit]=fstop[jj]
    npcombindstart[nfit]=0
    icombstart[nfit]=0
    fcombstart[nfit]=0
    jj=jj=1
  }
  if (jj <= nstopint) {
    nfit=nfit+1
    combint[nfit]=startint[j]
    icombstart[nfit]=1
    npcombindstart[nfit]=0
    if (npind[j]==1) {npcombindstart[nfit]=1}
    fcombstart[nfit]=fstart[j]
    if (abs(stopint[jj]-startint[j]) < .000001) {
      icombstop[nfit]=1
      npcombindstop[nfit]=0
      if (npindstop[jj]==1) {npcombindstop[nfit]=1}
      fcombstop[nfit]=fstop[jj]
      jj=jj+1
    } else {
      icombstop[nfit]=0
      fcombstop[nfit]=0
      npcombindstop[nfit]=0
    }
  } else {
    nfit=nfit+1
    combint[nfit]=startint[j]
    icombstart[nfit]=1
    npcombindstart[nfit]=0
    if (npind[j]==1) {npcombindstart[nfit]=1}
    fcombstart[nfit]=fstart[j]
    icombstop[nfit]=0
    fcombstop[nfit]=0
    npcombindstop[nfit]=0
  }
  j=j+1
}
j=jj
while (j <= nstopint) {
  nfit=nfit+1
  combint[nfit]=stopint[j]
  icombstop[nfit]=1
  npcombindstop[nfit]=0
  if (npindstop[j]==1) {npcombindstop[nfit]=1}
  fcombstop[nfit]=fstop[j]
  npcombindstart[nfit]=0
  icombstart[nfit]=0
  fcombstart[nfit]=0

  j=j+1
}

#cat('------------------------------------------','\n')
#cat('Number of combined TC intervals used=',nfit,'\n')
#cat('------------------------------------------','\n')
#j=1
#while (j <= nfit) {
#  cat('------------------------------------------','\n')
#  cat('j=',j,'\n')
#  cat('combint=',combint[j],'\n')
#  cat('icombstart=',icombstart[j],'\n')
#  cat('npcombindstart=',npcombindstart[j],'\n')
#  cat('fcombstart=',fcombstart[j],'\n')
#  cat('icombstop=',icombstop[j],'\n')
#  cat('npcombindstop=',npcombindstop[j],'\n')
#  cat('fcombstop=',fcombstop[j],'\n')
#  j=j+1
#}
#cat('------------------------------------------','\n')
#cat('Create Intervals','\n')
#cat('n=',n,'\n')

#----------------------------------------
# Create Intervals in Data Set
#----------------------------------------
newn=0
idprev=0
i=1
while (i <= n) {
  idcurr=id2[i]

#  cat('i=',i,'\n')
#  cat('----------------------------','\n')
#  cat('----------------------------','\n')
#  cat('idprev=',idprev,'\n')
#  cat('idcurr=',idcurr,'\n')
#  cat('start=',start2[i],'\n')
#  cat('stop=',stop2[i],'\n')
#  cat('status=',status2[i],'\n')
#  cat('treatment=',cov2[i,1],'\n')

  if (idcurr != idprev) {
    idnum=1
    istart=0
    istop=0
    icontinue=0
    j=1
    jnp=0
    jnpstop=0
    lastj=0
    lastjstop=0
    tstartcurr=0
    tstart1curr=0
    tstart0curr=0
    tstartstar=0
    tstopcurr=0
    tstop1curr=0
    tstop0curr=0
    tstopstar=0
    stopprev=0
  } else {
    idnum=idnum+1
  }
  if (nfit==0) {icontinue=1}

  tempstop=stop2[i]-1
  while (abs(tempstop-stop2[i]) >= .000001) {
#    cat('----------------------------','\n')
#    cat('first enter while loop','\n')
#    cat('j=',j,'\n')
#    cat('icontinue=',icontinue,'\n')
#    cat('istart=',istart,'\n')
#    cat('stopprev=',stopprev,'\n')
#    cat('tempstop=',tempstop,'\n')
#    cat('----------------------------','\n')
    tempstop=0
    tempstart=0
    if (icontinue==0) {
      tempstop=min(stop2[i],combint[j])
      if (j==1) {
        tempstart=max(start2[i],0)
      } else {
        tempstart=max(start2[i],combint[j-1])
      }
    }
    if (icontinue==1 | j > nfit) {
      tempstop=min(stop2[i],maxfollow)
      tempstart=max(start2[i],stopprev)
    }
#    cat('just before test','\n')
#    cat('j=',j,'\n')
#    cat('icontinue=',icontinue,'\n')
#    cat('stopprev=',stopprev,'\n')
#    cat('tempstart=',tempstart,'\n')
#    cat('tempstop=',tempstop,'\n')
#    cat('----------------------------','\n')
    stopprev=tempstop
    if (tempstart+.000001 < tempstop) {
      newn=newn+1
      id[newn]=id2[i]
      start[newn]=tempstart
      stop[newn]=tempstop
      status[newn]=status2[i]
      if (abs(tempstop-stop2[i]) >= .000001) {status[newn]=0}
      tti[newn]=tti2[i]
      tts[newn]=tts2[i]

      k=1
      while (k <= ncov) {
        cov[newn,k]=cov2[i,k]
        k=k+1
      }

#      cat('Add to data_adj','\n')
#      cat('newn=',newn,'\n')
#      cat('id=',id[newn],'\n')
#      cat('start=',start[newn],'\n')
#      cat('stop=',stop[newn],'\n')
#      cat('status=',status[newn],'\n')
#      cat('treatment=',cov[newn,1],'\n')
#      cat('----------------------------','\n')

      jj=1
      while (jj <= jnp) {
        tstartp[newn,jj]=tstartp[newn-1,jj]
        jj=jj+1
      }
      if (tstartstar != 0) {
        tstartp[newn,jnp]=0
        tstartstar=0
      }
      jj=1
      while (jj <= jnpstop) {
        tstopp[newn,jj]=tstopp[newn-1,jj]
        jj=jj+1
      }
      if (tstopstar != 0) {
        tstopp[newn,jnpstop]=0
        tstopstar=0
      }
    }

    if (j > 1 & icontinue==0 & istart==0) {
    if (icombstart[j-1]==1) {
      if ((tti2[i]==0 | tti2[i] > combint[j-1]-interval_width+.000001) & abs(tempstart-combint[j-1]) < .000001 & cov2[i,1]==0) {
        if (tti2[i] != 0 & tti2[i] <= combint[j-1]+.000001 & cov2[i,1]==0) {
          tstartcurr=tstartcurr+1-fcombstart[j-1]
          tstart1curr=tstart1curr+(1-fcombstart[j-1])*combint[j-1]
          if (npcombindstart[j-1]==0) {
            tstart0curr=tstart0curr+1-fcombstart[j-1]
          } else {
            if (tempstart+.000001 < tempstop) {
              jnp=jnp+1
              tstartp[newn,jnp]=1-fcombstart[j-1]
            } else {
              tstartstar=1-fcombstart[j-1]
            }
          }
        } else {
          tstartcurr=tstartcurr-fcombstart[j-1]
          tstart1curr=tstart1curr-(fcombstart[j-1])*combint[j-1]
          if (npcombindstart[j-1]==0) {
            tstart0curr=tstart0curr-fcombstart[j-1]
          } else {
            if (tempstart+.000001 < tempstop) {
              jnp=jnp+1
              tstartp[newn,jnp]=-1*fcombstart[j-1]
            } else {
              tstartstar=-1*fcombstart[j-1]
            }
          }
        }
      }
    }
    }

    if (lastj >= 1 & icontinue==0 & istart==1) {
    if (icombstart[lastj]==1) {
      if ((tti2[i]==0 | tti2[i] > combint[lastj]-interval_width+.000001) & abs(tempstart-combint[lastj]) < .000001 & cov2[i,1]==1) {
        if (tti2[i] !=0 & tti2[i] <= combint[lastj]+.000001 & cov2[i,1]==1) {
          tstartcurr=tstartcurr+1-fcombstart[lastj]
          tstart1curr=tstart1curr+(1-fcombstart[lastj])*combint[lastj]
          if (npcombindstart[lastj]==0) {
            tstart0curr=tstart0curr+1-fcombstart[lastj]
          } else {
            if (tempstart+.000001 < tempstop) {
              jnp=jnp+1
              tstartp[newn,jnp]=1-fcombstart[lastj]
            } else {
              tstartstar=1-fcombstart[lastj]
            }
          }
        } else {
          tstartcurr=tstartcurr-fcombstart[lastj]
          tstart1curr=tstart1curr-(fcombstart[lastj])*combint[lastj]
          if (npcombindstart[lastj]==0) {
            tstart0curr=tstart0curr-fcombstart[lastj]
          } else {
            if (tempstart+.000001 < tempstop) {
              jnp=jnp+1
              tstartp[newn,jnp]=-1*fcombstart[lastj]
            } else {
              tstartstar=-1*fcombstart[lastj]
            }
          }
        }
      }
    }
    }

    if (lastj >= 1 & icontinue==1 & istart==0) {
    if (icombstart[lastj]==1) {
      if ((tti2[i]==0 | tti2[i] > combint[lastj]-interval_width+.000001) & abs(tempstart-combint[lastj]) < .000001 & cov2[i,1]==0) {
        if (tti2[i] !=0 & tti2[i] <= combint[lastj]+.000001 & cov2[i,1]==0) {
          tstartcurr=tstartcurr+1-fcombstart[lastj]
          tstart1curr=tstart1curr+(1-fcombstart[lastj])*combint[lastj]
          if (npcombindstart[lastj]==0) {
            tstart0curr=tstart0curr+1-fcombstart[lastj]
          } else {
            if (tempstart+.000001 < tempstop) {
              jnp=jnp+1
              tstartp[newn,jnp]=1-fcombstart[lastj]
            } else {
              tstartstar=1-fcombstart[lastj]
            }
          }
        } else {
          tstartcurr=tstartcurr-fcombstart[lastj]
          tstart1curr=tstart1curr-(fcombstart[lastj])*combint[lastj]
          if (npcombindstart[lastj]==0) {
            tstart0curr=tstart0curr-fcombstart[lastj]
          } else {
            if (tempstart+.000001 < tempstop) {
              jnp=jnp+1
              tstartp[newn,jnp]=-1*fcombstart[lastj]
            } else {
              tstartstar=-1*fcombstart[lastj]
            }
          }
        }
      }
    }
    }

    if (j > 1 & icontinue==0 & istop==0) {
    if (icombstop[j-1]==1) {
      if ((tts2[i]==0 | tts2[i] > combint[j-1]-interval_width+.000001) & abs(tempstart-combint[j-1]) < .000001 & cov2[i,1]==1) {
        if (tts2[i] != 0 & tts2[i] <= combint[j-1]+.000001 & cov2[i,1]==1) {
          tstopcurr=tstopcurr+1-fcombstop[j-1]
          tstop1curr=tstop1curr+(1-fcombstop[j-1])*combint[j-1]
          if (npcombindstop[j-1]==0) {
            tstop0curr=tstop0curr+1-fcombstop[j-1]
          } else {
            if (tempstart+.000001 < tempstop) {
              jnpstop=jnpstop+1
              tstopp[newn,jnpstop]=1-fcombstop[j-1]
            } else {
              tstopstar=1-fcombstop[j-1]
            }
          }
        } else {
          tstopcurr=tstopcurr-fcombstop[j-1]
          tstop1curr=tstop1curr-(fcombstop[j-1])*combint[j-1]
          if (npcombindstop[j-1]==0) {
            tstop0curr=tstop0curr-fcombstop[j-1]
          } else {
            if (tempstart+.000001 < tempstop) {
              jnpstop=jnpstop+1
              tstopp[newn,jnpstop]=-1*fcombstop[j-1]
            } else {
              tstopstar=-1*fcombstop[j-1]
            }
          }
        }
      }
    }
    }

    if (lastjstop >= 1 & icontinue==0 & istop==1) {
    if (icombstop[lastjstop]==1) {
      if ((tts2[i]==0 | tts2[i] > combint[lastjstop]-interval_width+.000001) & abs(tempstart-combint[lastjstop]) < .000001 & cov2[i,1]==0) {
        if (tts2[i] != 0 & tts2[i] <= combint[lastjstop]+.000001 & cov2[i,1]==0) {
          tstopcurr=tstopcurr+1-fcombstop[lastjstop]
          tstop1curr=tstop1curr+(1-fcombstop[lastjstop])*combint[lastjstop]
          if (npcombindstop[lastjstop]==0) {
            tstop0curr=tstop0curr+1-fcombstop[lastjstop]
          } else {
            if (tempstart+.000001 < tempstop) {
              jnpstop=jnpstop+1
              tstopp[newn,jnpstop]=1-fcombstop[lastjstop]
            } else {
              tstopstar=1-fcombstop[lastjstop]
            }
          }
        } else {
          tstopcurr=tstopcurr-fcombstop[lastjstop]
          tstop1curr=tstop1curr-(fcombstop[lastjstop])*combint[lastjstop]
          if (npcombindstop[lastjstop]==0) {
            tstop0curr=tstop0curr-fcombstop[lastjstop]
          } else {
            if (tempstart+.000001 < tempstop) {
              jnpstop=jnpstop+1
              tstopp[newn,jnpstop]=-1*fcombstop[lastjstop]
            } else {
              tstopstar=-1*fcombstop[lastjstop]
            }
          }
        }
      }
    }
    }

    if (lastjstop >= 1 & icontinue==1 & istop==0) {
    if (icombstop[lastjstop]==1) {
      if ((tts2[i]==0 | tts2[i] > combint[lastjstop]-interval_width+.000001) & abs(tempstart-combint[lastjstop]) < .000001 & cov2[i,1]==1) {
        if (tts2[i] != 0 & tts2[i] <= combint[lastjstop]+.000001 & cov2[i,1]==1) {
          tstopcurr=tstopcurr+1-fcombstop[lastjstop]
          tstop1curr=tstop1curr+(1-fcombstop[lastjstop])*combint[lastjstop]
          if (npcombindstop[lastjstop]==0) {
            tstop0curr=tstop0curr+1-fcombstop[lastjstop]
          } else {
            if (tempstart+.000001 < tempstop) {
              jnpstop=jnpstop+1
              tstopp[newn,jnpstop]=1-fcombstop[lastjstop]
            } else {
              tstopstar=1-fcombstop[lastjstop]
            }
          }
        } else {
          tstopcurr=tstopcurr-fcombstop[lastjstop]
          tstop1curr=tstop1curr-(fcombstop[lastjstop])*combint[lastjstop]
          if (npcombindstop[lastjstop]==0) {
            tstop0curr=tstop0curr-fcombstop[lastjstop]
          } else {
            if (tempstart+.000001 < tempstop) {
              jnpstop=jnpstop+1
              tstopp[newn,jnpstop]=-1*fcombstop[lastjstop]
            } else {
              tstopstar=-1*fcombstop[lastjstop]
            }
          }
        }
      }
    }
    }

    if (tempstart+.000001 < tempstop) {
      tstart[newn]=tstartcurr
      tstart0[newn]=tstart0curr
      tstart1[newn]=tstart1curr
      tstop[newn]=tstopcurr
      tstop0[newn]=tstop0curr
      tstop1[newn]=tstop1curr
    }
    if (istart==0 & icombstart[j]==1) {lastj=j}
    if (istop==0 & icombstop[j]==1) {lastjstop=j}

    if (abs(tempstop-stop2[i]) >= .000001) {
      if (istop==0 & tts2[i] != 0 & combint[j]-tempstop <= .000001 & tts2[i]-tempstop < .000001) {
        istop=1
        j=j+1
      } else {
        if (istart==0 & tti2[i] != 0 & combint[j]-tempstop <= .000001 & tti2[i]-tempstop < .000001) {
          istart=1
          j=j+1
        } else {
          if (j==nfit) {
            icontinue=1
          } else {
            if (istart==1 & istop==1) {
              icontinue=1
            } else {
              if ((istart==0 | istop==0) & (tts2[i]-tempstop > .000001 | tts2[i]==0 | tti2[i]-tempstop > .000001 | tti2[i]==0)) {
                j=j+1
              }
            }
          }
        }
      }
    }

  }

  idprev=idcurr
  i=i+1
}
n=newn
#cat('Maximum number of obs during interval addition=',n,'\n')

#-----------------------------------------
#Remove un-needed Intervals
#-----------------------------------------
newn=0
idprev=0
i=1
while (i <= n) {
  idcurr=id[i]
  if (idcurr != idprev) {ichange=1}
  if (ichange != 1) {
    if (tstart[i] != tstart[newn]) {ichange=1}
    if (tstop[i] != tstop[newn]) {ichange=1}
    if (cov[i,1] != cov[newn,1]) {ichange=1}
  }
  if (ichange==1) {
    newn=newn+1
    tstart[newn]=tstart[i]
    tstart0[newn]=tstart0[i]
    tstart1[newn]=tstart1[i]
    tstop[newn]=tstop[i]
    tstop0[newn]=tstop0[i]
    tstop1[newn]=tstop1[i]
    id[newn]=id[i]
    start[newn]=start[i]
    stop[newn]=stop[i]
    status[newn]=status[i]
    cov[newn,1]=cov[i,1]
    tti[newn]=tti[i]
    tts[newn]=tts[i]

    j=1
    while (j <= ncov) {
      cov[newn,j]=cov[i,j]
      j=j+1
    }

    tstartp[newn,]=tstartp[i,]
    tstopp[newn,]=tstopp[i,]
  } else {
    stop[newn]=stop[i]
    status[newn]=status[i]
  }
  idprev=idcurr
  i=i+1
}
n=newn
#cat('Final number of obs after interval adding=',n,'\n')

dataset=data.frame(cov=cov,id=id,start=start,stop=stop,status=status,
                   tstart=tstart,tstart1=tstart1,tstart0=tstart0,tstop=tstop,tstop1=tstop1,tstop0=tstop0)
dataset=dataset[1:n,]

#cat('New OBSNUM=',n,'\n')
#cat('New Dataset=','\n')
#print(dataset[1:50,])
#cat('Number of ITC starting intervals=',nitc_start,'\n')
#cat('Number of ITC stopping intervals=',nitc_stop,'\n')

tstartp=tstartp[1:n,1:nitc_start]
tstopp=tstopp[1:n,1:nitc_stop]

#cat('------------------------------------------------------','\n')
#cat('------------------------------------------------------','\n')
#cat('tstartp=','\n')
#print(tstartp)
#cat('------------------------------------------------------','\n')
#cat('tstopp=','\n')
#print(tstopp)
#cat('------------------------------------------------------','\n')

list(dataset,tstartp,tstopp,nstartint,startint,nstopint,stopint,nitc_start,itc_start_endpoint,nitc_stop,itc_stop_endpoint)
}
# END OF ADDTC FUNCTION
#-------------------------------------------------











itcfitter=function(dataset,ncov,cov_names,nitc_start,itc_start_endpoint,nitc_stop,itc_stop_endpoint,tstartp,tstopp) {
#------------------------------------------------------
# FIT INTERVAL TREATMENT CHOICE MODEL (ITC)
#------------------------------------------------------
#cat('----------------------------------------------------------','\n')
#cat('----------------------------------------------------------','\n')
#cat('START INTERVAL TREATMENT CHOICE MODEL SECTION','\n')
#cat(' USES ANALYSIS WITH START AND STOP INTERVALS','\n')
#cat('----------------------------------------------------------','\n')
#cat('----------------------------------------------------------','\n')
data3=dataset
if (nitc_start > 0 & nitc_stop > 0) {
  data3=transform(dataset,treatstartp=tstartp,treatstopp=tstopp)
} else {
  if (nitc_start > 0) {
    data3=transform(dataset,treatstart=tstartp)
  }
  if (nitc_stop > 0) {
    data3=transform(dataset,treatstopp=tstopp)
  }
}
#cat('data3=','\n')
#print(data3[1:50,])

num_terms=nitc_start+nitc_stop
#cat('num_terms=',num_terms,'\n')
#cat('ncov=',ncov,'\n')

if (nitc_start+nitc_stop > 0) {
  names=paste("data3[,",c(1:ncov,(ncov+11):(ncov+11+num_terms-1)),sep="")
  names=paste(names,"]",sep="")
  (formula=stats::as.formula(paste("Surv(start,stop,status)~",paste(names,collapse="+"),paste("+cluster(id)"))))
  cov_names3=c(cov_names,names(data3)[(ncov+11):(ncov+11+num_terms-1)])
} else {
  names=paste("data3[,",1:ncov,sep="")
  names=paste(names,"]",sep="")
  (formula=stats::as.formula(paste("Surv(start,stop,status)~",paste(names,collapse="+"),paste("+cluster(id)"))))
  cov_names3=cov_names
}

fit_itc=survival::coxph(formula,data=data3)

#cat('ITC Cox Model=','\n')
#print(fit_itc)
#cat('cov_names=','\n')
#print(cov_names3)
#cat('----------------------------------------------------','\n')
#cat('nitc_start=',nitc_start,'\n')
#cat('itc_start_endpoint=','\n')
#print(itc_start_endpoint)
#cat('----------------------------------------------------','\n')
#cat('nitc_stop=',nitc_stop,'\n')
#cat('itc_stop_endpoint=','\n')
#print(itc_stop_endpoint)
#cat('----------------------------------------------------','\n')

#cat('----------------------------------------------------------','\n')
#cat('ITC Cox Model Treatment Effect log(HR)=',curr,'\n')
#cat('----------------------------------------------------------','\n')
treatcoef_tc3=stats::coef(fit_itc)[1]
treatse_tc3=sqrt(stats::vcov(fit_itc)[1,1])
#-----------------------------------------------------------------------
#-----------------------------------------------------------------------
list(treatcoef_tc3,treatse_tc3,fit_itc,cov_names3)
}
# END OF ITCFITTER FUNCTION
#-----------------------------------------------------------------------








htcfitter=function(dataset,ncov,cov_names,nitc_start,itc_start_endpoint,nitc_stop,itc_stop_endpoint,tstartp,tstopp) {
#------------------------------------------------------
# FIT HYBRID TREATMENT CHOICE MODEL (HTC)
#------------------------------------------------------
#cat('----------------------------------------------------------','\n')
#cat('----------------------------------------------------------','\n')
#cat('START HYBRID TREATMENT CHOICE MODEL SECTION','\n')
#cat('----------------------------------------------------------','\n')
#cat('----------------------------------------------------------','\n')
data2=dataset
if (nitc_start > 0 & nitc_stop > 0) {
  data2=transform(dataset,treatstartp=tstartp,treatstopp=tstopp)
} else {
  if (nitc_start > 0) {
    data2=transform(dataset,treatstart=tstartp)
  }
  if (nitc_stop > 0) {
    data2=transform(dataset,treatstopp=tstopp)
  }
}

num_terms=nitc_start+nitc_stop
#cat('num_terms=',num_terms,'\n')
#cat('ncov=',ncov,'\n')

if (nitc_start+nitc_stop > 0) {
  names=paste("data2[,",c(1:ncov,(ncov+11):(ncov+11+num_terms-1)),sep="")
  names=paste(names,"]",sep="")
  (formula=stats::as.formula(paste("Surv(start,stop,status)~",paste(names,collapse="+"),paste("+tstart0+tstop0+cluster(id)"))))
  cov_names2=c(cov_names,names(data2)[(ncov+11):(ncov+11+num_terms-1)],"tstart0","tstop0")
} else {
  names=paste("data2[,",1:ncov,sep="")
  names=paste(names,"]",sep="")
  (formula=stats::as.formula(paste("Surv(start,stop,status)~",paste(names,collapse="+"),paste("+cluster(id)"))))
  cov_names2=cov_names
}

fit_htc=survival::coxph(formula,data=data2)

#cat('HTC Cox Model =','\n')
#print(fit_htc)
#cat('cov_names=','\n')
#print(cov_names2)
#cat('----------------------------------------------------','\n')
#cat('nitc_start=',nitc_start,'\n')
#cat('itc_start_endpoint=','\n')
#print(itc_start_endpoint)
#cat('----------------------------------------------------','\n')
#cat('nitc_stop=',nitc_stop,'\n')
#cat('itc_stop_endpoint=','\n')
#print(itc_stop_endpoint)
#cat('----------------------------------------------------','\n')

#cat('----------------------------------------------------------','\n')
#cat('HTC Cox Model Treatment Effect log(HR)=',curr,'\n')
#cat('----------------------------------------------------------','\n')

treatcoef_tc2=stats::coef(fit_htc)[1]
treatse_tc2=sqrt(stats::vcov(fit_htc)[1,1])
#-----------------------------------------------------------------------
#-----------------------------------------------------------------------
list(treatcoef_tc2,treatse_tc2,fit_htc,cov_names2)
}
# END OF HTCFITTER FUNCTION
#-----------------------------------------------------------------------









ptcfitter=function(dataset,ncov,cov_names) {
#------------------------------------------------------
# FIT PARAMETRIC TC MODEL (PTC)
#------------------------------------------------------
#cat('----------------------------------------------------------','\n')
#cat('----------------------------------------------------------','\n')
#cat('START PARAMETRIC TC MODEL SECTION','\n')
#cat('----------------------------------------------------------','\n')
#cat('----------------------------------------------------------','\n')

names=paste("dataset[,",1:ncov,sep="")
names=paste(names,"]",sep="")
(formula=stats::as.formula(paste("Surv(start,stop,status)~",paste(names,collapse="+"),paste("+tstart+tstart1+tstop+tstop1+cluster(id)"))))
cov_names1=c(cov_names,"tstart","tstart1","tstop","tstop1")
fit_ptc=survival::coxph(formula,data=dataset)

#cat('PTC model is','\n')
#print(fit_ptc)
#cat('cov_names=','\n')
#print(cov_names1)

treatcoef_tc1=stats::coef(fit_ptc)[1]
treatse_tc1=sqrt(stats::vcov(fit_ptc)[1,1])
#-----------------------------------------------------------------------
#-----------------------------------------------------------------------
list(treatcoef_tc1,treatse_tc1,fit_ptc,cov_names1)
}
# END OF PTCFITTER FUNCTION
#-----------------------------------------------------------------------









coxfitter=function(dataset1,ncov,cov_names) {
#------------------------------------------
#Ordinary Cox model
#------------------------------------------

names=paste("dataset1[,",1:ncov,sep="")
names=paste(names,"]",sep="")
(formula=stats::as.formula(paste("Surv(start,stop,status)~",paste(names,collapse="+"),paste("+cluster(id)"))))
fit_cox=survival::coxph(formula,data=dataset1)


#cat('Ordinary Cox Model=','\n')
#print(fit_cox)
#cat('cov_names=','\n')
#print(cov_names)
#cat('--------------------------------------','\n')

treatcoef_cox1=stats::coef(fit_cox)[1]
treatse_cox1=sqrt(stats::vcov(fit_cox)[1,1])

list(treatcoef_cox1,treatse_cox1,fit_cox)
}
# END OF COXFITTER FUNCTION
#-----------------------------------------------------------------------







#------------------------------------------------------
# INTERVAL TREATMENT CHOICE MODEL (ITC)
#------------------------------------------------------
itc=function(dataset,ncov,cov_names,maxfollow,nmaxint,interval_width,min_exp_events,min_future_events,nitc_fixed,n_start_fixed,n_stop_fixed,interval_stop_beginning) {
  z=treatinit(dataset,nmaxint)
  dataset=z[[1]]
  obsnum=z[[2]]
  nperson=z[[3]]
  maxobs=z[[4]]
  numevents=z[[5]]
  start_times=z[[6]]
  stop_times=z[[7]]
  followup=z[[8]]
  tti=z[[9]]
  tts=z[[10]]
  medianfollowup=z[[11]]
  z=itcadd(dataset,nmaxint,interval_width,min_exp_events,nitc_fixed,n_start_fixed,n_stop_fixed,interval_stop_beginning)
  nitc_start=z[[1]]
  itc_start_endpoint=as.numeric(z[[2]])
  nitc_stop=z[[3]]
  itc_stop_endpoint=as.numeric(z[[4]])
  z=addtc(dataset,ncov,maxfollow,start_times,stop_times,min_future_events,numevents,nperson,nmaxint,maxobs,interval_width,
         nitc_start,itc_start_endpoint,nitc_stop,itc_stop_endpoint,tti,tts,followup)
  dataset=z[[1]]
  tstartp=z[[2]]
  tstopp=z[[3]]
  nstartint=z[[4]]
  startint=z[[5]][1:nstartint]
  nstopint=z[[6]]
  stopint=z[[7]][1:nstopint]
  nitc_start=z[[8]]
  if (nitc_start > 0) {
    itc_start_endpoint=z[[9]][1:nitc_start]
  } else {
    itc_start_endpoint=list()
  }
  nitc_stop=z[[10]]
  if (nitc_stop > 0) {
    itc_stop_endpoint=z[[11]][1:nitc_stop]
  } else {
    itc_stop_endpoint=list()
  }
  z=itcfitter(dataset,ncov,cov_names,nitc_start,itc_start_endpoint,nitc_stop,itc_stop_endpoint,tstartp,tstopp)
  treatcoef_tc3=z[[1]]
  treatse_tc3=z[[2]]
  fit_itc=z[[3]]
  cov_names3=z[[4]]
  list(fit_itc,nitc_start,itc_start_endpoint,nitc_stop,itc_stop_endpoint,nstartint,startint,nstopint,stopint,cov_names3,nperson,numevents,medianfollowup)
}
#------------------------------------------------------
# END: INTERVAL TREATMENT CHOICE MODEL FUNCTION
#------------------------------------------------------




#------------------------------------------------------
# HYBRID TREATMENT CHOICE MODEL (HTC)
#------------------------------------------------------
htc=function(dataset,ncov,cov_names,maxfollow,nmaxint,interval_width,min_exp_events,min_future_events,nitc_fixed,n_start_fixed,n_stop_fixed,interval_stop_beginning) {
  z=treatinit(dataset,nmaxint)
  dataset=z[[1]]
  obsnum=z[[2]]
  nperson=z[[3]]
  maxobs=z[[4]]
  numevents=z[[5]]
  start_times=z[[6]]
  stop_times=z[[7]]
  followup=z[[8]]
  tti=z[[9]]
  tts=z[[10]]
  medianfollowup=z[[11]]
  z=itcadd(dataset,nmaxint,interval_width,min_exp_events,nitc_fixed,n_start_fixed,n_stop_fixed,interval_stop_beginning)
  nitc_start=z[[1]]
  itc_start_endpoint=as.numeric(z[[2]])
  nitc_stop=z[[3]]
  itc_stop_endpoint=as.numeric(z[[4]])
  z=addtc(dataset,ncov,maxfollow,start_times,stop_times,min_future_events,numevents,nperson,nmaxint,maxobs,interval_width,
         nitc_start,itc_start_endpoint,nitc_stop,itc_stop_endpoint,tti,tts,followup)
  dataset=z[[1]]
  tstartp=z[[2]]
  tstopp=z[[3]]
  nstartint=z[[4]]
  startint=z[[5]][1:nstartint]
  nstopint=z[[6]]
  stopint=z[[7]][1:nstopint]
  nitc_start=z[[8]]
  if (nitc_start > 0) {
    itc_start_endpoint=z[[9]][1:nitc_start]
  } else {
    itc_start_endpoint=list()
  }
  nitc_stop=z[[10]]
  if (nitc_stop > 0) {
    itc_stop_endpoint=z[[11]][1:nitc_stop]
  } else {
    itc_stop_endpoint=list()
  }
  z=htcfitter(dataset,ncov,cov_names,nitc_start,itc_start_endpoint,nitc_stop,itc_stop_endpoint,tstartp,tstopp)
  treatcoef_tc2=z[[1]]
  treatse_tc2=z[[2]]
  fit_htc=z[[3]]
  cov_names2=z[[4]]
  list(fit_htc,nitc_start,itc_start_endpoint,nitc_stop,itc_stop_endpoint,nstartint,startint,nstopint,stopint,cov_names2,nperson,numevents,medianfollowup)
}
#------------------------------------------------------
# END: HYBRID TREATMENT CHOICE MODEL FUNCTION
#------------------------------------------------------





#------------------------------------------------------
# PARAMETRIC TREATMENT CHOICE MODEL (PTC)
#------------------------------------------------------
ptc=function(dataset,ncov,cov_names,maxfollow,nmaxint,interval_width,min_future_events) {
  z=treatinit(dataset,nmaxint)
  dataset=z[[1]]
  obsnum=z[[2]]
  nperson=z[[3]]
  maxobs=z[[4]]
  numevents=z[[5]]
  start_times=z[[6]]
  stop_times=z[[7]]
  followup=z[[8]]
  tti=z[[9]]
  tts=z[[10]]
  medianfollowup=z[[11]]
  nitc_start=0
  nitc_stop=0
  itc_start_endpoint=list()
  itc_stop_endpoint=list()
  z=addtc(dataset,ncov,maxfollow,start_times,stop_times,min_future_events,numevents,nperson,nmaxint,maxobs,interval_width,
         nitc_start,itc_start_endpoint,nitc_stop,itc_stop_endpoint,tti,tts,followup)
  dataset=z[[1]]
  tstartp=z[[2]]
  tstopp=z[[3]]
  nstartint=z[[4]]
  startint=z[[5]][1:nstartint]
  nstopint=z[[6]]
  stopint=z[[7]][1:nstopint]
  z=ptcfitter(dataset,ncov,cov_names)
  treatcoef_tc1=z[[1]]
  treatse_tc1=z[[2]]
  fit_ptc=z[[3]]
  cov_names1=z[[4]]
  list(fit_ptc,nstartint,startint,nstopint,stopint,cov_names1,nperson,numevents,medianfollowup)
}
#------------------------------------------------------
# END: PARAMETRIC TREATMENT CHOICE MODEL FUNCTION
#------------------------------------------------------


#------------------------------------------------------
# ORDINARY COX MODEL (COX)
#------------------------------------------------------
cox=function(dataset1,ncov,cov_names,nmaxint) {
  z=treatinit(dataset1,nmaxint)
  nperson=z[[3]]
  numevents=z[[5]]
  medianfollowup=z[[11]]
  z=coxfitter(dataset1,ncov,cov_names)
  treatcoef_cox1=z[[1]]
  treatse_cox1=z[[2]]
  fit_cox=z[[3]]
  list(fit_cox,nperson,numevents,medianfollowup)
}
#------------------------------------------------------
# END: ORDINARY COX MODEL FUNCTION
#------------------------------------------------------




tc=function(type='PTC',dataset,cov_names,maxfollow=100,nmaxint=80,interval_width=.1,min_exp_events=50,min_future_events=50,nitc_fixed=0,n_start_fixed=10,n_stop_fixed=10,interval_stop_beginning=1.1) {
#-----------------------------------------------------------------------
# Function to fit indicated TC model
#-----------------------------------------------------------------------
#cat('----------------------------------','\n')
#cat('start tc with type=',type,'\n')
#library(survival)
ncov=length(cov_names)
#cat('ncov=',ncov,'\n')
if (type=="ITC") {
#  cat('call itc','\n')
  z=itc(dataset,ncov,cov_names,maxfollow,nmaxint,interval_width,min_exp_events,min_future_events,nitc_fixed,n_start_fixed,n_stop_fixed,interval_stop_beginning)
  fit_itc=z[[1]]
  nitc_start=z[[2]]
  itc_start_endpoint=z[[3]]
  nitc_stop=z[[4]]
  itc_stop_endpoint=z[[5]]
  nstartint=z[[6]]
  startint=z[[7]]
  nstopint=z[[8]]
  stopint=z[[9]]
  cov_names3=z[[10]]
  nperson=z[[11]]
  numevents=z[[12]]
  medianfollowup=z[[13]]
  treatcoef_tc3=stats::coef(fit_itc)[1]
  treatse_tc3=sqrt(stats::vcov(fit_itc)[1,1])

  list(fit_itc,nitc_start,itc_start_endpoint,nitc_stop,itc_stop_endpoint,nstartint,startint,nstopint,stopint,cov_names3,nperson,numevents,medianfollowup)

} else {
if (type=="HTC") {
  z=htc(dataset,ncov,cov_names,maxfollow,nmaxint,interval_width,min_exp_events,min_future_events,nitc_fixed,n_start_fixed,n_stop_fixed,interval_stop_beginning)
  fit_htc=z[[1]]
  nitc_start=z[[2]]
  itc_start_endpoint=z[[3]]
  nitc_stop=z[[4]]
  itc_stop_endpoint=z[[5]]
  nstartint=z[[6]]
  startint=z[[7]]
  nstopint=z[[8]]
  stopint=z[[9]]
  cov_names2=z[[10]]
  nperson=z[[11]]
  numevents=z[[12]]
  medianfollowup=z[[13]]
  treatcoef_tc2=stats::coef(fit_htc)[1]
  treatse_tc2=sqrt(stats::vcov(fit_htc)[1,1])

  list(fit_htc,nitc_start,itc_start_endpoint,nitc_stop,itc_stop_endpoint,nstartint,startint,nstopint,stopint,cov_names2,nperson,numevents,medianfollowup)

} else {
if (type=="PTC") {
  z=ptc(dataset,ncov,cov_names,maxfollow,nmaxint,interval_width,min_future_events)
  fit_ptc=z[[1]]
  nstartint=z[[2]]
  startint=z[[3]]
  nstopint=z[[4]]
  stopint=z[[5]]
  cov_names1=z[[6]]
  nperson=z[[7]]
  numevents=z[[8]]
  medianfollowup=z[[9]]
  nitc_start=0
  itc_start_endpoint=list()
  nitc_stop=0
  itc_stop_endpoint=list()
  treatcoef_tc1=stats::coef(fit_ptc)[1]
  treatse_tc1=sqrt(stats::vcov(fit_ptc)[1,1])

  list(fit_ptc,nitc_start,itc_start_endpoint,nitc_stop,itc_stop_endpoint,nstartint,startint,nstopint,stopint,cov_names1,nperson,numevents,medianfollowup)

} else {
if (type=="COX") {
  nmaxint=0
  z=cox(dataset,ncov,cov_names,nmaxint)
  fit_cox=z[[1]]
  nperson=z[[2]]
  numevents=z[[3]]
  medianfollowup=z[[4]]
  nstartint=0
  startint=list()
  nstopint=0
  stopint=list()
  nitc_start=0
  itc_start_endpoint=list()
  nitc_stop=0
  itc_stop_endpoint=list()
  treatcoef_cox1=stats::coef(fit_cox)[1]
  treatse_cox1=sqrt(stats::vcov(fit_cox)[1,1])

  list(fit_cox,nitc_start,itc_start_endpoint,nitc_stop,itc_stop_endpoint,nstartint,startint,nstopint,stopint,cov_names,nperson,numevents,medianfollowup)
}
}
}
}
}
#-----------------------------------------------------------------------
#  END TC Function to fit indicated TC model
#-----------------------------------------------------------------------
#-----------------------------------------------------------------------

#---------------------------------------------
# END: functions of Package tccox
#---------------------------------------------
