library(shiny)
TAN_findIGcutoff<-function(n=30,alloc.ratio=1,
                           prior.mean1=1/3,prior.k1=1,
                           prior.mean2=1/3,prior.k2=1,
                           prior.a1=1,prior.b1=1,
                           prior.a2=1,prior.b2=1,
                           ssq1=1,ssq2=1,
                           control.mean=0,
                           rate=0.8,CT=0.25,
                           seed.num=369,nsim=10000,
                           stop.criterion=10^-3,direction='Greater'){
  ####prior####
  ###mu1~N(prior.mean1,prior.k1*1/var)
  #### 1/var~gamma(prior.a1,rate=prior.b1)
  ####
  cutoff<-NA
  set.seed(seed.num)
  S=1
  n1=round(n*alloc.ratio/(1+alloc.ratio))
  n2=n-n1
  upper=abs(control.mean)+1
  lower=-abs(control.mean)-1
  ssq1=n1*ssq1 ###transform to sum of squares n
  ssq2=n2*ssq2 ###transform to sum of squares n
  ########control.group####
  mupost_2=(prior.k2*prior.mean2+n2*control.mean)/(prior.k2+n2)
  kpost_2=prior.k2+n2
  apost_2=prior.a2+n2/2
  bpost_2=prior.b2+1/2*ssq2+prior.k2*n2*(control.mean-prior.mean2)^2/(2*(prior.k2+n2))
  y=rt(nsim,df=2*apost_2)*sqrt(bpost_2/(apost_2*kpost_2))+mupost_2

  while(S==1){
    temp=(upper+lower)/2
    mupost_temp=(prior.k1*prior.mean1+n1*temp)/(prior.k1+n1)
    kpost_temp=prior.k1+n1
    apost_temp=prior.a1+n1/2
    bpost_temp=prior.b1+1/2*ssq1+prior.k1*n1*(temp-prior.mean1)^2/(2*(prior.k1+n1))
    xtemp=rt(nsim,df=2*apost_temp)*sqrt(bpost_temp/(apost_temp*kpost_temp))+mupost_temp

    mupost_upper=(prior.k1*prior.mean1+n1*upper)/(prior.k1+n1)
    kpost_upper=prior.k1+n1
    apost_upper=prior.a1+n1/2
    bpost_upper=prior.b1+1/2*ssq1+prior.k1*n1*(upper-prior.mean1)^2/(2*(prior.k1+n1))
    xupper=rt(nsim,df=2*apost_upper)*sqrt(bpost_upper/(apost_upper*kpost_upper))+mupost_upper

    mupost_lower=(prior.k1*prior.mean1+n1*lower)/(prior.k1+n1)
    kpost_lower=prior.k1+n1
    apost_lower=prior.a1+n1/2
    bpost_lower=prior.b1+1/2*ssq1+prior.k1*n1*(lower-prior.mean1)^2/(2*(prior.k1+n1))
    xlower=rt(nsim,df=2*apost_lower)*sqrt(bpost_lower/(apost_lower*kpost_lower))+mupost_lower
    if(direction=='Greater'){
      probupper=sum(xupper>=y+CT)/nsim
      problower=sum(xlower>=y+CT)/nsim
      probtemp=sum(xtemp>=y+CT)/nsim
      if(problower>=rate){
        if(abs(problower-rate)<=stop.criterion){
          S=0
          cutoff=lower
        }
        upper=lower
        lower=-2*abs(lower)

      }else if(probupper<=rate){
        if(abs(probupper-rate)<=stop.criterion){
          S=0
          cutoff=upper
        }
        upper=2*abs(upper)
        lower=upper
      }else{
        if(abs(probtemp-rate)<=stop.criterion){
          S=0
          cutoff=temp
        }
        if(probtemp<rate){lower=temp}
        if(probtemp>rate){upper=temp}
      }
    }
    if(direction=='Less'){
      probupper=sum(xupper<=y+CT)/nsim
      problower=sum(xlower<=y+CT)/nsim
      probtemp=sum(xtemp<=y+CT)/nsim
      if(problower<=rate){
        if(abs(problower-rate)<=stop.criterion){
          S=0
          cutoff=lower
        }
        upper=lower
        lower=-2*abs(lower)

      }else if(probupper>=rate){
        if(abs(probupper-rate)<=stop.criterion){
          S=0
          cutoff=upper
        }
        upper=2*abs(upper)
        lower=upper
      }else{
        if(abs(probtemp-rate)<=stop.criterion){
          S=0
          cutoff=temp
        }
        if(probtemp>rate){lower=temp}
        if(probtemp<rate){upper=temp}
      }
    }
  }

  return(cutoff)

}

TAN_Normal_Cutoff<-function(n=30,alloc.ratio=1,
                            prior.mean1=1/3,prior.sd1=1,prior.k1=1,
                            prior.mean2=1/3,prior.sd2=1,prior.k2=1,
                            prior.a1=1,prior.b1=1,
                            prior.a2=1,prior.b2=1,
                            sd1=1,sd2=1,
                            ssq1=1,ssq2=1,
                            control.mean=0,
                            CT1.go=0.25,
                            false.go.CT1=TRUE,FGR.CT1=0.25,
                            CT1.nogo=0.25,
                            false.nogo.CT1=TRUE,FNGR.CT1=0.25,
                            CT2.go=0.3,
                            false.go.CT2=TRUE, FGR.CT2=0.5,
                            CT2.nogo=0.3,
                            false.nogo.CT2=TRUE, FNGR.CT2=0.5,
                            method='Bayesian',direction='Greater',
                            fix.var=TRUE,seed.num=369,nsim=10000,stop.criterion=10^-3,noninfo=TRUE,
                            logic.go='and',logic.nogo='or'){
  flag=rep(0,4)
  overlap.flag=0
  n1=round(n*alloc.ratio/(1+alloc.ratio))
  n2=n-n1
  sd.n1=sd1/sqrt(n1)
  sd.n2=sd2/sqrt(n2)
  sd.n=sqrt(sd1^2/n1+sd2^2/n2)
  est1.go=NA
  est2.go=NA
  est1.nogo=NA
  est2.nogo=NA
  if(is.na(CT1.go)){
    false.go.CT1=FALSE
  }
  if(is.na(CT1.nogo)){
    false.nogo.CT1=FALSE
  }
  if(is.na(CT2.go)){
    false.go.CT2=FALSE
  }
  if(is.na(CT2.nogo)){
    false.nogo.CT2=FALSE
  }

  if(direction=='Greater'){
    if(method=='Bayesian'&fix.var==TRUE){
      sd.post=sqrt(1/((1/prior.sd1)^2+(1/sd.n1)^2)+1/((1/prior.sd2)^2+(1/sd.n2)^2))
      control.post=(sd.n2^2/(prior.sd2^2+sd.n2^2))*prior.mean2+(prior.sd2^2/(prior.sd2^2+sd.n2^2))*control.mean
      if(false.go.CT1==TRUE){
        temp=qnorm(1-FGR.CT1,mean=CT1.go+control.post,sd=sd.post)
        est1.go=(temp-(sd.n1^2/(prior.sd1^2+sd.n1^2))*prior.mean1)/(prior.sd1^2/(prior.sd1^2+sd.n1^2))-control.mean
      }
      if(false.go.CT2==TRUE){
        temp=qnorm(1-FGR.CT2,mean=CT2.go+control.post,sd=sd.post)
        est2.go=(temp-(sd.n1^2/(prior.sd1^2+sd.n1^2))*prior.mean1)/(prior.sd1^2/(prior.sd1^2+sd.n1^2))-control.mean
      }

      if(false.nogo.CT1==TRUE){
        temp=qnorm(FNGR.CT1,mean=CT1.nogo+control.post,sd=sd.post)
        est1.nogo=(temp-(sd.n1^2/(prior.sd1^2+sd.n1^2))*prior.mean1)/(prior.sd1^2/(prior.sd1^2+sd.n1^2))-control.mean

      }
      if(false.nogo.CT2==TRUE){
        temp=qnorm(FNGR.CT2,mean=CT2.nogo+control.post,sd=sd.post)
        est2.nogo=(temp-(sd.n1^2/(prior.sd1^2+sd.n1^2))*prior.mean1)/(prior.sd1^2/(prior.sd1^2+sd.n1^2))-control.mean
      }



    }
    if(method=='Bayesian'&fix.var==FALSE){
      if(noninfo==TRUE){
        prior.mean1=0
        prior.k1=0
        prior.a1=-1/2
        prior.b1=0
        prior.mean2=0
        prior.k2=0
        prior.a2=-1/2
        prior.b2=0
      }
      if(false.go.CT1==TRUE){
        mean1.go=TAN_findIGcutoff(n=n,alloc.ratio=alloc.ratio,
                                  prior.mean1=prior.mean1,prior.k1=prior.k1,
                                  prior.mean2=prior.mean2,prior.k2=prior.k2,
                                  prior.a1=prior.a1,prior.b1=prior.b1,
                                  prior.a2=prior.a2,prior.b2=prior.b2,
                                  ssq1=ssq1,ssq2=ssq2,
                                  control.mean=control.mean,
                                  rate=1-FGR.CT1,CT=CT1.go,
                                  seed.num=seed.num ,nsim=nsim,
                                  stop.criterion=stop.criterion,direction="Greater")
        est1.go=mean1.go-control.mean
      }
      if(false.go.CT2==TRUE){
        mean2.go=TAN_findIGcutoff(n=n,alloc.ratio=alloc.ratio,
                                  prior.mean1=prior.mean1,prior.k1=prior.k1,
                                  prior.mean2=prior.mean2,prior.k2=prior.k2,
                                  prior.a1=prior.a1,prior.b1=prior.b1,
                                  prior.a2=prior.a2,prior.b2=prior.b2,
                                  ssq1=ssq1,ssq2=ssq2,
                                  control.mean=control.mean,
                                  rate=1-FGR.CT2,CT=CT2.go,
                                  seed.num=seed.num ,nsim=nsim,
                                  stop.criterion=stop.criterion,direction="Greater")
        est2.go=mean2.go-control.mean
      }

      if(false.nogo.CT1==TRUE){
        mean1.nogo=TAN_findIGcutoff(n=n,alloc.ratio=alloc.ratio,
                                    prior.mean1=prior.mean1,prior.k1=prior.k1,
                                    prior.mean2=prior.mean2,prior.k2=prior.k2,
                                    prior.a1=prior.a1,prior.b1=prior.b1,
                                    prior.a2=prior.a2,prior.b2=prior.b2,
                                    ssq1=ssq1,ssq2=ssq2,
                                    control.mean=control.mean,
                                    rate=1-FNGR.CT1,CT=CT1.nogo,
                                    seed.num=seed.num ,nsim=nsim,
                                    stop.criterion=stop.criterion,direction="Less")
        est1.nogo=mean1.nogo-control.mean
      }
      if(false.nogo.CT2==TRUE){
        mean2.nogo=TAN_findIGcutoff(n=n,alloc.ratio=alloc.ratio,
                                    prior.mean1=prior.mean1,prior.k1=prior.k1,
                                    prior.mean2=prior.mean2,prior.k2=prior.k2,
                                    prior.a1=prior.a1,prior.b1=prior.b1,
                                    prior.a2=prior.a2,prior.b2=prior.b2,
                                    ssq1=ssq1,ssq2=ssq2,
                                    control.mean=control.mean,
                                    rate=1-FNGR.CT2,CT=CT2.nogo,
                                    seed.num=seed.num ,nsim=nsim,
                                    stop.criterion=stop.criterion,direction="Less")
        est2.nogo=mean2.nogo-control.mean
      }

    }
    if(method=='Frequentist'){
      if(false.go.CT1==TRUE){
        #est1.go=qnorm(1-FGR.CT1,mean=CT1-control.mean,sd=sd.n)
        est1.go=qnorm(1-FGR.CT1,mean=CT1.go,sd=sd.n)
      }
      if(false.go.CT2==TRUE){
        #est2.go=qnorm(1-FGR.CT2,mean=CT2-control.mean,sd=sd.n)
        est2.go=qnorm(1-FGR.CT2,mean=CT2.go,sd=sd.n)
      }

      if(false.nogo.CT1==TRUE){
        #est1.nogo=qnorm(FNGR.CT1,mean=CT1-control.mean,sd=sd.n)
        est1.nogo=qnorm(FNGR.CT1,mean=CT1.nogo,sd=sd.n)

      }
      if(false.nogo.CT2==TRUE){
        #est2.nogo=qnorm(FNGR.CT2,mean=CT2-control.mean,sd=sd.n)
        est2.nogo=qnorm(FNGR.CT2,mean=CT2.nogo,sd=sd.n)
      }


    }
    if(any(is.na(c(est1.go,est2.go)))){logic.go='and'}
    if(any(is.na(c(est1.nogo,est2.nogo)))){logic.nogo='and'}
    if(logic.go=='and'){
      go_cutoff=max(est1.go,est2.go,na.rm=TRUE)
    }
    if(logic.go=='or'){
      go_cutoff=min(est1.go,est2.go,na.rm=TRUE)
    }
    if(logic.nogo=='and')
    {
      nogo_cutoff=min(est1.nogo,est2.nogo,na.rm=TRUE)
    }
    if(logic.nogo=='or')
    {
      nogo_cutoff=max(est1.nogo,est2.nogo,na.rm=TRUE)
    }


    if(go_cutoff>=nogo_cutoff){return(list(cutoff=c(go_cutoff,nogo_cutoff),flag=flag,overlap=overlap.flag))}
    else{
      overlap.flag=1
      return(list(cutoff=c(go_cutoff,nogo_cutoff),flag=flag,overlap=overlap.flag))}
    if(all(flag==0)==FALSE){return(list(cutoff=c(NA,NA),flag=flag,overlap=overlap.flag))}

  }

  if(direction=='Less'){

    if(method=='Bayesian'&fix.var==TRUE){
      sd.post=sqrt(1/((1/prior.sd1)^2+(1/sd.n1)^2)+1/((1/prior.sd2)^2+(1/sd.n2)^2))
      control.post=(sd.n2^2/(prior.sd2^2+sd.n2^2))*prior.mean2+(prior.sd2^2/(prior.sd2^2+sd.n2^2))*control.mean
      if(false.go.CT1==TRUE){
        temp=qnorm(FGR.CT1,mean=CT1.go+control.post,sd=sd.post)
        est1.go=(temp-(sd.n1^2/(prior.sd1^2+sd.n1^2))*prior.mean1)/(prior.sd1^2/(prior.sd1^2+sd.n1^2))-control.mean
      }
      if(false.go.CT2==TRUE){
        temp=qnorm(FGR.CT2,mean=CT2.go+control.post,sd=sd.post)
        est2.go=(temp-(sd.n1^2/(prior.sd1^2+sd.n1^2))*prior.mean1)/(prior.sd1^2/(prior.sd1^2+sd.n1^2))-control.mean
      }

      if(false.nogo.CT1==TRUE){
        temp=qnorm(1-FNGR.CT1,mean=CT1.nogo+control.post,sd=sd.post)
        est1.nogo=(temp-(sd.n1^2/(prior.sd1^2+sd.n1^2))*prior.mean1)/(prior.sd1^2/(prior.sd1^2+sd.n1^2))-control.mean

      }
      if(false.nogo.CT2==TRUE){
        temp=qnorm(1-FNGR.CT2,mean=CT2.nogo+control.post,sd=sd.post)
        est2.nogo=(temp-(sd.n1^2/(prior.sd1^2+sd.n1^2))*prior.mean1)/(prior.sd1^2/(prior.sd1^2+sd.n1^2))-control.mean
      }



    }

    if(method=='Bayesian'&fix.var==FALSE){
      if(noninfo==TRUE){
        prior.mean1=0
        prior.k1=0
        prior.a1=-1/2
        prior.b1=0
        prior.mean2=0
        prior.k2=0
        prior.a2=-1/2
        prior.b2=0
      }
      if(false.go.CT1==TRUE){
        mean1.go=TAN_findIGcutoff(n=n,alloc.ratio=alloc.ratio,
                                  prior.mean1=prior.mean1,prior.k1=prior.k1,
                                  prior.mean2=prior.mean2,prior.k2=prior.k2,
                                  prior.a1=prior.a1,prior.b1=prior.b1,
                                  prior.a2=prior.a2,prior.b2=prior.b2,
                                  ssq1=ssq1,ssq2=ssq2,
                                  control.mean=control.mean,
                                  rate=1-FGR.CT1,CT=CT1.go,
                                  seed.num=seed.num ,nsim=nsim,
                                  stop.criterion=stop.criterion,direction="Less")
        est1.go=mean1.go-control.mean
      }
      if(false.go.CT2==TRUE){
        mean2.go=TAN_findIGcutoff(n=n,alloc.ratio=alloc.ratio,
                                  prior.mean1=prior.mean1,prior.k1=prior.k1,
                                  prior.mean2=prior.mean2,prior.k2=prior.k2,
                                  prior.a1=prior.a1,prior.b1=prior.b1,
                                  prior.a2=prior.a2,prior.b2=prior.b2,
                                  ssq1=ssq1,ssq2=ssq2,
                                  control.mean=control.mean,
                                  rate=1-FGR.CT2,CT=CT2.go,
                                  seed.num=seed.num ,nsim=nsim,
                                  stop.criterion=stop.criterion,direction="Less")
        est2.go=mean2.go-control.mean
      }

      if(false.nogo.CT1==TRUE){
        mean1.nogo=TAN_findIGcutoff(n=n,alloc.ratio=alloc.ratio,
                                    prior.mean1=prior.mean1,prior.k1=prior.k1,
                                    prior.mean2=prior.mean2,prior.k2=prior.k2,
                                    prior.a1=prior.a1,prior.b1=prior.b1,
                                    prior.a2=prior.a2,prior.b2=prior.b2,
                                    ssq1=ssq1,ssq2=ssq2,
                                    control.mean=control.mean,
                                    rate=1-FNGR.CT1,CT=CT1.nogo,
                                    seed.num=seed.num ,nsim=nsim,
                                    stop.criterion=stop.criterion,direction="Greater")
        est1.nogo=mean1.nogo-control.mean
      }
      if(false.nogo.CT2==TRUE){
        mean2.nogo=TAN_findIGcutoff(n=n,alloc.ratio=alloc.ratio,
                                    prior.mean1=prior.mean1,prior.k1=prior.k1,
                                    prior.mean2=prior.mean2,prior.k2=prior.k2,
                                    prior.a1=prior.a1,prior.b1=prior.b1,
                                    prior.a2=prior.a2,prior.b2=prior.b2,
                                    ssq1=ssq1,ssq2=ssq2,
                                    control.mean=control.mean,
                                    rate=1-FNGR.CT2,CT=CT2.nogo,
                                    seed.num=seed.num ,nsim=nsim,
                                    stop.criterion=stop.criterion,direction="Greater")
        est2.nogo=mean2.nogo-control.mean
      }

    }

    if(method=='Frequentist'){
      if(false.go.CT1==TRUE){
        #est1.go=qnorm(FGR.CT1,mean=CT1-control.mean,sd=sd.n)
        est1.go=qnorm(FGR.CT1,mean=CT1.go,sd=sd.n)
      }
      if(false.go.CT2==TRUE){
        #est2.go=qnorm(FGR.CT2,mean=CT2-control.mean,sd=sd.n)
        est2.go=qnorm(FGR.CT2,mean=CT2.go,sd=sd.n)
      }

      if(false.nogo.CT1==TRUE){
        #est1.nogo=qnorm(1-FNGR.CT1,mean=CT1-control.mean,sd=sd.n)
        est1.nogo=qnorm(1-FNGR.CT1,mean=CT1.nogo,sd=sd.n)

      }
      if(false.nogo.CT2==TRUE){
        #est2.nogo=qnorm(1-FNGR.CT2,mean=CT2-control.mean,sd=sd.n)
        est2.nogo=qnorm(1-FNGR.CT2,mean=CT2.nogo,sd=sd.n)
      }


    }

    if(any(is.na(c(est1.go,est2.go)))){logic.go='and'}
    if(any(is.na(c(est1.nogo,est2.nogo)))){logic.nogo='and'}
    if(logic.go=='and'){
      go_cutoff=min(est1.go,est2.go,na.rm=TRUE)
    }
    if(logic.go=='or'){
      go_cutoff=max(est1.go,est2.go,na.rm=TRUE)
    }
    if(logic.nogo=='and')
    {
      nogo_cutoff=max(est1.nogo,est2.nogo,na.rm=TRUE)
    }
    if(logic.nogo=='or')
    {
      nogo_cutoff=min(est1.nogo,est2.nogo,na.rm=TRUE)
    }

    if(go_cutoff<=nogo_cutoff){return(list(cutoff=c(go_cutoff,nogo_cutoff),flag=flag,overlap=overlap.flag))}
    else{
      overlap.flag=1
      return(list(cutoff=c(go_cutoff,nogo_cutoff),flag=flag,overlap=overlap.flag))}
    if(all(flag==0)==FALSE){return(list(cutoff=c(NA,NA),flag=flag,overlap=overlap.flag))}


  }

}

TAN_Bayesian<-function(n=c(10,100),alloc.ratio=1,
                       prior.mean1=1/3,prior.sd1=1,prior.k1=1,
                       prior.mean2=1/3,prior.sd2=1,prior.k2=1,
                       prior.a1=1,prior.b1=1,
                       prior.a2=1,prior.b2=1,
                       control.mean=0,
                       mean=c(0.0,0.9),
                       sd1=1,sd2=1,
                       CT1.go=0.25,
                       false.go.CT1=TRUE,FGR.CT1=0.25,
                       CT1.nogo=0.25,
                       false.nogo.CT1=TRUE,FNGR.CT1=0.25,
                       CT2.go=0.3,
                       false.go.CT2=TRUE, FGR.CT2=0.5,
                       CT2.nogo=0.3,
                       false.nogo.CT2=TRUE, FNGR.CT2=0.5,
                       overlap.option='GO', direction='Greater',
                       fix.var=TRUE, noninfo=TRUE,
                       seed.num=369,nsim=10000,n_repeat=1000,
                       logic.go='and',logic.nogo='or'){
  set.seed(seed.num)
  
  go_prob<-matrix(NA,ncol=length(n),nrow=length(mean))
  nogo_prob<-matrix(NA,ncol=length(n),nrow=length(mean))
  inconclusive_prob<-matrix(NA,ncol=length(n),nrow=length(mean))
  overlap.flag=rep(0,length(n))
  
  for (ii in 1:length(mean)){
    for (ii2 in 1:length(n)){
      
      treatment.mean=control.mean+mean[ii]
      i1=round(n[ii2]*alloc.ratio/(1+alloc.ratio))
      i2=n[ii2]-i1
        
      x1=matrix(rnorm(n_repeat*i1,treatment.mean,sd1),nrow=n_repeat)
      x2=matrix(rnorm(n_repeat*i2,control.mean,sd2),nrow=n_repeat)
        
      mean_x1=apply(x1,1,mean)
      mean_x2=apply(x2,1,mean)
      if(fix.var==FALSE){
        ssq_x1=apply(x1,1,function(x){sum(x^2)})
        ssq_x2=apply(x2,1,function(x){sum(x^2)})
      }
      
      if(is.na(CT1.go)) false.go.CT1=FALSE
      if(is.na(CT1.nogo)) false.nogo.CT1=FALSE
      if(is.na(CT2.go)) false.go.CT2=FALSE
      if(is.na(CT2.nogo)) false.nogo.CT2=FALSE
      
      if(false.go.CT1==TRUE) pp_go_1 = rep(0,n_repeat)
      if(false.nogo.CT1==TRUE) pp_nogo_1 = rep(0,n_repeat)
      if(false.go.CT2==TRUE) pp_go_2 = rep(0,n_repeat)
      if(false.nogo.CT2==TRUE) pp_nogo_2 = rep(0,n_repeat)
      
      for (i in 1:n_repeat){
          
        if(fix.var==TRUE){
            
          sd.post_1=sqrt(1/((1/prior.sd1)^2+i1/sd1^2))
          sd.post_2=sqrt(1/((1/prior.sd2)^2+i2/sd2^2))
          mupost_1=(1/(1+prior.sd1^2/sd1^2*i1))*prior.mean1+(1/(1+sd1^2/i1/prior.sd1^2))*mean_x1[i]
          mupost_2=(1/(1+prior.sd2^2/sd2^2*i2))*prior.mean2+(1/(1+sd2^2/i2/prior.sd2^2))*mean_x2[i]
            
          if(direction=='Greater'){ # closed form
            if(false.go.CT1==TRUE) pp_go_1[i] = 1-pnorm(CT1.go,mean=mupost_1-mupost_2,sd=(sd.post_1^2+sd.post_2^2)^(1/2))
            if(false.nogo.CT1==TRUE) pp_nogo_1[i] = 1-pnorm(CT1.nogo,mean=mupost_1-mupost_2,sd=(sd.post_1^2+sd.post_2^2)^(1/2))
            if(false.go.CT2==TRUE) pp_go_2[i] = 1-pnorm(CT2.go,mean=mupost_1-mupost_2,sd=(sd.post_1^2+sd.post_2^2)^(1/2))
            if(false.nogo.CT2==TRUE) pp_nogo_2[i] = 1-pnorm(CT2.nogo,mean=mupost_1-mupost_2,sd=(sd.post_1^2+sd.post_2^2)^(1/2))
          }
          if(direction=='Less'){
            if(false.go.CT1==TRUE) pp_go_1[i] = pnorm(CT1.go,mean=mupost_1-mupost_2,sd=(sd.post_1^2+sd.post_2^2)^(1/2))
            if(false.nogo.CT1==TRUE) pp_nogo_1[i] = pnorm(CT1.nogo,mean=mupost_1-mupost_2,sd=(sd.post_1^2+sd.post_2^2)^(1/2))
            if(false.go.CT2==TRUE) pp_go_2[i] = pnorm(CT2.go,mean=mupost_1-mupost_2,sd=(sd.post_1^2+sd.post_2^2)^(1/2))
            if(false.nogo.CT2==TRUE) pp_nogo_2[i] = pnorm(CT2.nogo,mean=mupost_1-mupost_2,sd=(sd.post_1^2+sd.post_2^2)^(1/2))
          }
        }
          
        if(fix.var==FALSE){
            
          if(noninfo==TRUE){
            prior.mean1=0
            prior.k1=0
            prior.a1=-1/2
            prior.b1=0
            prior.mean2=0
            prior.k2=0
            prior.a2=-1/2
            prior.b2=0
          }
            
          mupost_1=(prior.k1*prior.mean1+i1*mean_x1[i])/(prior.k1+i1)
          kpost_1=prior.k1+i1
          apost_1=prior.a1+i1/2
          bpost_1=prior.b1+1/2*ssq_x1[i]+prior.k1*i1*(mean_x1[i]-prior.mean1)^2/(2*(prior.k1+i1))
            
          mupost_2=(prior.k2*prior.mean2+i2*mean_x2[i])/(prior.k2+i2)
          kpost_2=prior.k2+i2
          apost_2=prior.a2+i2/2
          bpost_2=prior.b2+1/2*ssq_x2[i]+prior.k2*i2*(mean_x2[i]-prior.mean2)^2/(2*(prior.k2+i2))
            
          temp1 = rt(nsim,df=2*apost_1)*sqrt(bpost_1/(apost_1*kpost_1))+mupost_1
          temp2 = rt(nsim,df=2*apost_2)*sqrt(bpost_2/(apost_2*kpost_2))+mupost_2
            
          if(direction=='Greater'){
            if(false.go.CT1==TRUE) pp_go_1[i] = mean(temp1-temp2>=CT1.go)
            if(false.nogo.CT1==TRUE) pp_nogo_1[i] = mean(temp1-temp2>=CT1.nogo)
            if(false.go.CT2==TRUE) pp_go_2[i] = mean(temp1-temp2>=CT2.go)
            if(false.nogo.CT2==TRUE) pp_nogo_2[i] = mean(temp1-temp2>=CT2.nogo)
          }
          if(direction=='Less'){
            if(false.go.CT1==TRUE) pp_go_1[i] = mean(temp1-temp2<=CT1.go)
            if(false.nogo.CT1==TRUE) pp_nogo_1[i] = mean(temp1-temp2<=CT1.nogo)
            if(false.go.CT2==TRUE) pp_go_2[i] = mean(temp1-temp2<=CT2.go)
            if(false.nogo.CT2==TRUE) pp_nogo_2[i] = mean(temp1-temp2<=CT2.nogo)
          }
        }
      }
      
      if(false.go.CT1==TRUE & false.go.CT2==TRUE){
        if(logic.go=='and')
          go_prob[ii,ii2] = mean(pp_go_1>=1-FGR.CT1 & pp_go_2>=1-FGR.CT2)
        if(logic.go=='or')
          go_prob[ii,ii2] = mean(pp_go_1>=1-FGR.CT1 | pp_go_2>=1-FGR.CT2)
      } else if(false.go.CT1==TRUE){
        go_prob[ii,ii2] = mean(pp_go_1>=1-FGR.CT1)
      } else if(false.go.CT2==TRUE){
        go_prob[ii,ii2] = mean(pp_go_2>=1-FGR.CT2)
      }
      
      if(false.nogo.CT1==TRUE & false.nogo.CT2==TRUE){
        if(logic.nogo=='and')
          nogo_prob[ii,ii2] = mean(pp_nogo_1<FNGR.CT1 & pp_nogo_2<FNGR.CT2)
        if(logic.nogo=='or')
          nogo_prob[ii,ii2] = mean(pp_nogo_1<FNGR.CT1 | pp_nogo_2<FNGR.CT2)
      } else if(false.nogo.CT1==TRUE){
        nogo_prob[ii,ii2] = mean(pp_nogo_1<FNGR.CT1)
      } else if(false.nogo.CT2==TRUE){
        nogo_prob[ii,ii2] = mean(pp_nogo_2<FNGR.CT2)
      }
      
      if (go_prob[ii,ii2]+nogo_prob[ii,ii2]>1){
        overlap.flag[ii2] = 1
        if (overlap.option=='GO'){
          nogo_prob[ii,ii2] = 1-go_prob[ii,ii2]
        } else if (overlap.option=='NOGO'){
          go_prob[ii,ii2] = 1-nogo_prob[ii,ii2]
        }
      }
      inconclusive_prob[ii,ii2]=1-go_prob[ii,ii2]-nogo_prob[ii,ii2]
    }
  }
  
  return(list(go_prob,nogo_prob,inconclusive_prob,overlap.flag))
}

Fix_SS_TAN_Normal_Prob<-function(n=100,alloc.ratio=1,
                                 prior.mean1=1/3,prior.sd1=1,prior.k1=1,
                                 prior.mean2=1/3,prior.sd2=1,prior.k2=1,
                                 prior.a1=1,prior.b1=1,
                                 prior.a2=1,prior.b2=1,
                                 control.mean=0,
                                 mean=c(0.0,0.9),
                                 sd1=1,sd2=1,
                                 ssq1=1,ssq2=2,
                                 CT1.go=0.25,
                                 false.go.CT1=TRUE,FGR.CT1=0.25,
                                 CT1.nogo=0.25,
                                 false.nogo.CT1=TRUE,FNGR.CT1=0.25,
                                 CT2.go=0.3,
                                 false.go.CT2=TRUE, FGR.CT2=0.5,
                                 CT2.nogo=0.3,
                                 false.nogo.CT2=TRUE, FNGR.CT2=0.5,
                                 overlap.option='GO',plot.figure=TRUE,
                                 method='Bayesian',direction='Greater',
                                 fix.var=TRUE,noninfo=TRUE,
                                 seed.num=369,nsim=10000,stop.criterion=10^-3,
                                 logic.go='and',logic.nogo='or'){
  meanseq=round(seq(min(mean),max(mean),length=10),3)
  go_prob<-rep(NA,length(mean))
  nogo_prob<-rep(NA,length(mean))
  inconclusive_prob<-rep(NA,length(mean))
  go_prob_plot<-rep(NA,length(meanseq))
  nogo_prob_plot<-rep(NA,length(meanseq))
  inconclusive_prob_plot<-rep(NA,length(meanseq))
  index=1
  unsatisfied.flag=0
  overlap.flag=0
  
  n1=round(n*alloc.ratio/(1+alloc.ratio))
  n2=n-n1
  if (abs(n*alloc.ratio/(1+alloc.ratio)-n1)>10^(-3))
    stop ("Sample size for each arm must be an integer.")
  
  
  sd.n=sqrt(sd1^2/n1+sd2^2/n2)
  
  temp=TAN_Normal_Cutoff(n=n,alloc.ratio=alloc.ratio,
                         prior.mean1=prior.mean1,prior.sd1=prior.sd1,prior.k1=prior.k1,
                         prior.mean2=prior.mean2,prior.sd2=prior.sd2,prior.k2=prior.k2,
                         prior.a1=prior.a1,prior.b1=prior.b1,
                         prior.a2=prior.a2,prior.b2=prior.b2,
                         sd1=sd1,sd2=sd2,control.mean=control.mean,
                         ssq1=ssq1,ssq2=ssq2,
                         CT1.go=CT1.go,
                         false.go.CT1 = false.go.CT1, FGR.CT1=FGR.CT1,
                         CT1.nogo=CT1.nogo,
                         false.nogo.CT1 = false.nogo.CT1, FNGR.CT1 = FNGR.CT1,
                         CT2.go=CT2.go,
                         false.go.CT2 = false.go.CT2, FGR.CT2=FGR.CT2,
                         CT2.nogo=CT2.nogo,
                         false.nogo.CT2 = false.nogo.CT2, FNGR.CT2 = FNGR.CT2,
                         method=method,direction=direction,
                         fix.var=fix.var,noninfo=noninfo,
                         seed.num=seed.num,nsim=nsim,stop.criterion=stop.criterion,
                         logic.go = logic.go,logic.nogo=logic.nogo)

  ###TAN####
  true_go_cutoff<-ifelse(temp$overlap==0,temp$cutoff[1],(overlap.option=='GO')*temp$cutoff[1]+(overlap.option=='NOGO')*(temp$cutoff[2]))
  true_nogo_cutoff<-ifelse(temp$overlap==0,temp$cutoff[2],(overlap.option=='GO')*temp$cutoff[1]+(overlap.option=='NOGO')*(temp$cutoff[2]))

  if(method=='Frequentist'){
    
    if(direction=='Greater'){
      if(all(temp$flag==0)&temp$overlap==0){
        go_prob=1-pnorm(temp$cutoff[1],mean=mean,sd=sd.n)
        nogo_prob=pnorm(temp$cutoff[2],mean=mean,sd=sd.n)
        inconclusive_prob=1-go_prob-nogo_prob
        go_prob_plot=1-pnorm(temp$cutoff[1],mean=meanseq,sd=sd.n)
        nogo_prob_plot=pnorm(temp$cutoff[2],mean=meanseq,sd=sd.n)
        inconclusive_prob_plot=1-go_prob_plot-nogo_prob_plot
      }
      if(all(temp$flag==0)&temp$overlap==1&overlap.option=="GO"){
        go_prob=1-pnorm(temp$cutoff[1],mean=mean,sd=sd.n)
        nogo_prob=pnorm(temp$cutoff[1],mean=mean,sd=sd.n)
        inconclusive_prob=1-go_prob-nogo_prob
        
        go_prob_plot=1-pnorm(temp$cutoff[1],mean=meanseq,sd=sd.n)
        nogo_prob_plot=pnorm(temp$cutoff[1],mean=meanseq,sd=sd.n)
        inconclusive_prob_plot=1-go_prob_plot-nogo_prob_plot
        
        overlap.flag=1
      }
      if(all(temp$flag==0)&temp$overlap==1&overlap.option=="NOGO"){
        go_prob=1-pnorm(temp$cutoff[2],mean=mean,sd=sd.n)
        nogo_prob=pnorm(temp$cutoff[2],mean=mean,sd=sd.n)
        inconclusive_prob=1-go_prob-nogo_prob
        
        go_prob_plot=1-pnorm(temp$cutoff[2],mean=meanseq,sd=sd.n)
        nogo_prob_plot=pnorm(temp$cutoff[2],mean=meanseq,sd=sd.n)
        inconclusive_prob_plot=1-go_prob_plot-nogo_prob_plot
        
        overlap.flag=1
      }
      if(all(temp$flag==0)==FALSE){
        unsatisfied.flag=1
        nogo_prob=rep(NA,length(mean))
        go_prob=rep(NA,length(mean))
        inconclusive_prob=rep(NA,length(mean))
      }
      
    }
    
    if(direction=='Less'){
      if(all(temp$flag==0)&temp$overlap==0){
        go_prob=pnorm(temp$cutoff[1],mean=mean,sd=sd.n)
        nogo_prob=1-pnorm(temp$cutoff[2],mean=mean,sd=sd.n)
        inconclusive_prob=1-go_prob-nogo_prob
        
        go_prob_plot=pnorm(temp$cutoff[1],mean=meanseq,sd=sd.n)
        nogo_prob_plot=1-pnorm(temp$cutoff[2],mean=meanseq,sd=sd.n)
        inconclusive_prob_plot=1-go_prob_plot-nogo_prob_plot
      }
      if(all(temp$flag==0)&temp$overlap==1&overlap.option=="GO"){
        go_prob=pnorm(temp$cutoff[1],mean=mean,sd=sd.n)
        nogo_prob=1-pnorm(temp$cutoff[1],mean=mean,sd=sd.n)
        inconclusive_prob=1-go_prob-nogo_prob
        
        go_prob_plot=pnorm(temp$cutoff[1],mean=meanseq,sd=sd.n)
        nogo_prob_plot=1-pnorm(temp$cutoff[1],mean=meanseq,sd=sd.n)
        inconclusive_prob_plot=1-go_prob_plot-nogo_prob_plot
        
        overlap.flag=1
      }
      if(all(temp$flag==0)&temp$overlap==1&overlap.option=="NOGO"){
        go_prob_plot=pnorm(temp$cutoff[2],mean=meanseq,sd=sd.n)
        nogo_prob_plot=1-pnorm(temp$cutoff[2],mean=meanseq,sd=sd.n)
        inconclusive_prob_plot=1-go_prob-nogo_prob_plot
        
        go_prob=pnorm(temp$cutoff[2],mean=mean,sd=sd.n)
        nogo_prob=1-pnorm(temp$cutoff[2],mean=mean,sd=sd.n)
        inconclusive_prob=1-go_prob-nogo_prob
        
        overlap.flag=1
      }
      if(all(temp$flag==0)==FALSE){
        unsatisfied.flag=1
        nogo_prob=rep(NA,length(mean))
        go_prob=rep(NA,length(mean))
        inconclusive_prob=rep(NA,length(mean))
      }
    }
  }

  if(method=='Bayesian'){
    
    n_repeat=1000
    
    temp2 = TAN_Bayesian(n=n,alloc.ratio=alloc.ratio,
                        prior.mean1=prior.mean1,prior.sd1=prior.sd1,prior.k1=prior.k1,
                        prior.mean2=prior.mean2,prior.sd2=prior.sd2,prior.k2=prior.k2,
                        prior.a1=prior.a1,prior.b1=prior.b1,
                        prior.a2=prior.a2,prior.b2=prior.b2,
                        control.mean=control.mean,
                        mean=mean,
                        sd1=sd1,sd2=sd2,
                        CT1.go=CT1.go,
                        false.go.CT1=false.go.CT1,FGR.CT1=FGR.CT1,
                        CT1.nogo=CT1.nogo,
                        false.nogo.CT1=false.nogo.CT1,FNGR.CT1=FNGR.CT1,
                        CT2.go=CT2.go,
                        false.go.CT2=false.go.CT2, FGR.CT2=FGR.CT2,
                        CT2.nogo=CT2.nogo,
                        false.nogo.CT2=false.nogo.CT2, FNGR.CT2=FNGR.CT2,
                        overlap.option=overlap.option, direction=direction,
                        fix.var=fix.var, noninfo=noninfo,
                        seed.num=seed.num,nsim=nsim,n_repeat=n_repeat,
                        logic.go=logic.go,logic.nogo=logic.nogo)
    
    go_prob = as.numeric(temp2[[1]])
    nogo_prob = as.numeric(temp2[[2]])
    inconclusive_prob = as.numeric(temp2[[3]])
    
    # plot
    temp3 = TAN_Bayesian(n=n,alloc.ratio=alloc.ratio,
                         prior.mean1=prior.mean1,prior.sd1=prior.sd1,prior.k1=prior.k1,
                         prior.mean2=prior.mean2,prior.sd2=prior.sd2,prior.k2=prior.k2,
                         prior.a1=prior.a1,prior.b1=prior.b1,
                         prior.a2=prior.a2,prior.b2=prior.b2,
                         control.mean=control.mean,
                         mean=meanseq,
                         sd1=sd1,sd2=sd2,
                         CT1.go=CT1.go,
                         false.go.CT1=false.go.CT1,FGR.CT1=FGR.CT1,
                         CT1.nogo=CT1.nogo,
                         false.nogo.CT1=false.nogo.CT1,FNGR.CT1=FNGR.CT1,
                         CT2.go=CT2.go,
                         false.go.CT2=false.go.CT2, FGR.CT2=FGR.CT2,
                         CT2.nogo=CT2.nogo,
                         false.nogo.CT2=false.nogo.CT2, FNGR.CT2=FNGR.CT2,
                         overlap.option=overlap.option, direction=direction,
                         fix.var=fix.var, noninfo=noninfo,
                         seed.num=seed.num,nsim=nsim,n_repeat=n_repeat,
                         logic.go=logic.go,logic.nogo=logic.nogo)
    
    go_prob_plot = as.numeric(temp3[[1]])
    nogo_prob_plot = as.numeric(temp3[[2]])
    inconclusive_prob_plot = as.numeric(temp3[[3]])
    
  }

###################################################################################
  if(plot.figure==TRUE){

    delta=meanseq
    p_nogo=nogo_prob_plot*100
    p_grey=inconclusive_prob_plot*100
    p_go=go_prob_plot*100

    if(true_go_cutoff!=true_nogo_cutoff){

    par(mfrow=c(1,1),pty="m",bty="l",cex=1.4)


    plot(delta,p_go,xlab=expression(paste("True difference in mean",sep="")),
         ylab="Prob of GO/ NOGO/ Inconclusive (%)",xlim=range(delta),
         ylim=c(0,100),type="n",axes=F)

    axis(1, at=meanseq, labels=meanseq)
    axis(2, at=seq(0,100,10),labels=T)
    #text(rep(-0.6,6),seq(0,100,20),seq(0,100,20),adj=1,xpd=T)
    box()
    #abline(h=seq(0,100,10),lty=3,col=rgb(0.8,0.8,0.8),lwd=1)

    points(delta,p_nogo,type="b", pch=16,col=rgb(0.9,0,0),lwd=3,lty=1)
    points(delta,p_grey,type="b",pch=4,col=rgb(0.9,0.6,0),lwd=3,lty=6)
    points(delta,p_go,type="b",pch=2,col=rgb(0,0.7,0),lwd=3,lty=2)

    text(min(mean),140,bquote(n==~.(n)),xpd=T,adj=0,cex=0.8,col=rgb(0, 0,0))
    if(overlap.flag==0&unsatisfied.flag==0){
      if(direction=='Greater'){
        text(min(mean),130,bquote(GO~symbol("\336")~Observed~difference~'in'~mean~symbol("\263")~.(round(temp$cutoff[1],digits=3))),xpd=T,adj=0,cex=0.8,col=rgb(0, 0.7,0))
        text(min(mean),120,bquote(NOGO~symbol("\336")~Observed~difference~'in'~mean~symbol("\243")~.(round(temp$cutoff[2],digits=3))),xpd=T,adj=0,cex=0.8,col=rgb(0.9, 0,0))
        text(min(mean),110,bquote(Inconclusive~symbol("\336")~anything~between~GO~and~NOGO),xpd=T,adj=0,cex=0.8,col=rgb(0.9,0.6,0))
      }
      if(direction=='Less'){
        text(min(mean),130,bquote(GO~symbol("\336")~Observed~difference~'in'~mean~symbol("\243")~.(round(temp$cutoff[1],digits=3))),xpd=T,adj=0,cex=0.8,col=rgb(0, 0.7,0))
        text(min(mean),120,bquote(NOGO~symbol("\336")~Observed~difference~'in'~mean~symbol("\263")~.(round(temp$cutoff[2],digits=3))),xpd=T,adj=0,cex=0.8,col=rgb(0.9, 0,0))
        text(min(mean),110,bquote(Inconclusive~symbol("\336")~anything~between~GO~and~NOGO),xpd=T,adj=0,cex=0.8,col=rgb(0.9,0.6,0))
      }
    }
    if(overlap.flag==1&unsatisfied.flag==0&overlap.option=='GO'){
      if(direction=="Greater"){
        text(min(mean),130,bquote(GO~symbol("\336")~Observed~difference~'in'~mean~symbol("\263")~.(round(temp$cutoff[1],digits=3))),xpd=T,adj=0,cex=0.8,col=rgb(0, 0.7,0))
        text(min(mean),120,bquote(NOGO~symbol("\336")~Observed~difference~'in'~mean~symbol("\243")~.(round(temp$cutoff[1],digits=3))),xpd=T,adj=0,cex=0.8,col=rgb(0.9, 0,0))
        #text(min(mean),115,paste0('Warning: GO and NOGO zones are overlaped at some sample sizes'),xpd=T,adj=0,cex=0.5,col=rgb(0.9,0.6,0))
        #text(min(mean),110,paste0('The zones are classfied by cutoff of ', overlap.option),xpd=T,adj=0,cex=0.5,col=rgb(0.9,0.6,0))

      }
      if(direction=="Less"){
        text(min(mean),130,bquote(GO~symbol("\336")~Observed~difference~'in'~mean~symbol("\243")~.(round(temp$cutoff[1],digits=3))),xpd=T,adj=0,cex=0.8,col=rgb(0, 0.7,0))
        text(min(mean),120,bquote(NOGO~symbol("\336")~Observed~difference~'in'~mean~symbol("\263")~.(round(temp$cutoff[1],digits=3))),xpd=T,adj=0,cex=0.8,col=rgb(0.9, 0,0))
        #text(min(mean),115,paste0('Warning: GO and NOGO zones are overlaped at some sample sizes'),xpd=T,adj=0,cex=0.5,col=rgb(0.9,0.6,0))
        #text(min(mean),110,paste0('The zones are classfied by cutoff of ', overlap.option),xpd=T,adj=0,cex=0.5,col=rgb(0.9,0.6,0))

      }
    }

    if(overlap.flag==1&unsatisfied.flag==0&overlap.option=='NOGO'){
      if(direction=="Greater"){
        text(min(mean),130,bquote(GO~symbol("\336")~Observed~difference~'in'~mean~symbol("\263")~.(round(temp$cutoff[2],digits=3))),xpd=T,adj=0,cex=0.8,col=rgb(0, 0.7,0))
        text(min(mean),120,bquote(NOGO~symbol("\336")~Observed~difference~'in'~mean~symbol("\243")~.(round(temp$cutoff[2],digits=3))),xpd=T,adj=0,cex=0.8,col=rgb(0.9, 0,0))
        #text(min(mean),115,paste0('Warning: GO and NOGO zones are overlaped at some sample sizes'),xpd=T,adj=0,cex=0.5,col=rgb(0.9,0.6,0))
        #text(min(mean),110,paste0('The zones are classfied by cutoff of ', overlap.option),xpd=T,adj=0,cex=0.5,col=rgb(0.9,0.6,0))
      }
      if(direction=="Less"){
        text(min(mean),130,bquote(GO~symbol("\336")~Observed~difference~'in'~mean~symbol("\243")~.(round(temp$cutoff[2],digits=3))),xpd=T,adj=0,cex=0.8,col=rgb(0, 0.7,0))
        text(min(mean),120,bquote(NOGO~symbol("\336")~Observed~difference~'in'~mean~symbol("\263")~.(round(temp$cutoff[2],digits=3))),xpd=T,adj=0,cex=0.8,col=rgb(0.9, 0,0))
        #text(min(mean),115,paste0('Warning: GO and NOGO zones are overlaped at some sample sizes'),xpd=T,adj=0,cex=0.5,col=rgb(0.9,0.6,0))
        #text(min(mean),110,paste0('The zones are classfied by cutoff of ', overlap.option),xpd=T,adj=0,cex=0.5,col=rgb(0.9,0.6,0))
      }
    }
    if(unsatisfied.flag==1){
      text(min(mean),130,paste0('Warning: We could not find classification of zones to satisfy your both GO and NOGO criterions at some sample sizes'),xpd=T,adj=0,cex=0.5,col=rgb(0.9,0,0))
      text(min(mean),120,paste0('Please check and modify your desicion rule or the range of sample size!'),xpd=T,adj=0,cex=0.5,col=rgb(0.9,0,0))

      }

    }
  }
###################%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  if(plot.figure==TRUE){

    delta=meanseq
    p_nogo=nogo_prob_plot*100
    p_grey=inconclusive_prob_plot*100
    p_go=go_prob_plot*100

    if(true_go_cutoff==true_nogo_cutoff){

      par(mfrow=c(1,1),pty="m",bty="l",cex=1.4)


      plot(delta,p_go,xlab=expression(paste("True difference in mean",sep="")),
           ylab="Prob of GO/ NOGO/ Inconclusive (%)",xlim=range(delta),
           ylim=c(0,100),type="n",axes=F)

      axis(1, at=meanseq, labels=meanseq)
      axis(2, at=seq(0,100,10),labels=T)
      #text(rep(-0.6,6),seq(0,100,20),seq(0,100,20),adj=1,xpd=T)
      box()
      #abline(h=seq(0,100,10),lty=3,col=rgb(0.8,0.8,0.8),lwd=1)

      points(delta,p_nogo,type="b", pch=16,col=rgb(0.9,0,0),lwd=3,lty=1)
      #points(delta,p_grey,type="b",pch=4,col=rgb(0.9,0.6,0),lwd=3,lty=6)
      points(delta,p_go,type="b",pch=2,col=rgb(0,0.7,0),lwd=3,lty=2)

      text(min(mean),140,bquote(n==~.(n)),xpd=T,adj=0,cex=0.8,col=rgb(0, 0,0))
      if(overlap.flag==0&unsatisfied.flag==0){
        if(direction=='Greater'){
          text(min(mean),130,bquote(GO~symbol("\336")~Observed~difference~'in'~mean~symbol("\263")~.(round(temp$cutoff[1],digits=3))),xpd=T,adj=0,cex=0.8,col=rgb(0, 0.7,0))
          text(min(mean),120,bquote(NOGO~symbol("\336")~Observed~difference~'in'~mean~symbol("\243")~.(round(temp$cutoff[2],digits=3))),xpd=T,adj=0,cex=0.8,col=rgb(0.9, 0,0))
          #text(min(mean),110,bquote(Inconclusive~symbol("\336")~anything~between~GO~and~NOGO),xpd=T,adj=0,cex=0.8,col=rgb(0.9,0.6,0))
        }
        if(direction=='Less'){
          text(min(mean),130,bquote(GO~symbol("\336")~Observed~difference~'in'~mean~symbol("\243")~.(round(temp$cutoff[1],digits=3))),xpd=T,adj=0,cex=0.8,col=rgb(0, 0.7,0))
          text(min(mean),120,bquote(NOGO~symbol("\336")~Observed~difference~'in'~mean~symbol("\263")~.(round(temp$cutoff[2],digits=3))),xpd=T,adj=0,cex=0.8,col=rgb(0.9, 0,0))
          #text(min(mean),110,bquote(Inconclusive~symbol("\336")~anything~between~GO~and~NOGO),xpd=T,adj=0,cex=0.8,col=rgb(0.9,0.6,0))
        }
      }
      if(overlap.flag==1&unsatisfied.flag==0&overlap.option=='GO'){
        if(direction=="Greater"){
          text(min(mean),130,bquote(GO~symbol("\336")~Observed~difference~'in'~mean~symbol("\263")~.(round(temp$cutoff[1],digits=3))),xpd=T,adj=0,cex=0.8,col=rgb(0, 0.7,0))
          text(min(mean),120,bquote(NOGO~symbol("\336")~Observed~difference~'in'~mean~symbol("\243")~.(round(temp$cutoff[1],digits=3))),xpd=T,adj=0,cex=0.8,col=rgb(0.9, 0,0))
          #text(min(mean),115,paste0('Warning: GO and NOGO zones are overlaped at some sample sizes'),xpd=T,adj=0,cex=0.5,col=rgb(0.9,0.6,0))
          #text(min(mean),110,paste0('The zones are classfied by cutoff of ', overlap.option),xpd=T,adj=0,cex=0.5,col=rgb(0.9,0.6,0))

        }
        if(direction=="Less"){
          text(min(mean),130,bquote(GO~symbol("\336")~Observed~difference~'in'~mean~symbol("\243")~.(round(temp$cutoff[1],digits=3))),xpd=T,adj=0,cex=0.8,col=rgb(0, 0.7,0))
          text(min(mean),120,bquote(NOGO~symbol("\336")~Observed~difference~'in'~mean~symbol("\263")~.(round(temp$cutoff[1],digits=3))),xpd=T,adj=0,cex=0.8,col=rgb(0.9, 0,0))
          #text(min(mean),115,paste0('Warning: GO and NOGO zones are overlaped at some sample sizes'),xpd=T,adj=0,cex=0.5,col=rgb(0.9,0.6,0))
          #text(min(mean),110,paste0('The zones are classfied by cutoff of ', overlap.option),xpd=T,adj=0,cex=0.5,col=rgb(0.9,0.6,0))

        }
      }

      if(overlap.flag==1&unsatisfied.flag==0&overlap.option=='NOGO'){
        if(direction=="Greater"){
          text(min(mean),130,bquote(GO~symbol("\336")~Observed~difference~'in'~mean~symbol("\263")~.(round(temp$cutoff[2],digits=3))),xpd=T,adj=0,cex=0.8,col=rgb(0, 0.7,0))
          text(min(mean),120,bquote(NOGO~symbol("\336")~Observed~difference~'in'~mean~symbol("\243")~.(round(temp$cutoff[2],digits=3))),xpd=T,adj=0,cex=0.8,col=rgb(0.9, 0,0))
          #text(min(mean),115,paste0('Warning: GO and NOGO zones are overlaped at some sample sizes'),xpd=T,adj=0,cex=0.5,col=rgb(0.9,0.6,0))
          #text(min(mean),110,paste0('The zones are classfied by cutoff of ', overlap.option),xpd=T,adj=0,cex=0.5,col=rgb(0.9,0.6,0))
        }
        if(direction=="Less"){
          text(min(mean),130,bquote(GO~symbol("\336")~Observed~difference~'in'~mean~symbol("\243")~.(round(temp$cutoff[2],digits=3))),xpd=T,adj=0,cex=0.8,col=rgb(0, 0.7,0))
          text(min(mean),120,bquote(NOGO~symbol("\336")~Observed~difference~'in'~mean~symbol("\263")~.(round(temp$cutoff[2],digits=3))),xpd=T,adj=0,cex=0.8,col=rgb(0.9, 0,0))
          #text(min(mean),115,paste0('Warning: GO and NOGO zones are overlaped at some sample sizes'),xpd=T,adj=0,cex=0.5,col=rgb(0.9,0.6,0))
          #text(min(mean),110,paste0('The zones are classfied by cutoff of ', overlap.option),xpd=T,adj=0,cex=0.5,col=rgb(0.9,0.6,0))
        }
      }
      if(unsatisfied.flag==1){
        text(min(mean),130,paste0('Warning: We could not find classification of zones to satisfy your both GO and NOGO criterions at some sample sizes'),xpd=T,adj=0,cex=0.5,col=rgb(0.9,0,0))
        text(min(mean),120,paste0('Please check and modify your desicion rule or the range of sample size!'),xpd=T,adj=0,cex=0.5,col=rgb(0.9,0,0))

      }

    }
  }

################@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

  return(list(go_prob=go_prob,nogo_prob=nogo_prob,
              inconclusive_prob=inconclusive_prob,
              overlap.flag=overlap.flag,overlap.option=overlap.option,
              unsatisfied.flag=unsatisfied.flag,cutoff=temp$cutoff,true_cutoff=c(true_go_cutoff,true_nogo_cutoff))
  )

}



Vary_SS_TAN_Normal_Prob<-function(nmin=250,nmax=300,
                                  prior.mean1=1/3,prior.sd1=1,prior.k1=1,
                                  prior.mean2=1/3,prior.sd2=1,prior.k2=1,
                                  prior.a1=1,prior.b1=1,
                                  prior.a2=1,prior.b2=1,
                                  mean=0.3,alloc.ratio=1,
                                  sd1=1,sd2=1,control.mean=0,
                                  ssq1=1,ssq2=2,
                                  CT1.go=0.25,
                                  false.go.CT1=TRUE,FGR.CT1=0.25,
                                  CT1.nogo=0.25,
                                  false.nogo.CT1=TRUE,FNGR.CT1=0.25,
                                  CT2.go=0.3,
                                  false.go.CT2=TRUE, FGR.CT2=0.5,
                                  CT2.nogo=0.3,
                                  false.nogo.CT2=TRUE, FNGR.CT2=0.5,
                                  overlap.option='GO',plot.cutoff=TRUE,plot.prob=TRUE,
                                  method="Bayesian",direction="Greater",
                                  fix.var=TRUE,noninfo=TRUE,
                                  seed.num=369,
                                  nsim=10000,stop.criterion=10^-3,
                                  logic.go='and',logic.nogo='or'){
  nseq=unique(round(c(seq(nmin,(nmax+nmin)/2,length=6)[-6],seq((nmax+nmin)/2,nmax,length=6))))
  go_prob<-matrix(NA,ncol=length(nseq),nrow=length(mean))
  nogo_prob<-matrix(NA,ncol=length(nseq),nrow=length(mean))
  inconclusive_prob<-matrix(NA,ncol=length(nseq),nrow=length(mean))
  go_cutoff<-rep(NA,length(nseq))
  nogo_cutoff<-rep(NA,length(nseq))
  index=1
  n_unsatisfied=NA
  n_overlap=NA
  ncore <- detectCores()
  cl<-makeCluster(ncore)
  registerDoParallel(cl)
  results<-foreach(i = nseq,.export=c('TAN_Normal_Cutoff','TAN_findIGcutoff'),.packages=c('mvtnorm'),.combine=rbind) %dopar% {


    temp=TAN_Normal_Cutoff(n=i,alloc.ratio=alloc.ratio,
                           prior.mean1=prior.mean1,prior.sd1=prior.sd1,prior.k1=prior.k1,
                           prior.mean2=prior.mean2,prior.sd2=prior.sd2,prior.k2=prior.k2,
                           prior.a1=prior.a1,prior.b1=prior.b1,
                           prior.a2=prior.a2,prior.b2=prior.b2,
                           sd1=sd1,sd2=sd2,control.mean=control.mean,
                           ssq1=ssq1,ssq2=ssq2,
                           CT1.go=CT1.go,
                           false.go.CT1 = false.go.CT1, FGR.CT1=FGR.CT1,
                           CT1.nogo=CT1.nogo,
                           false.nogo.CT1 = false.nogo.CT1, FNGR.CT1 = FNGR.CT1,
                           CT2.go=CT2.go,
                           false.go.CT2 = false.go.CT2, FGR.CT2=FGR.CT2,
                           CT2.nogo=CT2.nogo,
                           false.nogo.CT2 = false.nogo.CT2, FNGR.CT2 = FNGR.CT2,
                           method=method,direction=direction,
                           fix.var=fix.var,noninfo=noninfo,
                           seed.num=seed.num,nsim=nsim,stop.criterion=stop.criterion,
                           logic.go=logic.go,logic.nogo=logic.nogo)
    c(temp$overlap,temp$cutoff)
  }
  
  stopCluster(cl)
  go_cutoff<-results[,2]
  nogo_cutoff<-results[,3]
  overlap<-results[,1]

  for(mean.index in 1:length(mean)){
    true_go_cutoff<-go_cutoff
    true_nogo_cutoff<-nogo_cutoff
    true_go_cutoff[overlap==1]<-(overlap.option=='GO')*go_cutoff[overlap==1]+(overlap.option=='NOGO')*nogo_cutoff[overlap==1]
    true_nogo_cutoff[overlap==1]<-(overlap.option=='GO')*go_cutoff[overlap==1]+(overlap.option=='NOGO')*nogo_cutoff[overlap==1]
  }

  if(method=='Frequentist'){
    
    for(mean.index in 1:length(mean)){
      mean_temp=mean[mean.index]
      i1=round(nseq*alloc.ratio/(1+alloc.ratio))
      i2=nseq-i1
      sd.i=sqrt(sd1^2/i1+sd2^2/i2)
      
      if(direction=='Greater'){
        go_prob[mean.index,]=1-pnorm(true_go_cutoff,mean=mean_temp,sd=sd.i)
        nogo_prob[mean.index,]=pnorm(true_nogo_cutoff,mean=mean_temp,sd=sd.i)
      }
      if(direction=='Less'){
        go_prob[mean.index,]=pnorm(true_go_cutoff,mean=mean_temp,sd=sd.i)
        nogo_prob[mean.index,]=1-pnorm(true_nogo_cutoff,mean=mean_temp,sd=sd.i)
      }
      inconclusive_prob[mean.index,]=1-go_prob[mean.index,]-nogo_prob[mean.index,]
    }
    
  }

  if(method=='Bayesian'){
    
    n_repeat=1000
    
    temp2 = TAN_Bayesian(n=nseq,alloc.ratio=alloc.ratio,
                         prior.mean1=prior.mean1,prior.sd1=prior.sd1,prior.k1=prior.k1,
                         prior.mean2=prior.mean2,prior.sd2=prior.sd2,prior.k2=prior.k2,
                         prior.a1=prior.a1,prior.b1=prior.b1,
                         prior.a2=prior.a2,prior.b2=prior.b2,
                         control.mean=control.mean,
                         mean=mean,
                         sd1=sd1,sd2=sd2,
                         CT1.go=CT1.go,
                         false.go.CT1=false.go.CT1,FGR.CT1=FGR.CT1,
                         CT1.nogo=CT1.nogo,
                         false.nogo.CT1=false.nogo.CT1,FNGR.CT1=FNGR.CT1,
                         CT2.go=CT2.go,
                         false.go.CT2=false.go.CT2, FGR.CT2=FGR.CT2,
                         CT2.nogo=CT2.nogo,
                         false.nogo.CT2=false.nogo.CT2, FNGR.CT2=FNGR.CT2,
                         overlap.option=overlap.option, direction=direction,
                         fix.var=fix.var, noninfo=noninfo,
                         seed.num=seed.num,nsim=nsim,n_repeat=n_repeat,
                         logic.go=logic.go,logic.nogo=logic.nogo)
    
    go_prob = temp2[[1]]
    nogo_prob = temp2[[2]]
    inconclusive_prob = temp2[[3]]
    
  }

  n_overlap=n_overlap[-1]
  n_unsatisfied=n_unsatisfied[-1]
  ####plot figure
  if(plot.prob==TRUE){
    for(j in 1:length(mean)){
      p_go=go_prob[j,]*100
      p_nogo=nogo_prob[j,]*100
      class(p_nogo)
      dim(p_nogo)
      p_grey=100-p_go-p_nogo
      par(mfrow=c(1,1),pty="m",bty="l",cex=1.4)
      cum_p_nogo=p_nogo

      cum_p_grey=p_nogo+p_grey
      cum_p_go=p_nogo+p_grey+p_go
      delta=nseq
      plot(delta,cum_p_go,xlab="Sample size",
           ylab="Prob of GO/ NOGO/ Inconclusive (%)",xlim=c(nmin,nmax),
           ylim=c(0,100),type="n",axes=F,pch=NA)
      axis(1, at=nseq,labels=nseq)
      axis(2, at=seq(0,100,10),labels=T)
      #text(rep(-0.6,6),seq(0,100,20),seq(0,100,20),adj=1,xpd=T)
      box()
      #abline(h=seq(0,100,10),lty=3,col=rgb(0.8,0.8,0.8),lwd=1)
      polygon(c(delta,rev(delta)),c(cum_p_nogo,rev(rep(0,length(delta)))),col=rgb(0.9,0,0),border=NA)
      polygon(c(delta,rev(delta)),c(cum_p_grey,rev(cum_p_nogo)),col=rgb(0.9,0.6,0),border=NA)
      polygon(c(delta,rev(delta)),c(cum_p_go,rev(cum_p_grey)),col=rgb(0,0.7,0),border=NA)
      text((nmin+nmax)/2,120,paste0("True difference in mean=",round(mean[j],3),", SDs in group 1&2=",'(',round(sd1,3),',',round(sd2,3),')'),xpd=T,adj=0.5,cex=0.8)
      if(length(n_overlap)!=0){
        text(nmin,115,paste0('Warning: GO and NOGO cut-offs are overlaped at some sample sizes'),xpd=T,adj=0,cex=0.5,col=rgb(0.9,0.6,0))
        #text(nmin,110,paste(n_overlap),xpd=T,adj=0,cex=0.5,col=rgb(0.9,0.6,0))
        #print(n_overlap)
        text(nmin,110,paste0('The zones are classfied by cutoff of ', overlap.option),xpd=T,adj=0,cex=0.5,col=rgb(0.9,0.6,0))
      }
      if(length(n_unsatisfied)!=0){
        text(nmin,115,paste0('Warning: We could not find cutoffs to satisfy your both GO and NOGO criterions.'),xpd=T,adj=0,cex=0.5,col=rgb(0.9,0,0))
        text(nmin,110,paste0('Please check and modify your desicion rule!'),xpd=T,adj=0,cex=0.5,col=rgb(0.9,0,0))

      }
    }
  }
  if(plot.cutoff==TRUE){
    par(mfrow=c(1,1),pty="m",bty="l",cex=1.4)
    ylim_max=max(c(go_cutoff,nogo_cutoff,round(CT1.go,3),round(CT2.go,3),round(CT1.nogo,3),round(CT2.nogo,3)),na.rm=TRUE)+0.1
    ylim_min=min(c(go_cutoff,nogo_cutoff,round(CT1.go,3),round(CT2.go,3),round(CT1.nogo,3),round(CT2.nogo,3)),na.rm=TRUE)-0.1
    plot(NA,NA,xlab='Sample size',ylab="Observed mean",xlim=c(max(0,range(nseq)[1]-10),range(nseq)[2]+10),ylim=c(ylim_min,ylim_max),type="n",axes=F,col=rgb(1,0,0),lty=1,lwd=2)
    axis(1, at=nseq,labels=nseq)
    axis(2, at=round(c(seq(ylim_min,ylim_max,round((ylim_max-ylim_min)/10,digits=2)),CT1.go,CT2.go,CT1.nogo,CT2.nogo),digits=2),labels=T)
    #text(rep(-0.6,6),seq(0,100,20),seq(0,300,20),adj=1,xpd=T)
    box()
    lines(nseq,nogo_cutoff,col=rgb(0.9,0,0),lwd=2)
    lines(nseq,go_cutoff,col=rgb(0,0.7,0),lwd=2)
    legend('bottomright',legend=c("Cut off of GO","Cut off of NOGO"),
           col=c(rgb(0,0.7,0),rgb(0.9,0,0)),
           lwd=c(2,2),
           lty=c(1,1),cex=0.5)
    if(length(n_overlap)!=0){
      text(nmin,ylim_max+abs(ylim_max-ylim_min)/10,paste0('Warning: GO and NOGO zones are overlaped at some sample sizes'),xpd=T,adj=0,cex=0.5,col=rgb(0.9,0.6,0))
      #text(nmin,110,paste(n_overlap),xpd=T,adj=0,cex=0.5,col=rgb(0.9,0.6,0))
      #print(n_overlap)
      text(nmin,ylim_max+abs(ylim_max-ylim_min)/20,paste0('The zones are classfied by cutoff of ', overlap.option),xpd=T,adj=0,cex=0.5,col=rgb(0.9,0.6,0))
    }
    if(length(n_unsatisfied)!=0){
      text(nmin,ylim_max+abs(ylim_max-ylim_min)/10,paste0('Warning: We could not find classification of zones to satisfy your both GO and NOGO criterions at some sample sizes'),xpd=T,adj=0,cex=0.5,col=rgb(0.9,0,0))
      text(nmin,ylim_max+abs(ylim_max-ylim_min)/20,paste0('Please check and modify your desicion rule or the range of sample size!'),xpd=T,adj=0,cex=0.5,col=rgb(0.9,0,0))

    }


  }
  return(list(n_unsatisfied,n_overlap))


}


Interim_TAN<-function(interim_n=c(50,100,150),num_interim=3,
                      CT1.go=c(0.25,0.25,0.25),
                      false.go.CT1=c(TRUE,TRUE,TRUE),FGR.CT1=c(0.25,0.25,0.25),
                      CT1.nogo=c(0.25,0.25,0.25),
                      false.nogo.CT1=c(TRUE,TRUE,TRUE),FNGR.CT1=c(0.25,0.25,0.25),
                      CT2.go=c(0.25,0.25,0.25),
                      false.go.CT2=c(TRUE,TRUE,TRUE),FGR.CT2=c(0.25,0.25,0.25),
                      CT2.nogo=c(0.25,0.25,0.25),
                      false.nogo.CT2=c(TRUE,TRUE,TRUE),FNGR.CT2=c(0.25,0.25,0.25),
                      overlap.option=c('GO','GO','GO'),
                      method='Bayesian',direction=c("Greater",'Greater','Greater'),
                      task=c('Futility','Superiority','Futility and superiority'),
                      logic.go=c('and','and','or'),
                      logic.nogo=c('or','and','or'),
                      seed.num=369,nsim_IA=10000,stop.criterion=10^-3,
                      alloc.ratio=1,
                      prior.mean1=1/3,prior.sd1=1,prior.k1=1,
                      prior.mean2=1/3,prior.sd2=1,prior.k2=1,
                      prior.a1=1,prior.b1=1,
                      prior.a2=1,prior.b2=1,
                      control.mean=0,
                      sd1=1,sd2=1,
                      ssq1=1,ssq2=2,
                      mean=c(0.25),fix.var=TRUE,noninfo=TRUE,
                      nsim=10000){

  interim_n=sort(interim_n)
  
  temptable=c()
  for(meanindex in 1:length(mean)){

    set.seed(seed.num)
    diff_interim_n<-diff(interim_n)
    generate_n<-c(interim_n[1],diff_interim_n)
    
    mean_x1=matrix(NA,nrow=nsim_IA,ncol=num_interim)
    mean_x2=mean_x1
    if(method=='Bayesian' & fix.var==FALSE){
      ssq_x1=mean_x1
      ssq_x2=mean_x1
    }
    treatment.mean=control.mean+mean[meanindex]
    x1=c()
    x2=c()
    for(k in 1:num_interim){
      i1=round(generate_n[k]*alloc.ratio/(1+alloc.ratio))
      i2=generate_n[k]-i1
      
      if (abs(generate_n[k]*alloc.ratio/(1+alloc.ratio)-i1)>10^(-3))
        stop ("Sample size for each arm must be an integer in each interim.")
      
      x1=cbind(x1,matrix(rnorm(nsim_IA*i1,treatment.mean,sd1),nrow=nsim_IA))
      x2=cbind(x2,matrix(rnorm(nsim_IA*i2,control.mean,sd2),nrow=nsim_IA))
      
      mean_x1[,k]=apply(x1,1,mean)
      mean_x2[,k]=apply(x2,1,mean)
      if(method=='Bayesian' & fix.var==FALSE){
        ssq_x1[,k]=apply(x1,1,function(x){sum(x^2)})
        ssq_x2[,k]=apply(x2,1,function(x){sum(x^2)})
      }
    }


    #####
    go_matrix<-matrix(NA,nrow=nsim_IA,ncol=num_interim)
    nogo_matrix<-matrix(NA,nrow=nsim_IA,ncol=num_interim)
    inconclusive_matrix<-matrix(NA,nrow=nsim_IA,ncol=num_interim)
    table<-matrix(NA,ncol=num_interim+1,nrow=6)
    IA_go_matrix<-matrix(NA,nrow=nsim_IA,ncol=num_interim) ###whether continue to next stage
    
    if(method=='Frequentist'){
      
      go_cutoff<-rep(NA,num_interim)
      nogo_cutoff<-rep(NA,num_interim)
      true_go_cutoff<-rep(NA,num_interim)
      true_nogo_cutoff<-rep(NA,num_interim)
      overlap<-rep(NA,num_interim)
      for(i in 1:num_interim){
        temp=TAN_Normal_Cutoff(n=interim_n[i],alloc.ratio=alloc.ratio,
                               prior.mean1=prior.mean1,prior.sd1=prior.sd1,prior.k1=prior.k1,
                               prior.mean2=prior.mean2,prior.sd2=prior.sd2,prior.k2=prior.k2,
                               prior.a1=prior.a1,prior.b1=prior.b1,
                               prior.a2=prior.a2,prior.b2=prior.b2,
                               sd1=sd1,sd2=sd2,control.mean=control.mean,
                               ssq1=NA,ssq2=NA,
                               CT1.go=CT1.go[i],
                               false.go.CT1=false.go.CT1[i],FGR.CT1=FGR.CT1[i],
                               CT1.nogo=CT1.nogo[i],
                               false.nogo.CT1=false.nogo.CT1[i],FNGR.CT1=FNGR.CT1[i],
                               CT2.go=CT2.go[i],
                               false.go.CT2=false.go.CT2[i], FGR.CT2=FGR.CT2[i],
                               CT2.nogo=CT2.nogo[i],
                               false.nogo.CT2=false.nogo.CT2[i], FNGR.CT2=FNGR.CT2[i],
                               method=method,direction=direction[i],
                               fix.var=fix.var,noninfo=noninfo,
                               seed.num=seed.num,nsim=nsim,stop.criterion=stop.criterion,
                               logic.go=logic.go[i],logic.no=logic.nogo[i])
        go_cutoff[i]<-temp$cutoff[1]
        nogo_cutoff[i]<-temp$cutoff[2]
        overlap[i]<-temp$overlap
        ###TAN####
        true_go_cutoff[i]<-ifelse(overlap[i]==0,go_cutoff[i],(overlap.option[i]=='GO')*go_cutoff[i]+(overlap.option[i]=='NOGO')*(nogo_cutoff[i]))
        true_nogo_cutoff[i]<-ifelse(overlap[i]==0,nogo_cutoff[i],(overlap.option[i]=='GO')*(go_cutoff[i])+(overlap.option[i]=='NOGO')*(nogo_cutoff[i]))
        ####
      }
      
      for(j in 1:num_interim ){
        if(direction[j]=='Greater'){
          go_matrix[,j]<-mean_x1[,k]-mean_x2[,k]>=true_go_cutoff[j]
          nogo_matrix[,j]<-mean_x1[,k]-mean_x2[,k]<true_nogo_cutoff[j]
          inconclusive_matrix[,j]<-rep(1,nsim_IA)-go_matrix[,j]-nogo_matrix[,j]
        }
      }
      
      for(j in 1:num_interim ){
        if(direction[j]=='Less'){
          go_matrix[,j]<-mean_x1[,k]-mean_x2[,k]<=true_go_cutoff[j]
          nogo_matrix[,j]<-mean_x1[,k]-mean_x2[,k]>true_nogo_cutoff[j]
          inconclusive_matrix[,j]<-rep(1,nsim_IA)-go_matrix[,j]-nogo_matrix[,j]
        }
      }
    }
    
    if(method=='Bayesian'){
      
      for(j in 1:num_interim ){
      
        i1 = round(interim_n[j]*alloc.ratio/(1+alloc.ratio))
        i2 = interim_n[j]-i1

        if(is.na(CT1.go[j])) false.go.CT1[j]=FALSE
        if(is.na(CT1.nogo[j])) false.nogo.CT1[j]=FALSE
        if(is.na(CT2.go[j])) false.go.CT2[j]=FALSE
        if(is.na(CT2.nogo[j])) false.nogo.CT2[j]=FALSE
        
        if(false.go.CT1[j]==TRUE) pp_go_1 = rep(0,nsim_IA)
        if(false.nogo.CT1[j]==TRUE) pp_nogo_1 = rep(0,nsim_IA)
        if(false.go.CT2[j]==TRUE) pp_go_2 = rep(0,nsim_IA)
        if(false.nogo.CT2[j]==TRUE) pp_nogo_2 = rep(0,nsim_IA)
        
        for (i in 1:nsim_IA){
          
          if(fix.var==TRUE){
            
            sd.post_1=sqrt(1/((1/prior.sd1)^2+i1/sd1^2))
            sd.post_2=sqrt(1/((1/prior.sd2)^2+i2/sd2^2))
            mupost_1=(1/(1+prior.sd1^2/sd1^2*i1))*prior.mean1+(1/(1+sd1^2/i1/prior.sd1^2))*mean_x1[i,j]
            mupost_2=(1/(1+prior.sd2^2/sd2^2*i2))*prior.mean2+(1/(1+sd2^2/i2/prior.sd2^2))*mean_x2[i,j]
            
            if(direction[j]=='Greater'){ # closed form
              if(false.go.CT1[j]==TRUE) pp_go_1[i] = 1-pnorm(CT1.go[j],mean=mupost_1-mupost_2,sd=(sd.post_1^2+sd.post_2^2)^(1/2))
              if(false.nogo.CT1[j]==TRUE) pp_nogo_1[i] = 1-pnorm(CT1.nogo[j],mean=mupost_1-mupost_2,sd=(sd.post_1^2+sd.post_2^2)^(1/2))
              if(false.go.CT2[j]==TRUE) pp_go_2[i] = 1-pnorm(CT2.go[j],mean=mupost_1-mupost_2,sd=(sd.post_1^2+sd.post_2^2)^(1/2))
              if(false.nogo.CT2[j]==TRUE) pp_nogo_2[i] = 1-pnorm(CT2.nogo[j],mean=mupost_1-mupost_2,sd=(sd.post_1^2+sd.post_2^2)^(1/2))
            }
            if(direction[j]=='Less'){
              if(false.go.CT1[j]==TRUE) pp_go_1[i] = pnorm(CT1.go[j],mean=mupost_1-mupost_2,sd=(sd.post_1^2+sd.post_2^2)^(1/2))
              if(false.nogo.CT1[j]==TRUE) pp_nogo_1[i] = pnorm(CT1.nogo[j],mean=mupost_1-mupost_2,sd=(sd.post_1^2+sd.post_2^2)^(1/2))
              if(false.go.CT2[j]==TRUE) pp_go_2[i] = pnorm(CT2.go[j],mean=mupost_1-mupost_2,sd=(sd.post_1^2+sd.post_2^2)^(1/2))
              if(false.nogo.CT2[j]==TRUE) pp_nogo_2[i] = pnorm(CT2.nogo[j],mean=mupost_1-mupost_2,sd=(sd.post_1^2+sd.post_2^2)^(1/2))
            }
            
          }
          
          if(fix.var==FALSE){
            
            if(noninfo==TRUE){
              prior.mean1=0
              prior.k1=0
              prior.a1=-1/2
              prior.b1=0
              prior.mean2=0
              prior.k2=0
              prior.a2=-1/2
              prior.b2=0
            }
            
            mupost_1=(prior.k1*prior.mean1+i1*mean_x1[i,j])/(prior.k1+i1)
            kpost_1=prior.k1+i1
            apost_1=prior.a1+i1/2
            bpost_1=prior.b1+1/2*ssq_x1[i,j]+prior.k1*i1*(mean_x1[i,j]-prior.mean1)^2/(2*(prior.k1+i1))
            
            mupost_2=(prior.k2*prior.mean2+i2*mean_x2[i,j])/(prior.k2+i2)
            kpost_2=prior.k2+i2
            apost_2=prior.a2+i2/2
            bpost_2=prior.b2+1/2*ssq_x2[i,j]+prior.k2*i2*(mean_x2[i,j]-prior.mean2)^2/(2*(prior.k2+i2))
            
            temp1 = rt(nsim,df=2*apost_1)*sqrt(bpost_1/(apost_1*kpost_1))+mupost_1
            temp2 = rt(nsim,df=2*apost_2)*sqrt(bpost_2/(apost_2*kpost_2))+mupost_2
          
            if(direction[j]=='Greater'){
              if(false.go.CT1[j]==TRUE) pp_go_1[i] = mean(temp1-temp2>=CT1.go[j])
              if(false.nogo.CT1[j]==TRUE) pp_nogo_1[i] = mean(temp1-temp2>=CT1.nogo[j])
              if(false.go.CT2[j]==TRUE) pp_go_2[i] = mean(temp1-temp2>=CT2.go[j])
              if(false.nogo.CT2[j]==TRUE) pp_nogo_2[i] = mean(temp1-temp2>=CT2.nogo[j])
            }
            if(direction[j]=='Less'){
              if(false.go.CT1[j]==TRUE) pp_go_1[i] = mean(temp1-temp2<=CT1.go[j])
              if(false.nogo.CT1[j]==TRUE) pp_nogo_1[i] = mean(temp1-temp2<=CT1.nogo[j])
              if(false.go.CT2[j]==TRUE) pp_go_2[i] = mean(temp1-temp2<=CT2.go[j])
              if(false.nogo.CT2[j]==TRUE) pp_nogo_2[i] = mean(temp1-temp2<=CT2.nogo[j])
            }
          }
          
        }
        
        if(false.go.CT1[j]==TRUE & false.go.CT2[j]==TRUE){
          if(logic.go[j]=='and')
            go_matrix[,j] = pp_go_1>=1-FGR.CT1[j] & pp_go_2>=1-FGR.CT2[j]
          if(logic.go[j]=='or')
            go_matrix[,j] = pp_go_1>=1-FGR.CT1[j] | pp_go_2>=1-FGR.CT2[j]
        } else if(false.go.CT1[j]==TRUE){
          go_matrix[,j] = pp_go_1>=1-FGR.CT1[j]
        } else if(false.go.CT2[j]==TRUE){
          go_matrix[,j] = pp_go_2>=1-FGR.CT2[j]
        }
        
        if(false.nogo.CT1[j]==TRUE & false.nogo.CT2[j]==TRUE){
          if(logic.nogo[j]=='and')
            nogo_matrix[,j] = pp_nogo_1<FNGR.CT1[j] & pp_nogo_2<FNGR.CT2[j]
          if(logic.nogo[j]=='or')
            nogo_matrix[,j] = pp_nogo_1<FNGR.CT1[j] | pp_nogo_2<FNGR.CT2[j]
        } else if(false.nogo.CT1[j]==TRUE){
          nogo_matrix[,j] = pp_nogo_1<FNGR.CT1[j]
        } else if(false.nogo.CT2[j]==TRUE){
          nogo_matrix[,j] = pp_nogo_2<FNGR.CT2[j]
        }
        
        if (overlap.option[j]=='GO'){
          nogo_matrix[nogo_matrix[,j]==go_matrix[,j] & nogo_matrix[,j]==1,j] = 0
        } else if (overlap.option[j]=='NOGO'){
          go_matrix[go_matrix[,j]==nogo_matrix[,j] & go_matrix[,j]==1,j] = 0
        }
        
        inconclusive_matrix[,j] = rep(1,nsim_IA)-go_matrix[,j]-nogo_matrix[,j]
      }
    }

    
    for(ii in 1:(num_interim)){
      if(task[ii]=='Futility'){
        IA_go_matrix[,ii]=inconclusive_matrix[,ii]+go_matrix[,ii]
      }
      if(task[ii]=='Superiority'){
        IA_go_matrix[,ii]=inconclusive_matrix[,ii]+nogo_matrix[,ii]
      }
      if(task[ii]=='Futility and superiority'){
        IA_go_matrix[,ii]=inconclusive_matrix[,ii]
      }
    }
    cum_IA_go_matrix<-t(apply(IA_go_matrix,1,cumprod))
    for(j in 1:(num_interim)){
      table[1,j]=interim_n[j]
      table[2,j]=mean[meanindex]
      table[3,j]=task[j]
      if(j==1){
        if(task[j]=='Superiority'|task[j]=='Futility and superiority'){
          table[4,j]=round(sum(go_matrix[,j]==1)/nsim_IA,3)}else{table[4,j]=0}
        if(task[j]=='Futility'|task[j]=='Futility and superiority'){
          table[6,j]=round(sum(nogo_matrix[,j]==1)/nsim_IA,3)
        }else{table[6,j]=0}
      }else{
        if(task[j]=='Superiority'|task[j]=='Futility and superiority'){
          table[4,j]=round(sum(go_matrix[,j]==1&cum_IA_go_matrix[,j-1]==1)/nsim_IA,3)}else{table[4,j]=0}
        if(task[j]=='Futility'|task[j]=='Futility and superiority'){
          table[6,j]=round(sum(nogo_matrix[,j]==1&cum_IA_go_matrix[,j-1]==1)/nsim_IA,3)
        }else{table[6,j]=0}
      }
      table[5,j]=round(sum(cum_IA_go_matrix[,j]==1)/nsim_IA,3)
      # if(task[j]=='Futility'){
      #   table[7,j]=ifelse(direction[j]=='Greater',HTML(paste0('-/<',round(true_nogo_cutoff[j],3))),HTML(paste0('-/>',round(true_nogo_cutoff[j],3))))
      # }
      # if(task[j]=='Superiority'){
      #   table[7,j]=ifelse(direction[j]=='Greater',HTML(paste0('>=',round(true_go_cutoff[j],3),'/-')),HTML(paste0('<=',round(true_go_cutoff[j],3),'/-')))
      # }
      # if(task[j]=='Futility and superiority'){
      #   table[7,j]=ifelse(direction[j]=='Greater',HTML(paste0('>=',round(true_go_cutoff[j],3),' /','<',round(true_nogo_cutoff[j],3))),HTML(paste0('<=',round(true_go_cutoff[j],3),' / ','>',round(true_nogo_cutoff[j],3))))
      # }
      # table[10,j]=ifelse(direction[j]=='Greater',HTML(paste0('>=',round(go_cutoff[j],3))),HTML(paste0('<=',round(go_cutoff[j],3))))
      # table[9,j]=ifelse(direction[j]=='Greater',HTML(paste0('<',round(nogo_cutoff[j],3))),HTML(paste0('>',round(nogo_cutoff[j],3))))
      # table[8,j]<-ifelse(overlap[j]==1,paste0('GO/NOGO zones overlapped, classified by criterion of ',overlap.option[j]),'None')

    }
    expectss<-round(sum(as.numeric(table[1,1:num_interim])*(c(as.numeric(table[4,1:num_interim-1])+as.numeric(table[6,1:num_interim-1]),as.numeric(table[5,num_interim-1])))),3)
    table[1,num_interim+1]=HTML(paste0(expectss,' (expected)'))
    table[2,num_interim+1]=mean[meanindex]
    table[3,num_interim+1]=''
    table[4,num_interim+1]=round(sum(as.numeric(table[4,1:num_interim])),3)
    table[5,num_interim+1]=round(as.numeric(table[5,num_interim]),3)
    table[6,num_interim+1]=round(sum(as.numeric(table[6,1:num_interim])),3)
    # table[7,num_interim+1]=''
    # table[8,num_interim+1]=''
    # table[9,num_interim+1]=''
    # table[10,num_interim+1]=''
    table<-as.table(table)
    tablecolname<-c(paste0('Interim analysis ',1:(num_interim-1)),'Final analysis',"Summary")
    tablerowname<-c('Sample size','True difference of mean','Task','Success','To next interim/final or inconclusive',
                    'Stop')
    table<-cbind(tablerowname,rep(meanindex,6),table)
    colnames(table)<-c(" ",'Setting',tablecolname)
    temptable=rbind(temptable,table)


  }
  return(temptable)
}
#Interim_TAN(interim_n = c(66,131,197),CT1=c(1,2,1),FGR.CT1=c(1-0.707,0.653,0.0238),FNGR.CT1=c(0.707,1-0.653,1-0.0238),CT2=c(NA,NA,NA,NA),method='Bayesian',direction=c('Greater','Less','Greater'),task=c('Futility','Superiority','Superiority'),mean=c(1,2),fix.var=FALSE,noninfo = FALSE)


#
# Fix_SS_TAN_Normal_Prob(n=90,alloc.ratio=1,
#                        prior.mean1=1/3,prior.sd1=1,prior.k1=1,
#                        prior.mean2=1/3,prior.sd2=1,prior.k2=1,
#                        prior.a1=1,prior.b1=1,
#                        prior.a2=1,prior.b2=1,
#                        control.mean=0,
#                        mean=c(-0.2,0.9),
#                        sd1=1,sd2=1,
#                        ssq1=1,ssq2=2,
#                        CT1.go=0.25,
#                        false.go.CT1=TRUE,FGR.CT1=0.7,
#                        CT1.nogo=0.25,
#                        false.nogo.CT1=TRUE,FNGR.CT1=0.1,
#                        CT2.go=0.3,
#                        false.go.CT2=FALSE, FGR.CT2=0.5,
#                        CT2.nogo=0.3,
#                        false.nogo.CT2=FALSE, FNGR.CT2=0.5,
#                        overlap.option='GO',plot.figure=TRUE,
#                        method='Bayesian',direction='Greater',
#                        fix.var=TRUE,noninfo=TRUE,
#                        seed.num=369,nsim=10000,stop.criterion=10^-3,
#                        logic.go='and',logic.nogo='or')
