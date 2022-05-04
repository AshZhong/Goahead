library(shiny)
SAN_findIGcutoff<-function(n=30,
                           prior.mean=1/3,prior.k=1,
                           prior.a=1,prior.b=1,
                           ssq=1,
                           rate=0.8,CT=0.25,
                           seed.num=369,
                           stop.criterion=10^-3,direction='Greater'){
  ####prior####
  ###mu1~N(prior.mean1,prior.k1*1/var)
  #### 1/var~gamma(prior.a1,rate=prior.b1)
  ####
  cutoff<-NA
  set.seed(seed.num)
  S=1
  upper=abs(CT)+1
  lower=-abs(CT)-1
  ########control.group####
  ssq=ssq*n
  while(S==1){
    temp=(upper+lower)/2
    mupost_temp=(prior.k*prior.mean+n*temp)/(prior.k+n)
    kpost_temp=prior.k+n
    apost_temp=prior.a+n/2
    bpost_temp=prior.b+1/2*ssq+prior.k*n*(temp-prior.mean)^2/(2*(prior.k+n))

    mupost_upper=(prior.k*prior.mean+n*upper)/(prior.k+n)
    kpost_upper=prior.k+n
    apost_upper=prior.a+n/2
    bpost_upper=prior.b+1/2*ssq+prior.k*n*(upper-prior.mean)^2/(2*(prior.k+n))

    mupost_lower=(prior.k*prior.mean+n*lower)/(prior.k+n)
    kpost_lower=prior.k+n
    apost_lower=prior.a+n/2
    bpost_lower=prior.b+1/2*ssq+prior.k*n*(lower-prior.mean)^2/(2*(prior.k+n))
    
    if(direction=='Greater'){
      probupper=1-pt((CT-mupost_upper)/sqrt(bpost_upper/(apost_upper*kpost_upper)),df=2*apost_upper)
      problower=1-pt((CT-mupost_lower)/sqrt(bpost_lower/(apost_lower*kpost_lower)),df=2*apost_lower)
      probtemp=1-pt((CT-mupost_temp)/sqrt(bpost_temp/(apost_temp*kpost_temp)),df=2*apost_temp)
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
      probupper=pt((CT-mupost_upper)/sqrt(bpost_upper/(apost_upper*kpost_upper)),df=2*apost_upper)
      problower=pt((CT-mupost_lower)/sqrt(bpost_lower/(apost_lower*kpost_lower)),df=2*apost_lower)
      probtemp=pt((CT-mupost_temp)/sqrt(bpost_temp/(apost_temp*kpost_temp)),df=2*apost_temp)
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

SAN_Normal_Cutoff<-function(n=30,
                            prior.mean=1/3,prior.sd=1,
                            prior.k=1,prior.a=1,prior.b=1,
                            sd=1,ssq=1,
                            CT1.go=0.25,
                            false.go.CT1=TRUE,FGR.CT1=0.25,
                            CT1.nogo=0.25,
                            false.nogo.CT1=TRUE,FNGR.CT1=0.25,
                            CT2.go=0.3,
                            false.go.CT2=TRUE, FGR.CT2=0.5,
                            CT2.nogo=0.3,
                            false.nogo.CT2=TRUE, FNGR.CT2=0.5,
                            method='Bayesian',direction='Greater',
                            fix.var=TRUE,noninfo=TRUE,
                            seed.num=369,
                            stop.criterion=10^-3,
                            logic.go='and',logic.nogo='or'){
  flag=rep(0,4)
  overlap.flag=0
  est1.go=NA
  est2.go=NA
  est1.nogo=NA
  est2.nogo=NA
  flag=rep(0,4)
  overlap.flag=0
  sd.n=sd/sqrt(n)

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
    if(method=='Bayesian'&(fix.var==TRUE)){
      sd.post=sqrt(1/((1/prior.sd)^2+(1/sd.n)^2))
      if(false.go.CT1==TRUE){
        temp=qnorm(1-FGR.CT1,mean=CT1.go,sd=sd.post)
        est1.go=(temp-(sd.n^2/(prior.sd^2+sd.n^2))*prior.mean)/(prior.sd^2/(prior.sd^2+sd.n^2))
      }
      if(false.go.CT2==TRUE){
        temp=qnorm(1-FGR.CT2,mean=CT2.go,sd=sd.post)
        est2.go=(temp-(sd.n^2/(prior.sd^2+sd.n^2))*prior.mean)/(prior.sd^2/(prior.sd^2+sd.n^2))
      }

      if(false.nogo.CT1==TRUE){
        temp=qnorm(FNGR.CT1,mean=CT1.nogo,sd=sd.post)
        est1.nogo=(temp-(sd.n^2/(prior.sd^2+sd.n^2))*prior.mean)/(prior.sd^2/(prior.sd^2+sd.n^2))

      }
      if(false.nogo.CT2==TRUE){
        temp=qnorm(FNGR.CT2,mean=CT2.nogo,sd=sd.post)
        est2.nogo=(temp-(sd.n^2/(prior.sd^2+sd.n^2))*prior.mean)/(prior.sd^2/(prior.sd^2+sd.n^2))
      }



    }


    if(method=='Bayesian'&(fix.var==FALSE)){
      if(noninfo==TRUE){
        prior.k=0
        prior.a=-1/2
        prior.b=0
        prior.mean=0
      }
      if(false.go.CT1==TRUE){
        est1.go<-SAN_findIGcutoff(n=n,
                                  prior.mean=prior.mean,prior.k=prior.k ,
                                  prior.a=prior.a,prior.b=prior.b,
                                  ssq=ssq,
                                  rate=1-FGR.CT1,CT=CT1.go,
                                  seed.num=seed.num,
                                  stop.criterion=stop.criterion,direction='Greater')
      }
      if(false.go.CT2==TRUE){
        est2.go<-SAN_findIGcutoff(n=n,
                                  prior.mean=prior.mean,prior.k=prior.k ,
                                  prior.a=prior.a,prior.b=prior.b,
                                  ssq=ssq,
                                  rate=1-FGR.CT2,CT=CT2.go,
                                  seed.num=seed.num,
                                  stop.criterion=stop.criterion,direction='Greater')
      }

      if(false.nogo.CT1==TRUE){
        est1.nogo<-SAN_findIGcutoff(n=n,
                                    prior.mean=prior.mean,prior.k=prior.k ,
                                    prior.a=prior.a,prior.b=prior.b,
                                    ssq=ssq,
                                    rate=1-FNGR.CT1,CT=CT1.nogo,
                                    seed.num=seed.num,
                                    stop.criterion=stop.criterion,direction='Less')

      }
      if(false.nogo.CT2==TRUE){
        est2.nogo<-SAN_findIGcutoff(n=n,
                                    prior.mean=prior.mean,prior.k=prior.k ,
                                    prior.a=prior.a,prior.b=prior.b,
                                    ssq=ssq,
                                    rate=1-FNGR.CT2,CT=CT2.nogo,
                                    seed.num=seed.num,
                                    stop.criterion=stop.criterion,direction='Less')
      }



    }


    if(method=='Frequentist'){
      if(false.go.CT1==TRUE){
        est1.go=qnorm(1-FGR.CT1,mean=CT1.go,sd=sd.n)

      }
      if(false.go.CT2==TRUE){
        est2.go=qnorm(1-FGR.CT2,mean=CT2.go,sd=sd.n)

      }

      if(false.nogo.CT1==TRUE){
        est1.nogo=qnorm(FNGR.CT1,mean=CT1.nogo,sd=sd.n)

      }
      if(false.nogo.CT2==TRUE){
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
    if(method=='Bayesian'&(fix.var==TRUE)){
      sd.post=sqrt(1/((1/prior.sd)^2+(1/sd.n)^2))
      if(false.go.CT1==TRUE){
        temp=qnorm(FGR.CT1,mean=CT1.go,sd=sd.post)
        est1.go=(temp-(sd.n^2/(prior.sd^2+sd.n^2))*prior.mean)/(prior.sd^2/(prior.sd^2+sd.n^2))
      }
      if(false.go.CT2==TRUE){
        temp=qnorm(FGR.CT2,mean=CT2.go,sd=sd.post)
        est2.go=(temp-(sd.n^2/(prior.sd^2+sd.n^2))*prior.mean)/(prior.sd^2/(prior.sd^2+sd.n^2))
      }

      if(false.nogo.CT1==TRUE){
        temp=qnorm(1-FNGR.CT1,mean=CT1.nogo,sd=sd.post)
        est1.nogo=(temp-(sd.n^2/(prior.sd^2+sd.n^2))*prior.mean)/(prior.sd^2/(prior.sd^2+sd.n^2))

      }
      if(false.nogo.CT2==TRUE){
        temp=qnorm(1-FNGR.CT2,mean=CT2.nogo,sd=sd.post)
        est2.nogo=(temp-(sd.n^2/(prior.sd^2+sd.n^2))*prior.mean)/(prior.sd^2/(prior.sd^2+sd.n^2))
      }



    }

    if(method=='Bayesian'&(fix.var==FALSE)){
      if(noninfo==TRUE){
        prior.k=0
        prior.a=-1/2
        prior.b=0
        prior.mean=0
      }
      if(false.go.CT1==TRUE){
        est1.go<-SAN_findIGcutoff(n=n,
                                  prior.mean=prior.mean,prior.k=prior.k ,
                                  prior.a=prior.a,prior.b=prior.b,
                                  ssq=ssq,
                                  rate=1-FGR.CT1,CT=CT1.go,
                                  seed.num=seed.num,
                                  stop.criterion=stop.criterion,direction='Less')
      }
      if(false.go.CT2==TRUE){
        est2.go<-SAN_findIGcutoff(n=n,
                                  prior.mean=prior.mean,prior.k=prior.k ,
                                  prior.a=prior.a,prior.b=prior.b,
                                  ssq=ssq,
                                  rate=1-FGR.CT2,CT=CT2.go,
                                  seed.num=seed.num,
                                  stop.criterion=stop.criterion,direction='Less')
      }

      if(false.nogo.CT1==TRUE){
        est1.nogo<-SAN_findIGcutoff(n=n,
                                    prior.mean=prior.mean,prior.k=prior.k ,
                                    prior.a=prior.a,prior.b=prior.b,
                                    ssq=ssq,
                                    rate=1-FNGR.CT1,CT=CT1.nogo,
                                    seed.num=seed.num,
                                    stop.criterion=stop.criterion,direction='Greater')

      }
      if(false.nogo.CT2==TRUE){
        est2.nogo<-SAN_findIGcutoff(n=n,
                                    prior.mean=prior.mean,prior.k=prior.k ,
                                    prior.a=prior.a,prior.b=prior.b,
                                    ssq=ssq,
                                    rate=1-FNGR.CT2,CT=CT2.nogo,
                                    seed.num=seed.num,
                                    stop.criterion=stop.criterion,direction='Greater')
      }



    }

    if(method=='Frequentist'){
      if(false.go.CT1==TRUE){
        est1.go=qnorm(FGR.CT1,mean=CT1.go,sd=sd.n)

      }
      if(false.go.CT2==TRUE){
        est2.go=qnorm(FGR.CT2,mean=CT2.go,sd=sd.n)

      }

      if(false.nogo.CT1==TRUE){
        est1.nogo=qnorm(1-FNGR.CT1,mean=CT1.nogo,sd=sd.n)

      }
      if(false.nogo.CT2==TRUE){
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

SAN_Bayesian_Vary_Var<-function(n=c(10,100),
                                prior.mean=1/3,
                                prior.k=1,prior.a=1,prior.b=1,
                                mean=c(0.3,0.5), sd=1,
                                CT1.go=0.25,
                                false.go.CT1=TRUE,FGR.CT1=0.25,
                                CT1.nogo=0.25,
                                false.nogo.CT1=TRUE,FNGR.CT1=0.25,
                                CT2.go=0.3,
                                false.go.CT2=TRUE, FGR.CT2=0.5,
                                CT2.nogo=0.3,
                                false.nogo.CT2=TRUE, FNGR.CT2=0.5,
                                overlap.option='GO',direction='Greater',
                                noninfo=TRUE,
                                seed.num=369, n_repeat=1000,
                                logic.go='and',logic.nogo='or'){
  set.seed(seed.num)
  
  go_prob<-matrix(NA,ncol=length(n),nrow=length(mean))
  nogo_prob<-matrix(NA,ncol=length(n),nrow=length(mean))
  inconclusive_prob<-matrix(NA,ncol=length(n),nrow=length(mean))
  overlap.flag=rep(0,length(n))
  
  for (ii in 1:length(mean)){
    for (ii2 in 1:length(n)){
      x=matrix(rnorm(n_repeat*n[ii2],mean[ii],sd),nrow=n_repeat)
      mean_x=apply(x,1,mean)
      ssq_x=apply(x,1,function(x){sum(x^2)})
      
      if(is.na(CT1.go)) false.go.CT1=FALSE
      if(is.na(CT1.nogo)) false.nogo.CT1=FALSE
      if(is.na(CT2.go)) false.go.CT2=FALSE
      if(is.na(CT2.nogo)) false.nogo.CT2=FALSE
      
      if(false.go.CT1==TRUE) pp_go_1 = rep(0,n_repeat)
      if(false.nogo.CT1==TRUE) pp_nogo_1 = rep(0,n_repeat)
      if(false.go.CT2==TRUE) pp_go_2 = rep(0,n_repeat)
      if(false.nogo.CT2==TRUE) pp_nogo_2 = rep(0,n_repeat)
      
      for (i in 1:n_repeat){
        
        if(noninfo==TRUE){
          prior.mean=0
          prior.k=0
          prior.a=-1/2
          prior.b=0
        }
        
        mupost=(prior.k*prior.mean+n[ii2]*mean_x[i])/(prior.k+n[ii2])
        kpost=prior.k+n[ii2]
        apost=prior.a+n[ii2]/2
        bpost=prior.b+1/2*ssq_x[i]+prior.k*n[ii2]*(mean_x[i]-prior.mean)^2/(2*(prior.k+n[ii2]))
        
        if(direction=='Greater'){ # closed form. pmvt returns upper area
          if(false.go.CT1==TRUE) pp_go_1[i] = 1-pt((CT1.go-mupost)/sqrt(bpost/(apost*kpost)),df=2*apost)
          if(false.nogo.CT1==TRUE) pp_nogo_1[i] = 1-pt((CT1.nogo-mupost)/sqrt(bpost/(apost*kpost)),df=2*apost)
          if(false.go.CT2==TRUE) pp_go_2[i] = 1-pt((CT2.go-mupost)/sqrt(bpost/(apost*kpost)),df=2*apost)
          if(false.nogo.CT2==TRUE) pp_nogo_2[i] = 1-pt((CT2.nogo-mupost)/sqrt(bpost/(apost*kpost)),df=2*apost)
        }
        if(direction=='Less'){
          if(false.go.CT1==TRUE) pp_go_1[i] = pt((CT1.go-mupost)/sqrt(bpost/(apost*kpost)),df=2*apost)
          if(false.nogo.CT1==TRUE) pp_nogo_1[i] = pt((CT1.nogo-mupost)/sqrt(bpost/(apost*kpost)),df=2*apost)
          if(false.go.CT2==TRUE) pp_go_2[i] = pt((CT2.go-mupost)/sqrt(bpost/(apost*kpost)),df=2*apost)
          if(false.nogo.CT2==TRUE) pp_nogo_2[i] = pt((CT2.nogo-mupost)/sqrt(bpost/(apost*kpost)),df=2*apost)
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

#### this has plot number 3
Fix_SS_SAN_Normal_Prob<-function(n=100,
                                 prior.mean=1/3,prior.sd=1,
                                 prior.k=1,prior.a=1,prior.b=1,
                                 sd=1,ssq=1,
                                 mean=c(0.0,0.9),
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
                                 seed.num=369,
                                 stop.criterion=10^-3,
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
  
  sd.n=sd/sqrt(n)
  temp=SAN_Normal_Cutoff(n=n,
                         prior.mean=prior.mean,prior.sd = prior.sd,
                         prior.k = prior.k,prior.a=prior.a,prior.b=prior.b,
                         sd=sd,ssq=ssq,
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
                         seed.num=seed.num,
                         stop.criterion=stop.criterion,
                         logic.go=logic.go,logic.nogo = logic.nogo)
  ###SAN####
  true_go_cutoff<-ifelse(temp$overlap==0,temp$cutoff[1],(overlap.option=='GO')*temp$cutoff[1]+(overlap.option=='NOGO')*(temp$cutoff[2]))
  true_nogo_cutoff<-ifelse(temp$overlap==0,temp$cutoff[2],(overlap.option=='GO')*temp$cutoff[1]+(overlap.option=='NOGO')*(temp$cutoff[2]))
  
  overlap.flag = temp$overlap
  
  if(method=='Frequentist' | (method=='Bayesian' & fix.var==TRUE) ){
  
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
      }
      if(all(temp$flag==0)&temp$overlap==1&overlap.option=="NOGO"){
        go_prob=1-pnorm(temp$cutoff[2],mean=mean,sd=sd.n)
        nogo_prob=pnorm(temp$cutoff[2],mean=mean,sd=sd.n)
        inconclusive_prob=1-go_prob-nogo_prob
  
        go_prob_plot=1-pnorm(temp$cutoff[2],mean=meanseq,sd=sd.n)
        nogo_prob_plot=pnorm(temp$cutoff[2],mean=meanseq,sd=sd.n)
        inconclusive_prob_plot=1-go_prob_plot-nogo_prob_plot
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
      }
      if(all(temp$flag==0)&temp$overlap==1&overlap.option=="NOGO"){
        go_prob_plot=pnorm(temp$cutoff[2],mean=meanseq,sd=sd.n)
        nogo_prob_plot=1-pnorm(temp$cutoff[2],mean=meanseq,sd=sd.n)
        inconclusive_prob_plot=1-go_prob-nogo_prob_plot
  
        go_prob=pnorm(temp$cutoff[2],mean=mean,sd=sd.n)
        nogo_prob=1-pnorm(temp$cutoff[2],mean=mean,sd=sd.n)
        inconclusive_prob=1-go_prob-nogo_prob
      }
      if(all(temp$flag==0)==FALSE){
        unsatisfied.flag=1
        nogo_prob=rep(NA,length(mean))
        go_prob=rep(NA,length(mean))
        inconclusive_prob=rep(NA,length(mean))
      }
    }
  }
  
  if(method=='Bayesian' & fix.var==FALSE ){
    
    n_repeat=1000
    
    temp2 = SAN_Bayesian_Vary_Var(n=n,
                                  prior.mean=prior.mean,
                                  prior.k=prior.k,prior.a=prior.a,prior.b=prior.b,
                                  sd=sd,mean=mean,
                                  CT1.go=CT1.go,false.go.CT1=false.go.CT1,
                                  FGR.CT1=FGR.CT1,
                                  CT1.nogo=CT1.nogo,
                                  false.nogo.CT1=false.nogo.CT1,FNGR.CT1=FNGR.CT1,
                                  CT2.go=CT2.go,false.go.CT2=false.go.CT2, FGR.CT2=FGR.CT2,
                                  CT2.nogo=CT2.nogo,false.nogo.CT2=false.nogo.CT2,FNGR.CT2=FNGR.CT2,
                                  overlap.option=overlap.option,direction=direction,
                                  noninfo=noninfo,
                                  seed.num=seed.num,n_repeat=n_repeat,
                                  logic.go=logic.go,logic.nogo=logic.nogo)
    
    go_prob = as.numeric(temp2[[1]])
    nogo_prob = as.numeric(temp2[[2]])
    inconclusive_prob = as.numeric(temp2[[3]])
    
    # plot
    temp3 = SAN_Bayesian_Vary_Var(n=n,
                                  prior.mean=prior.mean,
                                  prior.k=prior.k,prior.a=prior.a,prior.b=prior.b,
                                  sd=sd,mean=meanseq,
                                  CT1.go=CT1.go,false.go.CT1=false.go.CT1,
                                  FGR.CT1=FGR.CT1,
                                  CT1.nogo=CT1.nogo,
                                  false.nogo.CT1=false.nogo.CT1,FNGR.CT1=FNGR.CT1,
                                  CT2.go=CT2.go,false.go.CT2=false.go.CT2, FGR.CT2=FGR.CT2,
                                  CT2.nogo=CT2.nogo,false.nogo.CT2=false.nogo.CT2,FNGR.CT2=FNGR.CT2,
                                  overlap.option=overlap.option,direction=direction,
                                  noninfo=noninfo,
                                  seed.num=seed.num,n_repeat=n_repeat,
                                  logic.go=logic.go,logic.nogo=logic.nogo)
    
    go_prob_plot = as.numeric(temp3[[1]])
    nogo_prob_plot = as.numeric(temp3[[2]])
    inconclusive_prob_plot = as.numeric(temp3[[3]])
    
  }
  

########################################################################
  if(plot.figure==TRUE){

    delta=meanseq
    p_nogo=nogo_prob_plot*100
    p_grey=inconclusive_prob_plot*100
    p_go=go_prob_plot*100

    if(true_go_cutoff!=true_nogo_cutoff){

    par(mfrow=c(1,1),pty="m",bty="l",cex=1.4)

    plot(delta,p_go,xlab=expression(paste("True ",mean,sep="")),
         ylab="Prob of GO/ NOGO/ Inconclusive (%)",xlim=range(delta),
         ylim=c(0,100),type="n",axes=F)

    axis(1, at=meanseq, labels=meanseq)
    axis(2, at=seq(0,100,10),labels=T)
    #text(rep(-0.6,6),seq(0,100,20),seq(0,100,20),adj=1,xpd=T)
    box()
    #abline(h=seq(0,100,10),lty=3,col=rgb(0.8,0.8,0.8),lwd=1)

    points(delta,p_nogo,type="b",pch=16,col=rgb(0.9,0,0),lwd=3,lty=1)
    points(delta,p_grey,type="b",pch=4,col=rgb(0.9,0.6,0),lwd=3,lty=6)
    points(delta,p_go,type="b",pch=2,col=rgb(0,0.7,0),lwd=3,lty=2)

    text(min(mean),140,bquote(n==~.(n)),xpd=T,adj=0,cex=0.8,col=rgb(0, 0,0))
    if(overlap.flag==0&unsatisfied.flag==0){
      if(direction=='Greater'){
        text(min(mean),130,bquote(GO~symbol("\336")~Observed~average~effect~symbol("\263")~.(round(temp$cutoff[1],digits=3))),xpd=T,adj=0,cex=0.8,col=rgb(0, 0.7,0))
        text(min(mean),120,bquote(NOGO~symbol("\336")~Observed~average~effect~symbol("\243")~.(round(temp$cutoff[2],digits=3))),xpd=T,adj=0,cex=0.8,col=rgb(0.9, 0,0))
        text(min(mean),110,bquote(Inconclusive~symbol("\336")~anything~between~GO~and~NOGO),xpd=T,adj=0,cex=0.8,col=rgb(0.9,0.6,0))
      }
      if(direction=='Less'){
        text(min(mean),130,bquote(GO~symbol("\336")~Observed~average~effect~symbol("\243")~.(round(temp$cutoff[1],digits=3))),xpd=T,adj=0,cex=0.8,col=rgb(0, 0.7,0))
        text(min(mean),120,bquote(NOGO~symbol("\336")~Observed~average~effect~symbol("\263")~.(round(temp$cutoff[2],digits=3))),xpd=T,adj=0,cex=0.8,col=rgb(0.9, 0,0))
        text(min(mean),110,bquote(Inconclusive~symbol("\336")~anything~between~GO~and~NOGO),xpd=T,adj=0,cex=0.8,col=rgb(0.9,0.6,0))
      }
    }
    if(overlap.flag==1&unsatisfied.flag==0&overlap.option=='GO'){
      if(direction=="Greater"){
        text(min(mean),130,bquote(GO~symbol("\336")~Observed~average~effect~symbol("\263")~.(round(temp$cutoff[1],digits=3))),xpd=T,adj=0,cex=0.8,col=rgb(0, 0.7,0))
        text(min(mean),120,bquote(NOGO~symbol("\336")~Observed~average~effect~symbol("\243")~.(round(temp$cutoff[1],digits=3))),xpd=T,adj=0,cex=0.8,col=rgb(0.9, 0,0))
        #text(min(mean),115,paste0('Warning: GO and NOGO zones are overlaped at some sample sizes'),xpd=T,adj=0,cex=0.5,col=rgb(0.9,0.6,0))
        #text(min(mean),110,paste0('The zones are classfied by cutoff of ', overlap.option),xpd=T,adj=0,cex=0.5,col=rgb(0.9,0.6,0))

      }
      if(direction=="Less"){
        text(min(mean),130,bquote(GO~symbol("\336")~Observed~average~effect~symbol("\243")~.(round(temp$cutoff[1],digits=3))),xpd=T,adj=0,cex=0.8,col=rgb(0, 0.7,0))
        text(min(mean),120,bquote(NOGO~symbol("\336")~Observed~average~effect~symbol("\263")~.(round(temp$cutoff[1],digits=3))),xpd=T,adj=0,cex=0.8,col=rgb(0.9, 0,0))
        #text(min(mean),115,paste0('Warning: GO and NOGO zones are overlaped at some sample sizes'),xpd=T,adj=0,cex=0.5,col=rgb(0.9,0.6,0))
        #text(min(mean),110,paste0('The zones are classfied by cutoff of ', overlap.option),xpd=T,adj=0,cex=0.5,col=rgb(0.9,0.6,0))

      }
    }

    if(overlap.flag==1&unsatisfied.flag==0&overlap.option=='NOGO'){
      if(direction=="Greater"){
        text(min(mean),130,bquote(GO~symbol("\336")~Observed~average~effect~symbol("\263")~.(round(temp$cutoff[2],digits=3))),xpd=T,adj=0,cex=0.8,col=rgb(0, 0.7,0))
        text(min(mean),120,bquote(NOGO~symbol("\336")~Observed~average~effect~symbol("\243")~.(round(temp$cutoff[2],digits=3))),xpd=T,adj=0,cex=0.8,col=rgb(0.9, 0,0))
        #text(min(mean),115,paste0('Warning: GO and NOGO zones are overlaped at some sample sizes'),xpd=T,adj=0,cex=0.5,col=rgb(0.9,0.6,0))
        #text(min(mean),110,paste0('The zones are classfied by cutoff of ', overlap.option),xpd=T,adj=0,cex=0.5,col=rgb(0.9,0.6,0))
      }
      if(direction=="Less"){
        text(min(mean),130,bquote(GO~symbol("\336")~Observed~average~effect~symbol("\243")~.(round(temp$cutoff[2],digits=3))),xpd=T,adj=0,cex=0.8,col=rgb(0, 0.7,0))
        text(min(mean),120,bquote(NOGO~symbol("\336")~Observed~average~effect~symbol("\263")~.(round(temp$cutoff[2],digits=3))),xpd=T,adj=0,cex=0.8,col=rgb(0.9, 0,0))
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

###################################################################

  if(plot.figure==TRUE){

    delta=meanseq
    p_nogo=nogo_prob_plot*100
    p_grey=inconclusive_prob_plot*100
    p_go=go_prob_plot*100

    if(true_go_cutoff==true_nogo_cutoff){

    par(mfrow=c(1,1),pty="m",bty="l",cex=1.4)

    plot(delta,p_go,xlab=expression(paste("True ",mean,sep="")),
         ylab="Prob of GO/ NOGO/ Inconclusive (%)",xlim=range(delta),
         ylim=c(0,100),type="n",axes=F)

    axis(1, at=meanseq, labels=meanseq)
    axis(2, at=seq(0,100,10),labels=T)
    #text(rep(-0.6,6),seq(0,100,20),seq(0,100,20),adj=1,xpd=T)
    box()

    points(delta,p_nogo,type="b",pch=16,col=rgb(0.9,0,0),lwd=3,lty=1)
    #points(delta,p_grey,type="b",pch=4,col=rgb(0.9,0.6,0),lwd=3,lty=6)
    points(delta,p_go,type="b",pch=2,col=rgb(0,0.7,0),lwd=3,lty=2)

    text(min(mean),140,bquote(n==~.(n)),xpd=T,adj=0,cex=0.8,col=rgb(0, 0,0))
    if(overlap.flag==0&unsatisfied.flag==0){
      if(direction=='Greater'){
        text(min(mean),130,bquote(GO~symbol("\336")~Observed~average~effect~symbol("\263")~.(round(temp$cutoff[1],digits=3))),xpd=T,adj=0,cex=0.8,col=rgb(0, 0.7,0))
        text(min(mean),120,bquote(NOGO~symbol("\336")~Observed~average~effect~symbol("\243")~.(round(temp$cutoff[2],digits=3))),xpd=T,adj=0,cex=0.8,col=rgb(0.9, 0,0))
        #text(min(mean),110,bquote(Inconclusive~symbol("\336")~anything~between~GO~and~NOGO),xpd=T,adj=0,cex=0.8,col=rgb(0.9,0.6,0))
      }
      if(direction=='Less'){
        text(min(mean),130,bquote(GO~symbol("\336")~Observed~average~effect~symbol("\243")~.(round(temp$cutoff[1],digits=3))),xpd=T,adj=0,cex=0.8,col=rgb(0, 0.7,0))
        text(min(mean),120,bquote(NOGO~symbol("\336")~Observed~average~effect~symbol("\263")~.(round(temp$cutoff[2],digits=3))),xpd=T,adj=0,cex=0.8,col=rgb(0.9, 0,0))
        #text(min(mean),110,bquote(Inconclusive~symbol("\336")~anything~between~GO~and~NOGO),xpd=T,adj=0,cex=0.8,col=rgb(0.9,0.6,0))
      }
    }
    if(overlap.flag==1&unsatisfied.flag==0&overlap.option=='GO'){
      if(direction=="Greater"){
        text(min(mean),130,bquote(GO~symbol("\336")~Observed~average~effect~symbol("\263")~.(round(temp$cutoff[1],digits=3))),xpd=T,adj=0,cex=0.8,col=rgb(0, 0.7,0))
        text(min(mean),120,bquote(NOGO~symbol("\336")~Observed~average~effect~symbol("\243")~.(round(temp$cutoff[1],digits=3))),xpd=T,adj=0,cex=0.8,col=rgb(0.9, 0,0))
        #text(min(mean),115,paste0('Warning: GO and NOGO zones are overlaped at some sample sizes'),xpd=T,adj=0,cex=0.5,col=rgb(0.9,0.6,0))
        #text(min(mean),110,paste0('The zones are classfied by cutoff of ', overlap.option),xpd=T,adj=0,cex=0.5,col=rgb(0.9,0.6,0))

      }
      if(direction=="Less"){
        text(min(mean),130,bquote(GO~symbol("\336")~Observed~average~effect~symbol("\243")~.(round(temp$cutoff[1],digits=3))),xpd=T,adj=0,cex=0.8,col=rgb(0, 0.7,0))
        text(min(mean),120,bquote(NOGO~symbol("\336")~Observed~average~effect~symbol("\263")~.(round(temp$cutoff[1],digits=3))),xpd=T,adj=0,cex=0.8,col=rgb(0.9, 0,0))
        #text(min(mean),115,paste0('Warning: GO and NOGO zones are overlaped at some sample sizes'),xpd=T,adj=0,cex=0.5,col=rgb(0.9,0.6,0))
        #text(min(mean),110,paste0('The zones are classfied by cutoff of ', overlap.option),xpd=T,adj=0,cex=0.5,col=rgb(0.9,0.6,0))

      }
    }

    if(overlap.flag==1&unsatisfied.flag==0&overlap.option=='NOGO'){
      if(direction=="Greater"){
        text(min(mean),130,bquote(GO~symbol("\336")~Observed~average~effect~symbol("\263")~.(round(temp$cutoff[2],digits=3))),xpd=T,adj=0,cex=0.8,col=rgb(0, 0.7,0))
        text(min(mean),120,bquote(NOGO~symbol("\336")~Observed~average~effect~symbol("\243")~.(round(temp$cutoff[2],digits=3))),xpd=T,adj=0,cex=0.8,col=rgb(0.9, 0,0))
        #text(min(mean),115,paste0('Warning: GO and NOGO zones are overlaped at some sample sizes'),xpd=T,adj=0,cex=0.5,col=rgb(0.9,0.6,0))
        #text(min(mean),110,paste0('The zones are classfied by cutoff of ', overlap.option),xpd=T,adj=0,cex=0.5,col=rgb(0.9,0.6,0))
      }
      if(direction=="Less"){
        text(min(mean),130,bquote(GO~symbol("\336")~Observed~average~effect~symbol("\243")~.(round(temp$cutoff[2],digits=3))),xpd=T,adj=0,cex=0.8,col=rgb(0, 0.7,0))
        text(min(mean),120,bquote(NOGO~symbol("\336")~Observed~average~effect~symbol("\263")~.(round(temp$cutoff[2],digits=3))),xpd=T,adj=0,cex=0.8,col=rgb(0.9, 0,0))
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


###################################################################

  return(list(go_prob=go_prob,nogo_prob=nogo_prob,
              inconclusive_prob=inconclusive_prob,
              overlap.flag=overlap.flag,overlap.option=overlap.option,
              unsatisfied.flag=unsatisfied.flag,cutoff=temp$cutoff,true_cutoff=c(true_go_cutoff,true_nogo_cutoff))
  )

}




Vary_SS_SAN_Normal_Prob<-function(nmin=250,nmax=300,prior.mean=1/3,prior.sd=1,
                                  prior.k=1,prior.a=1,prior.b=1,
                                  mean=0.3,
                                  sd=1,ssq=1,
                                  CT1.go=0.25,
                                  false.go.CT1=TRUE,FGR.CT1=0.25,
                                  CT1.nogo=0.25,
                                  false.nogo.CT1=TRUE,FNGR.CT1=0.25,
                                  CT2.go=0.3,
                                  false.go.CT2=TRUE, FGR.CT2=0.5,
                                  CT2.nogo=0.3,
                                  false.nogo.CT2=TRUE, FNGR.CT2=0.5,
                                  overlap.option='GO',plot.cutoff=TRUE,plot.prob=TRUE,
                                  method="Bayesian",direction="Greater",fix.var=TRUE,noninfo=TRUE,
                                  seed.num=369,
                                  stop.criterion=10^-3,
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
  results<-foreach(i = nseq,.export=c('SAN_Normal_Cutoff','SAN_findIGcutoff'),.packages=c('mvtnorm'),.combine=rbind) %dopar% {


    temp=SAN_Normal_Cutoff(n=i,
                           prior.mean=prior.mean,prior.sd=prior.sd,
                           prior.k = prior.k,prior.a=prior.a,prior.b=prior.b,
                           sd=sd,ssq=ssq,
                           CT1.go=CT1.go,
                           false.go.CT1 = false.go.CT1, FGR.CT1=FGR.CT1,
                           CT1.nogo=CT1.nogo,
                           false.nogo.CT1 = false.nogo.CT1, FNGR.CT1 = FNGR.CT1,
                           CT2.go=CT2.go,
                           false.go.CT2 = false.go.CT2, FGR.CT2=FGR.CT2,
                           CT2.nogo=CT2.nogo,
                           false.nogo.CT2 = false.nogo.CT2, FNGR.CT2 = FNGR.CT2,
                           method=method,direction=direction,
                           fix.var = fix.var,noninfo=noninfo,
                           seed.num = seed.num,stop.criterion=stop.criterion,
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
  
  if(method=='Frequentist' | (method=='Bayesian' & fix.var==TRUE) ){
    
    for(mean.index in 1:length(mean)){
      mean_temp=mean[mean.index]
      sd.i=sqrt(sd^2/nseq)
      
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
  
  if(method=='Bayesian' & fix.var==FALSE ){
    
    n_repeat=1000
    
    temp2 = SAN_Bayesian_Vary_Var(n=nseq,
                                  prior.mean=prior.mean,
                                  prior.k=prior.k,prior.a=prior.a,prior.b=prior.b,
                                  mean=mean, sd=sd,
                                  CT1.go=CT1.go,
                                  false.go.CT1=false.go.CT1,FGR.CT1=FGR.CT1,
                                  CT1.nogo=CT1.nogo,
                                  false.nogo.CT1=false.nogo.CT1,FNGR.CT1=FNGR.CT1,
                                  CT2.go=CT2.go,
                                  false.go.CT2=false.go.CT2, FGR.CT2=FGR.CT2,
                                  CT2.nogo=CT2.nogo,
                                  false.nogo.CT2=false.nogo.CT2, FNGR.CT2=FNGR.CT2,
                                  overlap.option=overlap.option,direction=direction,
                                  noninfo=noninfo,
                                  seed.num=seed.num, n_repeat=n_repeat,
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
      text((nmin+nmax)/2,120,paste0("True mean=",round(mean[j],3),", sd=",round(sd,3)),xpd=T,adj=0.5,cex=0.8)
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

Interim_SAN<-function(interim_n=c(50,100,150),num_interim=3,
                      CT1.go=c(0.25,0.25,0.25),
                      false.go.CT1=c(TRUE,TRUE,TRUE),FGR.CT1=c(0.25,0.25,0.25),
                      CT1.nogo=c(0.25,0.25,0.25),
                      false.nogo.CT1=c(TRUE,TRUE,TRUE),FNGR.CT1=c(0.25,0.25,0.25),
                      CT2.go=c(0.25,0.25,0.25),
                      false.go.CT2=c(TRUE,TRUE,TRUE),FGR.CT2=c(0.25,0.25,0.25),
                      CT2.nogo=c(0.25,0.25,0.25),
                      false.nogo.CT2=c(TRUE,TRUE,TRUE),FNGR.CT2=c(0.25,0.25,0.25),
                      overlap.option=c('GO','GO','GO'),
                      method='Bayesian',direction="Greater",
                      task=c('Futility','Superiority','Futility and superiority'),
                      seed.num=369,nsim_IA=10000,stop.criterion=10^-3,
                      prior.mean = 1,prior.sd=1,
                      prior.k=1,prior.a=1,
                      prior.b=1,sd=1,ssq=1,
                      mean=c(0.25),fix.var=TRUE,noninfo=TRUE,
                      logic.go=c('and','or','and'),logic.nogo=c('or','and','and')){

  interim_n=sort(interim_n)

  ###SAN####
  temptable=c()
  for(meanindex in 1:length(mean)){
    set.seed(seed.num)
    
    diff_interim_n<-diff(interim_n)
    generate_n<-c(interim_n[1],diff_interim_n)
    
    mean_x=matrix(NA,nrow=nsim_IA,ncol=num_interim)
    if(method=='Bayesian' & fix.var==FALSE){
      ssq_x=mean_x
    }
    x=c()
    for(k in 1:num_interim){
      
      x=cbind(x,matrix(rnorm(nsim_IA*generate_n[k],mean[meanindex],sd),nrow=nsim_IA))
      
      mean_x[,k]=apply(x,1,mean)
      if(method=='Bayesian' & fix.var==FALSE){
        ssq_x[,k]=apply(x,1,function(x){sum(x^2)})
      }
    }
    
    #####
    go_matrix<-matrix(NA,nrow=nsim_IA,ncol=num_interim)
    nogo_matrix<-matrix(NA,nrow=nsim_IA,ncol=num_interim)
    inconclusive_matrix<-matrix(NA,nrow=nsim_IA,ncol=num_interim)
    table<-matrix(NA,ncol=num_interim+1,nrow=6)
    IA_go_matrix<-matrix(NA,nrow=nsim_IA,ncol=num_interim) ###whether continue to next stage
    
    if(method=='Frequentist' | (method=='Bayesian' & fix.var==TRUE) ){
      
      go_cutoff<-rep(NA,num_interim)
      nogo_cutoff<-rep(NA,num_interim)
      true_go_cutoff<-rep(NA,num_interim)
      true_nogo_cutoff<-rep(NA,num_interim)
      overlap<-rep(NA,num_interim)
      for(j in 1:num_interim ){
        temp<-SAN_Normal_Cutoff(n=interim_n[j],
                                prior.mean = prior.mean,prior.sd = prior.sd,
                                prior.k=prior.k,prior.a=prior.a,
                                prior.b=prior.b,sd=sd,ssq=NA, #ssq is unused
                                CT1.go=CT1.go[j],
                                false.go.CT1=false.go.CT1[j],FGR.CT1=FGR.CT1[j],
                                CT1.nogo=CT1.nogo[j],
                                false.nogo.CT1=false.nogo.CT1[j],FNGR.CT1=FNGR.CT1[j],
                                CT2.go=CT2.go[j],
                                false.go.CT2=false.go.CT2[j], FGR.CT2=FGR.CT2[j],
                                CT2.nogo=CT2.nogo[j],
                                false.nogo.CT2=false.nogo.CT2[j], FNGR.CT2=FNGR.CT2[j],
                                method=method,direction=direction[j],
                                fix.var=fix.var,noninfo=noninfo,
                                seed.num = seed.num,stop.criterion = stop.criterion,
                                logic.go=logic.go[j],logic.nogo=logic.nogo[j])
        go_cutoff[j]<-temp$cutoff[1]
        nogo_cutoff[j]<-temp$cutoff[2]
        overlap[j]<-temp$overlap
        ###SAN####
        true_go_cutoff[j]<-ifelse(overlap[j]==0,go_cutoff[j],(overlap.option[j]=='GO')*go_cutoff[j]+(overlap.option[j]=='NOGO')*(nogo_cutoff[j]))
        true_nogo_cutoff[j]<-ifelse(overlap[j]==0,nogo_cutoff[j],(overlap.option[j]=='GO')*(go_cutoff[j])+(overlap.option[j]=='NOGO')*(nogo_cutoff[j]))
        ####
      }
      
      for(j in 1:num_interim ){
        if(direction[j]=='Greater'){
          go_matrix[,j]<-mean_x[,j]>=true_go_cutoff[j]
          nogo_matrix[,j]<-mean_x[,j]<true_nogo_cutoff[j]
          inconclusive_matrix[,j]<-rep(1,nsim_IA)-go_matrix[,j]-nogo_matrix[,j]
        }
      
        if(direction[j]=='Less'){
          go_matrix[,j]<-mean_x[,j]<=true_go_cutoff[j]
          nogo_matrix[,j]<-mean_x[,j]>true_nogo_cutoff[j]
          inconclusive_matrix[,j]<-rep(1,nsim_IA)-go_matrix[,j]-nogo_matrix[,j]
        }
      }
    }
    
    if(method=='Bayesian' & fix.var==FALSE){
      
      for(j in 1:num_interim ){
        
        if(is.na(CT1.go[j])) false.go.CT1[j]=FALSE
        if(is.na(CT1.nogo[j])) false.nogo.CT1[j]=FALSE
        if(is.na(CT2.go[j])) false.go.CT2[j]=FALSE
        if(is.na(CT2.nogo[j])) false.nogo.CT2[j]=FALSE
        
        if(false.go.CT1[j]==TRUE) pp_go_1 = rep(0,nsim_IA)
        if(false.nogo.CT1[j]==TRUE) pp_nogo_1 = rep(0,nsim_IA)
        if(false.go.CT2[j]==TRUE) pp_go_2 = rep(0,nsim_IA)
        if(false.nogo.CT2[j]==TRUE) pp_nogo_2 = rep(0,nsim_IA)
        
        for (i in 1:nsim_IA){
            
          if(noninfo==TRUE){
            prior.mean=0
            prior.k=0
            prior.a=-1/2
            prior.b=0
          }
          
          mupost=(prior.k*prior.mean+interim_n[j]*mean_x[i,j])/(prior.k+interim_n[j])
          kpost=prior.k+interim_n[j]
          apost=prior.a+interim_n[j]/2
          bpost=prior.b+1/2*ssq_x[i,j]+prior.k*interim_n[j]*(mean_x[i,j]-prior.mean)^2/(2*(prior.k+interim_n[j]))
          
          if(direction[j]=='Greater'){ # closed form. pmvt returns upper area
            if(false.go.CT1[j]==TRUE) pp_go_1[i] = 1-pt((CT1.go[j]-mupost)/sqrt(bpost/(apost*kpost)),df=2*apost)
            if(false.nogo.CT1[j]==TRUE) pp_nogo_1[i] = 1-pt((CT1.nogo[j]-mupost)/sqrt(bpost/(apost*kpost)),df=2*apost)
            if(false.go.CT2[j]==TRUE) pp_go_2[i] = 1-pt((CT2.go[j]-mupost)/sqrt(bpost/(apost*kpost)),df=2*apost)
            if(false.nogo.CT2[j]==TRUE) pp_nogo_2[i] = 1-pt((CT2.nogo[j]-mupost)/sqrt(bpost/(apost*kpost)),df=2*apost)
          }
          if(direction[j]=='Less'){
            if(false.go.CT1[j]==TRUE) pp_go_1[i] = pt((CT1.go[j]-mupost)/sqrt(bpost/(apost*kpost)),df=2*apost)
            if(false.nogo.CT1[j]==TRUE) pp_nogo_1[i] = pt((CT1.nogo[j]-mupost)/sqrt(bpost/(apost*kpost)),df=2*apost)
            if(false.go.CT2[j]==TRUE) pp_go_2[i] = pt((CT2.go[j]-mupost)/sqrt(bpost/(apost*kpost)),df=2*apost)
            if(false.nogo.CT2[j]==TRUE) pp_nogo_2[i] = pt((CT2.nogo[j]-mupost)/sqrt(bpost/(apost*kpost)),df=2*apost)
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
    tablerowname<-c('Sample size','True mean','Task','Success','To next interim/final or inconclusive',
                    'Stop')
    table<-cbind(tablerowname,rep(meanindex,6),table)
    colnames(table)<-c(" ",'Setting',tablecolname)
    temptable=rbind(temptable,table)


  }
  return(temptable)
}

#Interim_SAN(interim_n = c(66,131,197),CT1.go=c(1,0.6,1),CT1.nogo=c(1,0.6,1),FGR.CT1=c(1-0.707,0.653,0.0238),FNGR.CT1=c(0.707,1-0.653,1-0.0238),CT2.go=c(NA,NA,NA),CT2.nogo=c(NA,NA,NA),method='Frequentist',direction=c(rep('Less',1),'Greater','Less'),task=c('Futility','Superiority','Superiority'),mean=c(1.1,2))

# Fix_SS_SAN_Normal_Prob(n=50,
#                                  prior.mean=1,prior.sd=1000,
#                                  prior.k=1,prior.a=1,prior.b=1,
#                                  sd=1,ssq=1,
#                                  mean=c(2,5),
#                                  CT1.go=4.8,
#                                  false.go.CT1=TRUE,FGR.CT1=0.8,
#                                  CT1.nogo=4.8,
#                                  false.nogo.CT1=TRUE,FNGR.CT1=0.8,
#                                  CT2.go=0.3,
#                                  false.go.CT2=FALSE, FGR.CT2=0.5,
#                                  CT2.nogo=0.3,
#                                  false.nogo.CT2=FALSE, FNGR.CT2=0.5,
#                                  overlap.option='GO',plot.figure=TRUE,
#                                  method='Frequentist',direction='Greater',
#                                  fix.var=TRUE,noninfo=TRUE,
#                                  seed.num=369,nsim=10000,
#                                  stop.criterion=10^-3,
#                                  logic.go='and',logic.nogo='or')
#
