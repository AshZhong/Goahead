library(shiny)
####Two arms Binomial####

findbetacutoff<-function(a1=1,b1=1,a2=1,b2=1,n1=50,n2=50,p2=0.3,CT=0.2,rate=0.8,nsim=10000,stop.criterion=10^-3,direction='Greater',seed.num=369){
  S=1
  p1=p2
  set.seed(seed.num)
  upper=1
  lower=0
  while(S==1){
    x1=rbeta(nsim,a1+n1,b1)
    x0=rbeta(nsim,a1,b1+n1)
    x=rbeta(nsim,a1+n1*p1,b1+n1*(1-p1))
    y=rbeta(nsim,a2+n2*p2,b2+n2*(1-p2))
    if(direction=='Greater'){
      prob1=sum(x1>=y+CT)/nsim
      prob0=sum(x0>=y+CT)/nsim
      prob=sum(x>=y+CT)/nsim
      if(prob1<=rate){
        p1=Inf
        S=0
      }
      if(prob0>=rate){
        p1=0
        S=0
      }
      if(abs(prob-rate)<=stop.criterion){
        S=0
      }else if(prob-rate>stop.criterion){
        upper=p1
        p1=(lower+p1)/2
      }else{
        lower=p1
        p1=(upper+p1)/2
      }
    }
    if(direction=='Less'){
      prob1=sum(x1<=y+CT)/nsim
      prob0=sum(x0<=y+CT)/nsim
      prob=sum(x<=y+CT)/nsim
      if(prob0<=rate){
        p1=-Inf
        S=0
      }
      if(prob1>=rate){
        p1=1
        S=0
      }
      if(abs(prob-rate)<=stop.criterion){
        S=0
      }else if(prob-rate>stop.criterion){
        lower=p1
        p1=(upper+p1)/2
      }else{
        upper=p1
        p1=(lower+p1)/2
      }
    }
  }
  return(p1)
}
findtabcutoff<-function(n1=50,n2=50,p2=0.3,CT=0.2,rate=0.8,stop.criterion=10^-5,direction='Greater'){
  upper=1
  lower=0
  S=1
  p1=p2
  if(direction=='Greater'){
    while(S==1){
      if((1-p2)<(CT+qnorm(rate)*sqrt(p2*(1-p2)/n2))){
        p1=Inf
        S=0
        break
      }else if((0-p2)>=(CT+qnorm(rate)*sqrt(p2*(1-p2)/n2))){
        p1=0
        S=0
        break
      }
      sd=sqrt(p1*(1-p1)/n1+p2*(1-p2)/n2)
      temp=p1-p2-CT-qnorm(rate)*sd
      if(abs(temp)<=stop.criterion){
        S=0
      }else if(temp>stop.criterion){
        upper=p1
        p1=(lower+p1)/2
      }else if(temp< -stop.criterion){
        lower=p1
        p1=(upper+p1)/2
      }
    }
  }
  if(direction=='Less'){
    while(S==1){
      if((1-p2)<=(CT+qnorm(1-rate)*sqrt(p2*(1-p2)/n2))){
        p1=1
        S=0
        break
      }else if((0-p2)>(CT+qnorm(1-rate)*sqrt(p2*(1-p2)/n2))){
        p1=-Inf
        S=0
        break
      }

      sd=sqrt(p1*(1-p1)/n1+p2*(1-p2)/n2)
      temp=p1-p2-CT-qnorm(1-rate)*sd
      if(abs(temp)<=stop.criterion){
        S=0
      }else if(temp>stop.criterion){
        upper=p1
        p1=(lower+p1)/2
      }else if(temp< -stop.criterion){
        lower=p1
        p1=(upper+p1)/2
      }

    }


  }
  return(p1)
}

TAB_Bin_Cutoff<-function(n=30,alloc.ratio=1,
                         prior.a1=1/3,prior.b1=1,
                         prior.a2=1/3,prior.b2=1,
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
                         seed.num=369,nsim=10000,stop.criterion=10^-3,
                         logic.go='and',logic.nogo='and'){
  flag=rep(0,4)
  overlap.flag=0
  n1=round(n*alloc.ratio/(1+alloc.ratio))
  n2=n-n1
  var2=control.mean*(1-control.mean)
  var11.go=(control.mean+CT1.go)*(1-control.mean-CT1.go)
  var12.go=(control.mean+CT2.go)*(1-control.mean-CT2.go)
  var11.nogo=(control.mean+CT1.nogo)*(1-control.mean-CT1.nogo)
  var12.nogo=(control.mean+CT2.nogo)*(1-control.mean-CT2.nogo)
  sd.n1.go=sqrt(var11.go/n1+var2/n2)  ###for CT1.go###
  sd.n2.go=sqrt(var12.go/n1+var2/n2)  ###for CT2.go###
  sd.n1.nogo=sqrt(var11.nogo/n1+var2/n2)  ###for CT1.go###
  sd.n2.nogo=sqrt(var12.nogo/n1+var2/n2)  ###for CT2.go###
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
    if(method=='Bayesian'){
      if(false.go.CT1==TRUE){
        p1.go<-findbetacutoff(a1=prior.a1,b1=prior.b1,
                              a2=prior.a2,b2=prior.b2,
                              n1=n1,n2=n2,p2=control.mean,CT=CT1.go,
                              rate=1-FGR.CT1,nsim=10000,
                              stop.criterion=10^-3,direction='Greater',seed.num=seed.num)
        est1.go=p1.go-control.mean
      }
      if(false.go.CT2==TRUE){
        p2.go<-findbetacutoff(a1=prior.a1,b1=prior.b1,
                              a2=prior.a2,b2=prior.b2,
                              n1=n1,n2=n2,p2=control.mean,CT=CT2.go,
                              rate=1-FGR.CT2,nsim=10000,
                              stop.criterion=10^-3,direction='Greater',seed.num=seed.num)
        est2.go=p2.go-control.mean
      }

      if(false.nogo.CT1==TRUE){
        p1.nogo<-findbetacutoff(a1=prior.a1,b1=prior.b1,
                                a2=prior.a2,b2=prior.b2,
                                n1=n1,n2=n2,p2=control.mean,CT=CT1.nogo,
                                rate=1-FNGR.CT1,nsim=nsim,
                                stop.criterion=stop.criterion,direction='Less',seed.num=seed.num)
        est1.nogo=p1.nogo-control.mean
      }
      if(false.nogo.CT2==TRUE){
        p2.nogo<-findbetacutoff(a1=prior.a1,b1=prior.b1,
                                a2=prior.a2,b2=prior.b2,
                                n1=n1,n2=n2,p2=control.mean,CT=CT2.nogo,
                                rate=1-FNGR.CT2,nsim=nsim,
                                stop.criterion=stop.criterion,direction='Less',seed.num=seed.num)
        est2.nogo=p2.nogo-control.mean
      }

    }
    if(method=='Frequentist'){
      if(false.go.CT1==TRUE){
        p1.go<-findtabcutoff(n1=n1,n2=n2,p2=control.mean,CT=CT1.go,rate=1-FGR.CT1,stop.criterion=10^-5,direction='Greater')
        est1.go=p1.go-control.mean
      }
      if(false.go.CT2==TRUE){
        p2.go<-findtabcutoff(n1=n1,n2=n2,p2=control.mean,CT=CT2.go,rate=1-FGR.CT2,stop.criterion=10^-5,direction='Greater')
        est2.go=p2.go-control.mean
      }

      if(false.nogo.CT1==TRUE){
        p1.nogo=findtabcutoff(n1=n1,n2=n2,p2=control.mean,CT=CT1.nogo,rate=1-FNGR.CT1,stop.criterion=10^-5,direction='Less')
        est1.nogo=p1.nogo-control.mean
      }
      if(false.nogo.CT2==TRUE){
        p2.nogo=findtabcutoff(n1=n1,n2=n2,p2=control.mean,CT=CT2.nogo,rate=1-FNGR.CT2,stop.criterion=10^-5,direction='Less')
        est2.nogo=p2.nogo-control.mean
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
    if(method=='Bayesian'){
      if(false.go.CT1==TRUE){
        p1.go<-findbetacutoff(a1=prior.a1,b1=prior.b1,
                              a2=prior.a2,b2=prior.b2,
                              n1=n1,n2=n2,p2=control.mean,CT=CT1.go,
                              rate=1-FGR.CT1,nsim=nsim,
                              stop.criterion=stop.criterion,direction=direction,seed.num=seed.num)
        est1.go=p1.go-control.mean
      }
      if(false.go.CT2==TRUE){
        p2.go<-findbetacutoff(a1=prior.a1,b1=prior.b1,
                              a2=prior.a2,b2=prior.b2,
                              n1=n1,n2=n2,p2=control.mean,CT=CT2.go,
                              rate=1-FGR.CT2,nsim=nsim,
                              stop.criterion=stop.criterion,direction=direction,seed.num=seed.num)
        est2.go=p2.go-control.mean
      }

      if(false.nogo.CT1==TRUE){
        p1.nogo<-findbetacutoff(a1=prior.a1,b1=prior.b1,
                                a2=prior.a2,b2=prior.b2,
                                n1=n1,n2=n2,p2=control.mean,CT=CT1.nogo,
                                rate=1-FNGR.CT1,nsim=nsim,
                                stop.criterion=stop.criterion,direction='Greater',seed.num=seed.num)
        est1.nogo=p1.nogo-control.mean
      }
      if(false.nogo.CT2==TRUE){
        p2.nogo<-findbetacutoff(a1=prior.a1,b1=prior.b1,
                                a2=prior.a2,b2=prior.b2,
                                n1=n1,n2=n2,p2=control.mean,CT=CT2.nogo,
                                rate=1-FNGR.CT2,nsim=nsim,
                                stop.criterion=stop.criterion,direction='Greater',seed.num=seed.num)
        est2.nogo=p2.nogo-control.mean
      }

    }

    if(method=='Frequentist'){
      if(false.go.CT1==TRUE){
        p1.go<-findtabcutoff(n1=n1,n2=n2,p2=control.mean,CT=CT1.go,rate=1-FGR.CT1,stop.criterion=10^-5,direction='Less')
        est1.go=p1.go-control.mean
      }
      if(false.go.CT2==TRUE){
        p2.go<-findtabcutoff(n1=n1,n2=n2,p2=control.mean,CT=CT2.go,rate=1-FGR.CT2,stop.criterion=10^-5,direction='Less')
        est2.go=p2.go-control.mean
      }

      if(false.nogo.CT1==TRUE){
        p1.nogo=findtabcutoff(n1=n1,n2=n2,p2=control.mean,CT=CT1.nogo,rate=1-FNGR.CT1,stop.criterion=10^-5,direction='Greater')
        est1.nogo=p1.nogo-control.mean
      }
      if(false.nogo.CT2==TRUE){
        p2.nogo=findtabcutoff(n1=n1,n2=n2,p2=control.mean,CT=CT2.nogo,rate=1-FNGR.CT2,stop.criterion=10^-5,direction='Greater')
        est2.nogo=p2.nogo-control.mean
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

TAB_Bayesian<-function(n=c(10,100),alloc.ratio=1,
                       prior.a1=1/3,prior.b1=1,
                       prior.a2=1/3,prior.b2=1,
                       control.mean=0,
                       mean=c(-0.5,0.8),
                       CT1.go=0.25,
                       false.go.CT1=TRUE,FGR.CT1=0.25,
                       CT1.nogo=0.25,
                       false.nogo.CT1=TRUE,FNGR.CT1=0.25,
                       CT2.go=0.3,
                       false.go.CT2=TRUE, FGR.CT2=0.5,
                       CT2.nogo=0.3,
                       false.nogo.CT2=TRUE, FNGR.CT2=0.5,
                       overlap.option='GO',
                       direction='Greater',
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
      
      sum_x1=rbinom(n_repeat,i1,treatment.mean)
      sum_x2=rbinom(n_repeat,i2,control.mean)
      
      if(is.na(CT1.go)) false.go.CT1=FALSE
      if(is.na(CT1.nogo)) false.nogo.CT1=FALSE
      if(is.na(CT2.go)) false.go.CT2=FALSE
      if(is.na(CT2.nogo)) false.nogo.CT2=FALSE
      
      if(false.go.CT1==TRUE) pp_go_1 = rep(0,n_repeat)
      if(false.nogo.CT1==TRUE) pp_nogo_1 = rep(0,n_repeat)
      if(false.go.CT2==TRUE) pp_go_2 = rep(0,n_repeat)
      if(false.nogo.CT2==TRUE) pp_nogo_2 = rep(0,n_repeat)
      
      for (i in 1:n_repeat){
        temp1 = rbeta(nsim,prior.a1+sum_x1[i],prior.b1+i1-sum_x1[i])
        temp2 = rbeta(nsim,prior.a2+sum_x2[i],prior.b2+i2-sum_x2[i])
        if(direction=='Greater'){
          if(false.go.CT1==TRUE) pp_go_1[i] = sum(temp1-temp2>=CT1.go)/nsim
          if(false.nogo.CT1==TRUE) pp_nogo_1[i] = sum(temp1-temp2>=CT1.nogo)/nsim
          if(false.go.CT2==TRUE) pp_go_2[i] = sum(temp1-temp2>=CT2.go)/nsim
          if(false.nogo.CT2==TRUE) pp_nogo_2[i] = sum(temp1-temp2>=CT2.nogo)/nsim
        }
        if(direction=='Less'){
          if(false.go.CT1==TRUE) pp_go_1[i] = sum(temp1-temp2<=CT1.go)/nsim
          if(false.nogo.CT1==TRUE) pp_nogo_1[i] = sum(temp1-temp2<=CT1.nogo)/nsim
          if(false.go.CT2==TRUE) pp_go_2[i] = sum(temp1-temp2<=CT2.go)/nsim
          if(false.nogo.CT2==TRUE) pp_nogo_2[i] = sum(temp1-temp2<=CT2.nogo)/nsim
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

###############################################################
Fix_SS_TAB_Bin_Prob<-function(n=100,alloc.ratio=1,
                              prior.a1=1/3,prior.b1=1,
                              prior.a2=1/3,prior.b2=1,
                              control.mean=0,
                              mean=c(-0.5,0.8),
                              CT1.go=0.25,
                              false.go.CT1=TRUE,FGR.CT1=0.25,
                              CT1.nogo=0.25,
                              false.nogo.CT1=TRUE,FNGR.CT1=0.25,
                              CT2.go=0.3,
                              false.go.CT2=TRUE, FGR.CT2=0.5,
                              CT2.nogo=0.3,
                              false.nogo.CT2=TRUE, FNGR.CT2=0.5,
                              overlap.option='GO',plot.figure=TRUE,
                              method='Bayesian',direction='Greater',seed.num=369,
                              nsim=nsim,stop.criterion=stop.criterion,
                              logic.go='and',logic.nogo='or'){
  minmean=ifelse(min(mean)+control.mean>=0,min(mean),-control.mean)
  maxmean=ifelse(max(mean)+control.mean<=1,max(mean),1-control.mean)
  meanseq=round(seq(minmean,maxmean,0.05),3)
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
  
  var1=(mean+control.mean)*(1-mean-control.mean)
  var2=control.mean*(1-control.mean)
  var1seq=(meanseq+control.mean)*(1-meanseq-control.mean)
  sd.n=sqrt(var1/n1+var2/n2)
  sd.nseq=sqrt(var1seq/n1+var2/n2)
  
  temp=TAB_Bin_Cutoff(n=n,alloc.ratio=alloc.ratio,
                      prior.a1=prior.a1,prior.b1=prior.b1,
                      prior.a2=prior.a2,prior.b2=prior.b2,
                      control.mean=control.mean,
                      CT1.go=CT1.go,
                      false.go.CT1 = false.go.CT1, FGR.CT1=FGR.CT1,
                      CT1.nogo = CT1.nogo,
                      false.nogo.CT1 = false.nogo.CT1, FNGR.CT1 = FNGR.CT1,
                      CT2.go=CT2.go,
                      false.go.CT2 = false.go.CT2, FGR.CT2=FGR.CT2,
                      CT2.nogo=CT2.nogo,
                      false.nogo.CT2 = false.nogo.CT2, FNGR.CT2 = FNGR.CT2,
                      method=method,direction=direction,seed.num=seed.num,
                      nsim=nsim,stop.criterion=stop.criterion,
                      logic.go=logic.go,logic.nogo = logic.nogo)
  ###TAB####
  true_go_cutoff<-ifelse(temp$overlap==0,temp$cutoff[1],(overlap.option=='GO')*temp$cutoff[1]+(overlap.option=='NOGO')*(temp$cutoff[2]))
  true_nogo_cutoff<-ifelse(temp$overlap==0,temp$cutoff[2],(overlap.option=='GO')*temp$cutoff[1]+(overlap.option=='NOGO')*(temp$cutoff[2]))

  if(method=='Frequentist'){
    
    if(direction=='Greater'){
      if(all(temp$flag==0)&temp$overlap==0){
        go_prob=1-pnorm(temp$cutoff[1],mean=mean,sd=sd.n)
        nogo_prob=pnorm(temp$cutoff[2],mean=mean,sd=sd.n)
        go_prob_plot=1-pnorm(temp$cutoff[1],mean=meanseq,sd=sd.nseq)
        nogo_prob_plot=pnorm(temp$cutoff[2],mean=meanseq,sd=sd.nseq)
        if(any(control.mean+temp$cutoff[1]>=1)){
          go_prob[control.mean+temp$cutoff[1]>=1]=0
          go_prob_plot[control.mean+temp$cutoff[1]>=1]=0
        }
        if(any(control.mean+temp$cutoff[1]<=0)){
          go_prob[control.mean+temp$cutoff[1]<=0]=1
          go_prob_plot[control.mean+temp$cutoff[1]<=0]=1
        }
        if(any(control.mean+temp$cutoff[2]<=0)){
          nogo_prob[control.mean+temp$cutoff[2]<=0]=0
          nogo_prob_plot[control.mean+temp$cutoff[2]<=0]=0
        }
        if(any(control.mean+temp$cutoff[2]>=1)){
          nogo_prob[control.mean+temp$cutoff[2]>=1]=1
          nogo_prob_plot[control.mean+temp$cutoff[2]>=1]=1
        }
        inconclusive_prob_plot=1-go_prob_plot-nogo_prob_plot
        inconclusive_prob=1-go_prob-nogo_prob
      }
      if(all(temp$flag==0)&temp$overlap==1&overlap.option=="GO"){
        go_prob=1-pnorm(temp$cutoff[1],mean=mean,sd=sd.n)
        nogo_prob=pnorm(temp$cutoff[1],mean=mean,sd=sd.n)
  
        go_prob_plot=1-pnorm(temp$cutoff[1],mean=meanseq,sd=sd.nseq)
        nogo_prob_plot=pnorm(temp$cutoff[1],mean=meanseq,sd=sd.nseq)
  
        if(any(control.mean+temp$cutoff[1]>=1)){
          go_prob[control.mean+temp$cutoff[1]>=1]=0
          go_prob_plot[control.mean+temp$cutoff[1]>=1]=0
        }
        if(any(control.mean+temp$cutoff[1]<=0)){
          go_prob[control.mean+temp$cutoff[1]<=0]=1
          go_prob_plot[control.mean+temp$cutoff[1]<=0]=1
        }
        if(any(control.mean+temp$cutoff[1]<=0)){
          nogo_prob[control.mean+temp$cutoff[1]<=0]=0
          nogo_prob_plot[control.mean+temp$cutoff[1]<=0]=0
        }
        if(any(control.mean+temp$cutoff[1]>=1)){
          nogo_prob[control.mean+temp$cutoff[1]>=1]=1
          nogo_prob_plot[control.mean+temp$cutoff[1]>=1]=1
        }
        inconclusive_prob_plot=1-go_prob_plot-nogo_prob_plot
        inconclusive_prob=1-go_prob-nogo_prob
  
        overlap.flag=1
      }
      if(all(temp$flag==0)&temp$overlap==1&overlap.option=="NOGO"){
        go_prob=1-pnorm(temp$cutoff[2],mean=mean,sd=sd.n)
        nogo_prob=pnorm(temp$cutoff[2],mean=mean,sd=sd.n)
  
        go_prob_plot=1-pnorm(temp$cutoff[2],mean=meanseq,sd=sd.nseq)
        nogo_prob_plot=pnorm(temp$cutoff[2],mean=meanseq,sd=sd.nseq)
  
        if(any(control.mean+temp$cutoff[2]>=1)){
          go_prob[control.mean+temp$cutoff[2]>=1]=0
          go_prob_plot[control.mean+temp$cutoff[2]>=1]=0
        }
        if(any(control.mean+temp$cutoff[2]<=0)){
          go_prob[control.mean+temp$cutoff[2]<=0]=1
          go_prob_plot[control.mean+temp$cutoff[2]<=0]=1
        }
        if(any(control.mean+temp$cutoff[2]<=0)){
          nogo_prob[control.mean+temp$cutoff[2]<=0]=0
          nogo_prob_plot[control.mean+temp$cutoff[2]<=0]=0
        }
        if(any(control.mean+temp$cutoff[2]>=1)){
          nogo_prob[control.mean+temp$cutoff[2]>=1]=1
          nogo_prob_plot[control.mean+temp$cutoff[2]>=1]=1
        }
        inconclusive_prob_plot=1-go_prob_plot-nogo_prob_plot
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
  
    if(direction=='Less'){
      if(all(temp$flag==0)&temp$overlap==0){
        go_prob=pnorm(temp$cutoff[1],mean=mean,sd=sd.n)
        nogo_prob=1-pnorm(temp$cutoff[2],mean=mean,sd=sd.n)
  
        go_prob_plot=pnorm(temp$cutoff[1],mean=meanseq,sd=sd.nseq)
        nogo_prob_plot=1-pnorm(temp$cutoff[2],mean=meanseq,sd=sd.nseq)
  
        if(any(control.mean+temp$cutoff[1]>=1)){
          go_prob[control.mean+temp$cutoff[1]>=1]=1
          go_prob_plot[control.mean+temp$cutoff[1]>=1]=1
        }
        if(any(control.mean+temp$cutoff[1]<=0)){
          go_prob[control.mean+temp$cutoff[1]<=0]=0
          go_prob_plot[control.mean+temp$cutoff[1]<=0]=0
        }
        if(any(control.mean+temp$cutoff[2]<=0)){
          nogo_prob[control.mean+temp$cutoff[2]<=0]=1
          nogo_prob_plot[control.mean+temp$cutoff[2]<=0]=1
        }
        if(any(control.mean+temp$cutoff[2]>=1)){
          nogo_prob[control.mean+temp$cutoff[2]>=1]=0
          nogo_prob_plot[control.mean+temp$cutoff[2]>=1]=0
        }
        inconclusive_prob_plot=1-go_prob_plot-nogo_prob_plot
        inconclusive_prob=1-go_prob-nogo_prob
  
      }
      if(all(temp$flag==0)&temp$overlap==1&overlap.option=="GO"){
        go_prob=pnorm(temp$cutoff[1],mean=mean,sd=sd.n)
        nogo_prob=1-pnorm(temp$cutoff[1],mean=mean,sd=sd.n)
  
  
        go_prob_plot=pnorm(temp$cutoff[1],mean=meanseq,sd=sd.nseq)
        nogo_prob_plot=1-pnorm(temp$cutoff[1],mean=meanseq,sd=sd.nseq)
  
        if(any(control.mean+temp$cutoff[1]>=1)){
          go_prob[control.mean+temp$cutoff[1]>=1]=1
          go_prob_plot[control.mean+temp$cutoff[1]>=1]=1
        }
        if(any(control.mean+temp$cutoff[1]<=0)){
          go_prob[control.mean+temp$cutoff[1]<=0]=0
          go_prob_plot[control.mean+temp$cutoff[1]<=0]=0
        }
        if(any(control.mean+temp$cutoff[1]<=0)){
          nogo_prob[control.mean+temp$cutoff[1]<=0]=1
          nogo_prob_plot[control.mean+temp$cutoff[1]<=0]=1
        }
        if(any(control.mean+temp$cutoff[1]>=1)){
          nogo_prob[control.mean+temp$cutoff[1]>=1]=0
          nogo_prob_plot[control.mean+temp$cutoff[1]>=1]=0
        }
        inconclusive_prob_plot=1-go_prob_plot-nogo_prob_plot
        inconclusive_prob=1-go_prob-nogo_prob
  
        overlap.flag=1
      }
      if(all(temp$flag==0)&temp$overlap==1&overlap.option=="NOGO"){
        go_prob=pnorm(temp$cutoff[2],mean=mean,sd=sd.n)
        nogo_prob=1-pnorm(temp$cutoff[2],mean=mean,sd=sd.n)
  
        go_prob_plot=pnorm(temp$cutoff[2],mean=meanseq,sd=sd.nseq)
        nogo_prob_plot=1-pnorm(temp$cutoff[2],mean=meanseq,sd=sd.nseq)
  
        if(any(control.mean+temp$cutoff[2]>=1)){
          go_prob[control.mean+temp$cutoff[2]>=1]=1
          go_prob_plot[control.mean+temp$cutoff[2]>=1]=1
        }
        if(any(control.mean+temp$cutoff[2]<=0)){
          go_prob[control.mean+temp$cutoff[2]<=0]=0
          go_prob_plot[control.mean+temp$cutoff[2]<=0]=0
        }
        if(any(control.mean+temp$cutoff[2]<=0)){
          nogo_prob[control.mean+temp$cutoff[2]<=0]=1
          nogo_prob_plot[control.mean+temp$cutoff[2]<=0]=1
        }
        if(any(control.mean+temp$cutoff[2]>=1)){
          nogo_prob[control.mean+temp$cutoff[2]>=1]=0
          nogo_prob_plot[control.mean+temp$cutoff[2]>=1]=0
        }
        inconclusive_prob_plot=1-go_prob_plot-nogo_prob_plot
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
    
    temp2 = TAB_Bayesian(n=n,alloc.ratio=alloc.ratio,
                           prior.a1=prior.a1,prior.b1=prior.b1,
                           prior.a2=prior.a2,prior.b2=prior.b2,
                           control.mean=control.mean,
                           mean=mean,
                           CT1.go=CT1.go,
                           false.go.CT1=false.go.CT1,FGR.CT1=FGR.CT1,
                           CT1.nogo=CT1.nogo,
                           false.nogo.CT1=false.nogo.CT1,FNGR.CT1=FNGR.CT1,
                           CT2.go=CT2.go,
                           false.go.CT2=false.go.CT2, FGR.CT2=FGR.CT2,
                           CT2.nogo=CT2.nogo,
                           false.nogo.CT2=false.nogo.CT2, FNGR.CT2=FNGR.CT2,
                           overlap.option=overlap.option,
                           direction=direction,
                           seed.num=seed.num,nsim=nsim,n_repeat=n_repeat,
                           logic.go=logic.go,logic.nogo=logic.nogo)
    
    go_prob = as.numeric(temp2[[1]])
    nogo_prob = as.numeric(temp2[[2]])
    inconclusive_prob = as.numeric(temp2[[3]])
    
    # plot
    temp3 = TAB_Bayesian(n=n,alloc.ratio=alloc.ratio,
                         prior.a1=prior.a1,prior.b1=prior.b1,
                         prior.a2=prior.a2,prior.b2=prior.b2,
                         control.mean=control.mean,
                         mean=meanseq,
                         CT1.go=CT1.go,
                         false.go.CT1=false.go.CT1,FGR.CT1=FGR.CT1,
                         CT1.nogo=CT1.nogo,
                         false.nogo.CT1=false.nogo.CT1,FNGR.CT1=FNGR.CT1,
                         CT2.go=CT2.go,
                         false.go.CT2=false.go.CT2, FGR.CT2=FGR.CT2,
                         CT2.nogo=CT2.nogo,
                         false.nogo.CT2=false.nogo.CT2, FNGR.CT2=FNGR.CT2,
                         overlap.option=overlap.option,
                         direction=direction,
                         seed.num=seed.num,nsim=nsim,n_repeat=n_repeat,
                         logic.go=logic.go,logic.nogo=logic.nogo)
    
    go_prob_plot = as.numeric(temp3[[1]])
    nogo_prob_plot = as.numeric(temp3[[2]])
    inconclusive_prob_plot = as.numeric(temp3[[3]])
  }

#######################################################################################
  if(plot.figure==TRUE){

    delta=meanseq
    p_nogo=nogo_prob_plot*100
    p_grey=inconclusive_prob_plot*100
    p_go=go_prob_plot*100

    if(true_go_cutoff!=true_nogo_cutoff){

    par(mfrow=c(1,1),pty="m",bty="l",cex=1.4)

    plot(delta,p_go,xlab=expression(paste("True difference in Response Rates",sep="")),
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

    text(minmean,140,bquote(n==~.(n)),xpd=T,adj=0,cex=0.8,col=rgb(0, 0,0))
    if(overlap.flag==0&unsatisfied.flag==0){
      if(direction=='Greater'){
        text(minmean,130,bquote(GO~symbol("\336")~Observed~difference~'in'~RR~symbol("\263")~.(round(temp$cutoff[1],digits=3))),xpd=T,adj=0,cex=0.8,col=rgb(0, 0.7,0))
        text(minmean,120,bquote(NOGO~symbol("\336")~Observed~difference~'in'~RR~symbol("\243")~.(round(temp$cutoff[2],digits=3))),xpd=T,adj=0,cex=0.8,col=rgb(0.9, 0,0))
        text(minmean,110,bquote(Inconclusive~symbol("\336")~anything~between~GO~and~NOGO),xpd=T,adj=0,cex=0.8,col=rgb(0.9,0.6,0))
      }
      if(direction=='Less'){
        text(minmean,130,bquote(GO~symbol("\336")~Observed~difference~'in'~RR~symbol("\243")~.(round(temp$cutoff[1],digits=3))),xpd=T,adj=0,cex=0.8,col=rgb(0, 0.7,0))
        text(minmean,120,bquote(NOGO~symbol("\336")~Observed~difference~'in'~RR~symbol("\263")~.(round(temp$cutoff[2],digits=3))),xpd=T,adj=0,cex=0.8,col=rgb(0.9, 0,0))
        text(minmean,110,bquote(Inconclusive~symbol("\336")~anything~between~GO~and~NOGO),xpd=T,adj=0,cex=0.8,col=rgb(0.9,0.6,0))
      }
    }
    if(overlap.flag==1&unsatisfied.flag==0&overlap.option=='GO'){
      if(direction=="Greater"){
        text(minmean,130,bquote(GO~symbol("\336")~Observed~difference~'in'~RR~symbol("\263")~.(round(temp$cutoff[1],digits=3))),xpd=T,adj=0,cex=0.8,col=rgb(0, 0.7,0))
        text(minmean,120,bquote(NOGO~symbol("\336")~Observed~difference~'in'~RR~symbol("\243")~.(round(temp$cutoff[1],digits=3))),xpd=T,adj=0,cex=0.8,col=rgb(0.9, 0,0))
        #text(minmean,115,paste0('Warning: GO and NOGO zones are overlaped at some sample sizes'),xpd=T,adj=0,cex=0.5,col=rgb(0.9,0.6,0))
        #text(minmean,110,paste0('The zones are classfied by cutoff of ', overlap.option),xpd=T,adj=0,cex=0.5,col=rgb(0.9,0.6,0))

      }
      if(direction=="Less"){
        text(minmean,130,bquote(GO~symbol("\336")~Observed~difference~'in'~RR~symbol("\243")~.(round(temp$cutoff[1],digits=3))),xpd=T,adj=0,cex=0.8,col=rgb(0, 0.7,0))
        text(minmean,120,bquote(NOGO~symbol("\336")~Observed~difference~'in'~RR~symbol("\263")~.(round(temp$cutoff[1],digits=3))),xpd=T,adj=0,cex=0.8,col=rgb(0.9, 0,0))
        #text(minmean,115,paste0('Warning: GO and NOGO zones are overlaped at some sample sizes'),xpd=T,adj=0,cex=0.5,col=rgb(0.9,0.6,0))
        #text(minmean,110,paste0('The zones are classfied by cutoff of ', overlap.option),xpd=T,adj=0,cex=0.5,col=rgb(0.9,0.6,0))

      }
    }

    if(overlap.flag==1&unsatisfied.flag==0&overlap.option=='NOGO'){
      if(direction=="Greater"){
        text(minmean,130,bquote(GO~symbol("\336")~Observed~difference~'in'~RR~symbol("\263")~.(round(temp$cutoff[2],digits=3))),xpd=T,adj=0,cex=0.8,col=rgb(0, 0.7,0))
        text(minmean,120,bquote(NOGO~symbol("\336")~Observed~difference~'in'~RR~symbol("\243")~.(round(temp$cutoff[2],digits=3))),xpd=T,adj=0,cex=0.8,col=rgb(0.9, 0,0))
        #text(minmean,115,paste0('Warning: GO and NOGO zones are overlaped at some sample sizes'),xpd=T,adj=0,cex=0.5,col=rgb(0.9,0.6,0))
        #text(minmean,110,paste0('The zones are classfied by cutoff of ', overlap.option),xpd=T,adj=0,cex=0.5,col=rgb(0.9,0.6,0))
      }
      if(direction=="Less"){
        text(minmean,130,bquote(GO~symbol("\336")~Observed~difference~'in'~RR~symbol("\243")~.(round(temp$cutoff[2],digits=3))),xpd=T,adj=0,cex=0.8,col=rgb(0, 0.7,0))
        text(minmean,120,bquote(NOGO~symbol("\336")~Observed~difference~'in'~RR~symbol("\263")~.(round(temp$cutoff[2],digits=3))),xpd=T,adj=0,cex=0.8,col=rgb(0.9, 0,0))
        #text(minmean,115,paste0('Warning: GO and NOGO zones are overlaped at some sample sizes'),xpd=T,adj=0,cex=0.5,col=rgb(0.9,0.6,0))
        #text(minmean,110,paste0('The zones are classfied by cutoff of ', overlap.option),xpd=T,adj=0,cex=0.5,col=rgb(0.9,0.6,0))
      }
    }
    if(unsatisfied.flag==1){
      text(minmean,130,paste0('Warning: We could not find classification of zones to satisfy your both GO and NOGO criterions at some sample sizes'),xpd=T,adj=0,cex=0.5,col=rgb(0.9,0,0))
      text(minmean,120,paste0('Please check and modify your desicion rule or the range of sample size!'),xpd=T,adj=0,cex=0.5,col=rgb(0.9,0,0))

      }

    }
  }
###########################################################################

  if(plot.figure==TRUE){

    delta=meanseq
    p_nogo=nogo_prob_plot*100
    p_grey=inconclusive_prob_plot*100
    p_go=go_prob_plot*100

    if(true_go_cutoff==true_nogo_cutoff){

    par(mfrow=c(1,1),pty="m",bty="l",cex=1.4)


    plot(delta,p_go,xlab=expression(paste("True difference in Response Rates",sep="")),
         ylab="Prob of GO/ NOGO/ Inconclusive (%)",xlim=range(delta),
         ylim=c(0,100),type="n",axes=F)

    axis(1, at=meanseq, labels=meanseq)
    axis(2, at=seq(0,100,10),labels=T)
    #text(rep(-0.6,6),seq(0,100,20),seq(0,100,20),adj=1,xpd=T)
    box()
    #abline(h=seq(0,100,10),lty=3,col=rgb(0.8,0.8,0.8),lwd=1)

    points(delta,p_nogo,type="b",pch=16,col=rgb(0.9,0,0),lwd=3,lty=1)
    #points(delta,p_grey,type="b",pch=4,col=rgb(0.9,0.6,0),lwd=3,lty=6)
    points(delta,p_go,type="b",pch=2,col=rgb(0,0.7,0),lwd=3,lty=2)

    text(minmean,140,bquote(n==~.(n)),xpd=T,adj=0,cex=0.8,col=rgb(0, 0,0))
    if(overlap.flag==0&unsatisfied.flag==0){
      if(direction=='Greater'){
        text(minmean,130,bquote(GO~symbol("\336")~Observed~difference~'in'~RR~symbol("\263")~.(round(temp$cutoff[1],digits=3))),xpd=T,adj=0,cex=0.8,col=rgb(0, 0.7,0))
        text(minmean,120,bquote(NOGO~symbol("\336")~Observed~difference~'in'~RR~symbol("\243")~.(round(temp$cutoff[2],digits=3))),xpd=T,adj=0,cex=0.8,col=rgb(0.9, 0,0))
        #text(minmean,110,bquote(Inconclusive~symbol("\336")~anything~between~GO~and~NOGO),xpd=T,adj=0,cex=0.8,col=rgb(0.9,0.6,0))
      }
      if(direction=='Less'){
        text(minmean,130,bquote(GO~symbol("\336")~Observed~difference~'in'~RR~symbol("\243")~.(round(temp$cutoff[1],digits=3))),xpd=T,adj=0,cex=0.8,col=rgb(0, 0.7,0))
        text(minmean,120,bquote(NOGO~symbol("\336")~Observed~difference~'in'~RR~symbol("\263")~.(round(temp$cutoff[2],digits=3))),xpd=T,adj=0,cex=0.8,col=rgb(0.9, 0,0))
        #text(minmean,110,bquote(Inconclusive~symbol("\336")~anything~between~GO~and~NOGO),xpd=T,adj=0,cex=0.8,col=rgb(0.9,0.6,0))
      }
    }
    if(overlap.flag==1&unsatisfied.flag==0&overlap.option=='GO'){
      if(direction=="Greater"){
        text(minmean,130,bquote(GO~symbol("\336")~Observed~difference~'in'~RR~symbol("\263")~.(round(temp$cutoff[1],digits=3))),xpd=T,adj=0,cex=0.8,col=rgb(0, 0.7,0))
        text(minmean,120,bquote(NOGO~symbol("\336")~Observed~difference~'in'~RR~symbol("\243")~.(round(temp$cutoff[1],digits=3))),xpd=T,adj=0,cex=0.8,col=rgb(0.9, 0,0))
        #text(minmean,115,paste0('Warning: GO and NOGO zones are overlaped at some sample sizes'),xpd=T,adj=0,cex=0.5,col=rgb(0.9,0.6,0))
        #text(minmean,110,paste0('The zones are classfied by cutoff of ', overlap.option),xpd=T,adj=0,cex=0.5,col=rgb(0.9,0.6,0))

      }
      if(direction=="Less"){
        text(minmean,130,bquote(GO~symbol("\336")~Observed~difference~'in'~RR~symbol("\243")~.(round(temp$cutoff[1],digits=3))),xpd=T,adj=0,cex=0.8,col=rgb(0, 0.7,0))
        text(minmean,120,bquote(NOGO~symbol("\336")~Observed~difference~'in'~RR~symbol("\263")~.(round(temp$cutoff[1],digits=3))),xpd=T,adj=0,cex=0.8,col=rgb(0.9, 0,0))
        #text(minmean,115,paste0('Warning: GO and NOGO zones are overlaped at some sample sizes'),xpd=T,adj=0,cex=0.5,col=rgb(0.9,0.6,0))
        #text(minmean,110,paste0('The zones are classfied by cutoff of ', overlap.option),xpd=T,adj=0,cex=0.5,col=rgb(0.9,0.6,0))

      }
    }

    if(overlap.flag==1&unsatisfied.flag==0&overlap.option=='NOGO'){
      if(direction=="Greater"){
        text(minmean,130,bquote(GO~symbol("\336")~Observed~difference~'in'~RR~symbol("\263")~.(round(temp$cutoff[2],digits=3))),xpd=T,adj=0,cex=0.8,col=rgb(0, 0.7,0))
        text(minmean,120,bquote(NOGO~symbol("\336")~Observed~difference~'in'~RR~symbol("\243")~.(round(temp$cutoff[2],digits=3))),xpd=T,adj=0,cex=0.8,col=rgb(0.9, 0,0))
        #text(minmean,115,paste0('Warning: GO and NOGO zones are overlaped at some sample sizes'),xpd=T,adj=0,cex=0.5,col=rgb(0.9,0.6,0))
        #text(minmean,110,paste0('The zones are classfied by cutoff of ', overlap.option),xpd=T,adj=0,cex=0.5,col=rgb(0.9,0.6,0))
      }
      if(direction=="Less"){
        text(minmean,130,bquote(GO~symbol("\336")~Observed~difference~'in'~RR~symbol("\243")~.(round(temp$cutoff[2],digits=3))),xpd=T,adj=0,cex=0.8,col=rgb(0, 0.7,0))
        text(minmean,120,bquote(NOGO~symbol("\336")~Observed~difference~'in'~RR~symbol("\263")~.(round(temp$cutoff[2],digits=3))),xpd=T,adj=0,cex=0.8,col=rgb(0.9, 0,0))
        #text(minmean,115,paste0('Warning: GO and NOGO zones are overlaped at some sample sizes'),xpd=T,adj=0,cex=0.5,col=rgb(0.9,0.6,0))
        #text(minmean,110,paste0('The zones are classfied by cutoff of ', overlap.option),xpd=T,adj=0,cex=0.5,col=rgb(0.9,0.6,0))
      }
    }
    if(unsatisfied.flag==1){
      text(minmean,130,paste0('Warning: We could not find classification of zones to satisfy your both GO and NOGO criterions at some sample sizes'),xpd=T,adj=0,cex=0.5,col=rgb(0.9,0,0))
      text(minmean,120,paste0('Please check and modify your desicion rule or the range of sample size!'),xpd=T,adj=0,cex=0.5,col=rgb(0.9,0,0))

      }

    }
  }

###########################################################################


  return(list(go_prob=go_prob,nogo_prob=nogo_prob,
              inconclusive_prob=inconclusive_prob,
              overlap.flag=overlap.flag,overlap.option=overlap.option,
              unsatisfied.flag=unsatisfied.flag,cutoff=temp$cutoff,true_cutoff=c(true_go_cutoff,true_nogo_cutoff))
  )

}



Vary_SS_TAB_Bin_Prob<-function(nmin=10,nmax=50,
                               prior.a1=1/3,prior.b1=1,
                               prior.a2=1/3,prior.b2=1,
                               mean=0.3,alloc.ratio=1,
                               control.mean=0,
                               CT1.go=0.25,
                               false.go.CT1=TRUE,FGR.CT1=0.25,
                               CT1.nogo=0.25,
                               false.nogo.CT1=TRUE,FNGR.CT1=0.25,
                               CT2.go=0.3,
                               false.go.CT2=TRUE, FGR.CT2=0.5,
                               CT2.nogo=0.3,
                               false.nogo.CT2=TRUE, FNGR.CT2=0.5,
                               overlap.option='GO',plot.cutoff=TRUE,plot.prob=TRUE,
                               method="Bayesian",direction="Greater",seed.num=369,
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
  results<-foreach(i = nseq,.export=c('TAB_Bin_Cutoff','findbetacutoff',"findtabcutoff"),.combine=rbind) %dopar% {
    i1=round(i*alloc.ratio/(1+alloc.ratio))
    i2=i-i1
    var1=(mean+control.mean)*(1-mean-control.mean)
    var2=control.mean*(1-control.mean)
    sd.i=sqrt(var1/i1+var2/i2)
    temp=TAB_Bin_Cutoff(n=i,alloc.ratio=alloc.ratio,
                        prior.a1=prior.a1,prior.b1=prior.b1,
                        prior.a2=prior.a2,prior.b2=prior.b2,
                        control.mean=control.mean,
                        CT1.go=CT1.go,
                        false.go.CT1 = false.go.CT1, FGR.CT1=FGR.CT1,
                        CT1.nogo=CT1.nogo,
                        false.nogo.CT1 = false.nogo.CT1, FNGR.CT1 = FNGR.CT1,
                        CT2.go = CT2.go,
                        false.go.CT2 = false.go.CT2, FGR.CT2=FGR.CT2,
                        CT2.nogo=CT2.nogo,
                        false.nogo.CT2 = false.nogo.CT2, FNGR.CT2 = FNGR.CT2,
                        method=method,direction=direction,seed.num = seed.num,
                        nsim=nsim,stop.criterion=stop.criterion,
                        logic.go = logic.go,logic.nogo = logic.nogo)
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
      var1=(mean_temp+control.mean)*(1-mean_temp-control.mean)
      var2=control.mean*(1-control.mean)
      sd.i=sqrt(var1/nseq+var2/nseq)
      
      if(direction=='Greater'){
        go_prob[mean.index,]=1-pnorm(true_go_cutoff,mean=mean_temp,sd=sd.i)
        nogo_prob[mean.index,]=pnorm(true_nogo_cutoff,mean=mean_temp,sd=sd.i)
        if(any(control.mean+true_go_cutoff>=1)){
          go_prob[mean.index,control.mean+true_go_cutoff>=1]=0
        }
        if(any(control.mean+true_go_cutoff<=0)){
          go_prob[mean.index,control.mean+cutoff<=0]=1
        }
        if(any(control.mean+true_nogo_cutoff<=0)){
          nogo_prob[mean.index,control.mean+true_nogo_cutoff<=0]=0
        }
        if(any(control.mean+true_nogo_cutoff>=1)){
          nogo_prob[mean.index,control.mean+true_nogo_cutoff>=1]=1
        }
      }
      if(direction=='Less'){
        go_prob[mean.index,]=pnorm(true_go_cutoff,mean=mean_temp,sd=sd.i)
        nogo_prob[mean.index,]=1-pnorm(true_nogo_cutoff,mean=mean_temp,sd=sd.i)
        if(any(control.mean+true_go_cutoff>=1)){
          go_prob[mean.index,control.mean+true_go_cutoff>=1]=1
        }
        if(any(control.mean+true_go_cutoff<=0)){
          go_prob[mean.index,control.mean+cutoff<=0]=0
        }
        if(any(control.mean+true_nogo_cutoff<=0)){
          nogo_prob[mean.index,control.mean+true_nogo_cutoff<=0]=1
        }
        if(any(control.mean+true_nogo_cutoff>=1)){
          nogo_prob[mean.index,control.mean+true_nogo_cutoff>=1]=0
        }
      }
      inconclusive_prob[mean.index,]=1-go_prob[mean.index,]-nogo_prob[mean.index,]
    }
  }
  
  if(method=='Bayesian'){
    
    n_repeat=1000
    
    temp2 = TAB_Bayesian(n=nseq,alloc.ratio=alloc.ratio,
                 prior.a1=prior.a1,prior.b1=prior.b1,
                 prior.a2=prior.a2,prior.b2=prior.b2,
                 control.mean=control.mean,
                 mean=mean,
                 CT1.go=CT1.go,
                 false.go.CT1=false.go.CT1,FGR.CT1=FGR.CT1,
                 CT1.nogo=CT1.nogo,
                 false.nogo.CT1=false.nogo.CT1,FNGR.CT1=FNGR.CT1,
                 CT2.go=CT2.go,
                 false.go.CT2=false.go.CT2, FGR.CT2=FGR.CT2,
                 CT2.nogo=CT2.nogo,
                 false.nogo.CT2=false.nogo.CT2, FNGR.CT2=FNGR.CT2,
                 overlap.option=overlap.option,
                 direction=direction,
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
      axis(1, at=nseq, labels=nseq)
      axis(2, at=seq(0,100,10),labels=T)
      #text(rep(-0.6,6),seq(0,100,20),seq(0,100,20),adj=1,xpd=T)
      box()
      #abline(h=seq(0,100,10),lty=3,col=rgb(0.8,0.8,0.8),lwd=1)
      polygon(c(delta,rev(delta)),c(cum_p_nogo,rev(rep(0,length(delta)))),col=rgb(0.9,0,0),border=NA)
      polygon(c(delta,rev(delta)),c(cum_p_grey,rev(cum_p_nogo)),col=rgb(0.9,0.6,0),border=NA)
      polygon(c(delta,rev(delta)),c(cum_p_go,rev(cum_p_grey)),col=rgb(0,0.7,0),border=NA)
      text((nmin+nmax)/2,120,paste0("True difference in RR=",round(mean[j],3),';Oberserved RR in control group=',round(control.mean,3)),xpd=T,adj=0.5,cex=0.8)
      if(any(overlap!=0)){
        text(nmin,115,paste0('Warning: GO and NOGO cut-offs are overlaped at some sample sizes'),xpd=T,adj=0,cex=0.5,col=rgb(0.9,0.6,0))
        #text(nmin,110,paste(n_overlap),xpd=T,adj=0,cex=0.5,col=rgb(0.9,0.6,0))
        #print(n_overlap)
        text(nmin,110,paste0('The zones are classfied by cutoff of ', overlap.option),xpd=T,adj=0,cex=0.5,col=rgb(0.9,0.6,0))
      }

    }
  }
  if(plot.cutoff==TRUE){
    par(mfrow=c(1,1),pty="m",bty="l",cex=1.4)
    ylim_max=max(c(go_cutoff,nogo_cutoff,round(CT1.go,3),round(CT1.nogo,3),round(CT2.go,3),round(CT2.nogo,3)),na.rm=TRUE)+0.1
    ylim_min=min(c(go_cutoff,nogo_cutoff,round(CT1.go,3),round(CT1.nogo,3),round(CT2.go,3),round(CT2.nogo,3)),na.rm=TRUE)-0.1
    plot(NA,NA,xlab='Sample size',ylab="Observed difference in RR",xlim=c(max(0,range(nseq)[1]-10),range(nseq)[2]+10),ylim=c(ylim_min,ylim_max),type="n",axes=F,col=rgb(1,0,0),lty=1,lwd=2)
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
    text(nmin,ylim_max+3*abs(ylim_max-ylim_min)/20,paste0('Oberserved RR in control group=',round(control.mean,3)),xpd=T,adj=0,cex=0.5,col=rgb(0.9,0.6,0))

    if(any(overlap!=0)){
      text(nmin,ylim_max+abs(ylim_max-ylim_min)/10,paste0('Warning: GO and NOGO zones are overlaped at some sample sizes'),xpd=T,adj=0,cex=0.5,col=rgb(0.9,0.6,0))
      #text(nmin,110,paste(n_overlap),xpd=T,adj=0,cex=0.5,col=rgb(0.9,0.6,0))
      #print(n_overlap)
      text(nmin,ylim_max+abs(ylim_max-ylim_min)/20,paste0('The zones are classfied by cutoff of ', overlap.option),xpd=T,adj=0,cex=0.5,col=rgb(0.9,0.6,0))
    }



  }
  return(list(overlap))


}

Interim_TAB<-function(num_interim=3,interim_n=c(50,100,150),
                      alloc.ratio = 1,
                      prior.a1=1/3,prior.b1=1,
                      prior.a2=1/3,prior.b2=1,
                      control.mean=0,mean=c(0.25),
                      CT1.go=c(0.25,0.25,0.25),
                      false.go.CT1=c(TRUE,TRUE,TRUE),FGR.CT1=c(0.25,0.25,0.25),
                      CT1.nogo=c(0.25,0.25,0.25),
                      false.nogo.CT1=c(TRUE,TRUE,TRUE),FNGR.CT1=c(0.25,0.25,0.25),
                      CT2.go=c(0.3,0.3,0.3),
                      false.go.CT2=c(TRUE,TRUE,TRUE),FGR.CT2=c(0.25,0.25,0.25),
                      CT2.nogo=c(0.3,0.3,0.3),
                      false.nogo.CT2=c(TRUE,TRUE,TRUE),FNGR.CT2=c(0.25,0.25,0.25),
                      overlap.option=c('GO','GO','GO'),
                      method='Bayesian',direction=c("Greater",'Greater','Greater'),
                      task=c('Futility','Superiority','Futility and superiority'),
                      nsim=10000,seed.num=369,nsim_IA=10000,stop.criterion=10^-3,
                      logic.go=c('and','and','or'),logic.nogo=c('and','or','and')){
  interim_n=sort(interim_n)

  temptable=c()
  for(meanindex in 1:length(mean)){
    
    set.seed(seed.num)
    diff_interim_n<-diff(interim_n)
    generate_n<-c(interim_n[1],diff_interim_n)
    
    sum_x1=matrix(NA,nrow=nsim_IA,ncol=num_interim)
    sum_x2=sum_x1
    treatment.mean=control.mean+mean[meanindex]
    x1=c()
    x2=c()
    for(k in 1:num_interim){
      i1=round(generate_n[k]*alloc.ratio/(1+alloc.ratio))
      i2=generate_n[k]-i1
      
      if (abs(generate_n[k]*alloc.ratio/(1+alloc.ratio)-i1)>10^(-3))
        stop ("Sample size for each arm must be an integer in each interim.")
      
      x1=cbind(x1,as.matrix(rbinom(nsim_IA,i1,treatment.mean)))
      x2=cbind(x2,as.matrix(rbinom(nsim_IA,i2,control.mean)))
      
      sum_x1[,k]=apply(x1,1,sum)
      sum_x2[,k]=apply(x2,1,sum)
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
        temp<-TAB_Bin_Cutoff(n=interim_n[i],alloc.ratio = alloc.ratio,
                             prior.a1=prior.a1,prior.b1 = prior.b1,
                             prior.a2=prior.a2,prior.b2 = prior.b2,
                             control.mean=control.mean,
                             CT1.go=CT1.go[i],
                             false.go.CT1=false.go.CT1[i],FGR.CT1=FGR.CT1[i],
                             CT1.nogo=CT1.nogo[i],
                             false.nogo.CT1=false.nogo.CT1[i],FNGR.CT1=FNGR.CT1[i],
                             CT2.go=CT2.go[i],
                             false.go.CT2=false.go.CT2[i], FGR.CT2=FGR.CT2[i],
                             CT2.nogo=CT2.nogo[i],
                             false.nogo.CT2=false.nogo.CT2[i], FNGR.CT2=FNGR.CT2[i],
                             method=method,direction=direction[i],
                             seed.num=seed.num,nsim=nsim,stop.criterion=stop.criterion,
                             logic.go=logic.go[i],logic.nogo=logic.nogo[i])
        go_cutoff[i]<-temp$cutoff[1]
        nogo_cutoff[i]<-temp$cutoff[2]
        overlap[i]<-temp$overlap
        ###TAB####
        true_go_cutoff[i]<-ifelse(overlap[i]==0,go_cutoff[i],(overlap.option[i]=='GO')*go_cutoff[i]+(overlap.option[i]=='NOGO')*(nogo_cutoff[i]))
        true_nogo_cutoff[i]<-ifelse(overlap[i]==0,nogo_cutoff[i],(overlap.option[i]=='GO')*(go_cutoff[i])+(overlap.option[i]=='NOGO')*(nogo_cutoff[i]))
        ####
      }
      
      for(j in 1:num_interim ){
        i1 = round(interim_n[j]*alloc.ratio/(1+alloc.ratio))
        i2 = interim_n[j]-i1
        if(direction[j]=='Greater'){
          go_matrix[,j]<-sum_x1[,j]/i1-sum_x2[,j]/i2>=(true_go_cutoff[j])
          nogo_matrix[,j]<-sum_x1[,j]/i1-sum_x2[,j]/i2<(true_nogo_cutoff[j])
          inconclusive_matrix[,j]<-rep(1,nsim_IA)-go_matrix[,j]-nogo_matrix[,j]
        }
        if(direction[j]=='Less'){
          go_matrix[,j]<-sum_x1[,j]/i1-sum_x2[,j]/i2<=(true_go_cutoff[j])
          nogo_matrix[,j]<-sum_x1[,j]/i1-sum_x2[,j]/i2>(true_nogo_cutoff[j])
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
          temp1 = rbeta(nsim,prior.a1+sum_x1[i,j],prior.b1+i1-sum_x1[i,j])
          temp2 = rbeta(nsim,prior.a2+sum_x2[i,j],prior.b2+i2-sum_x2[i,j])
          if(direction[j]=='Greater'){
            if(false.go.CT1[j]==TRUE) pp_go_1[i] = sum(temp1-temp2>=CT1.go[j])/nsim
            if(false.nogo.CT1[j]==TRUE) pp_nogo_1[i] = sum(temp1-temp2>=CT1.nogo[j])/nsim
            if(false.go.CT2[j]==TRUE) pp_go_2[i] = sum(temp1-temp2>=CT2.go[j])/nsim
            if(false.nogo.CT2[j]==TRUE) pp_nogo_2[i] = sum(temp1-temp2>=CT2.nogo[j])/nsim
          }
          if(direction[j]=='Less'){
            if(false.go.CT1[j]==TRUE) pp_go_1[i] = sum(temp1-temp2<=CT1.go[j])/nsim
            if(false.nogo.CT1[j]==TRUE) pp_nogo_1[i] = sum(temp1-temp2<=CT1.nogo[j])/nsim
            if(false.go.CT2[j]==TRUE) pp_go_2[i] = sum(temp1-temp2<=CT2.go[j])/nsim
            if(false.nogo.CT2[j]==TRUE) pp_nogo_2[i] = sum(temp1-temp2<=CT2.nogo[j])/nsim
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
      #   table[7,j]=ifelse(direction[j]=='Greater',HTML(paste0('-/<',round(true_nogo_cutoff[1,j],3))),HTML(paste0('-/>',round(true_nogo_cutoff[1,j],3))))
      # }
      # if(task[j]=='Superiority'){
      #   table[7,j]=ifelse(direction[j]=='Greater',HTML(paste0('>=',round(true_go_cutoff[1,j],3),'/-')),HTML(paste0('<=',round(true_go_cutoff[1,j],3),'/-')))
      # }
      # if(task[j]=='Futility and superiority'){
      #   table[7,j]=ifelse(direction[j]=='Greater',HTML(paste0('>=',round(true_go_cutoff[1,j],3),' /','<',round(true_nogo_cutoff[1,j],3))),HTML(paste0('<=',round(true_go_cutoff[1,j],3),' / ','>',round(true_nogo_cutoff[1,j],3))))
      # }
      # table[10,j]=ifelse(direction[j]=='Greater',HTML(paste0('>=',round(go_cutoff[1,j],3))),HTML(paste0('<=',round(go_cutoff[1,j],3))))
      # table[9,j]=ifelse(direction[j]=='Greater',HTML(paste0('<',round(nogo_cutoff[1,j],3))),HTML(paste0('>',round(nogo_cutoff[1,j],3))))
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
    tablerowname<-c('Sample size','True difference of RR','Task','Success','To next interim/final or inconclusive',
                    'Stop')
    table<-cbind(tablerowname,rep(meanindex,6),table)
    colnames(table)<-c(" ",'Setting',tablecolname)
    temptable=rbind(temptable,table)

  }
  return(temptable)
}

#Interim_TAB(direction=c('Greater','Greater','Greater'),logic.go=c('and','or','or'),logic.nogo=c('and','or','or'))

# Fix_SS_TAB_Bin_Prob(n=45,alloc.ratio=1,
#                               prior.a1=1/3,prior.b1=1,
#                               prior.a2=1/3,prior.b2=1,
#                               control.mean=0,
#                               mean=c(-0.5,0.8),
#                               CT1.go=0.15,
#                               false.go.CT1=TRUE,FGR.CT1=0.8,
#                               CT1.nogo=0.25,
#                               false.nogo.CT1=TRUE,FNGR.CT1=0.2,
#                               CT2.go=0.3,
#                               false.go.CT2=FALSE, FGR.CT2=0.5,
#                               CT2.nogo=0.3,
#                               false.nogo.CT2=FALSE, FNGR.CT2=0.5,
#                               overlap.option='GO',plot.figure=TRUE,
#                               method='Frequentist',direction='Less',seed.num=369,
#                               nsim=2000,stop.criterion=0.00001,
#                               logic.go='and',logic.nogo='or')
#
