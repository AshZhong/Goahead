library(shiny)

findSAScutoff<-function(n=50,CT=0.2,rate=0.8,direction='Greater'){
  lambdat=-log(CT)
  var=(lambdat)^2*CT^2/n
  sd=sqrt(var)
  if(direction=='Greater'){
    p=CT+qnorm(rate)*sd
  }
  if(direction=='Less'){
    p=CT-qnorm(rate)*sd
  }
  return(p)
}


SAS_Survival_Cutoff<-function(npatients=126,
                              a=1,b=1,
                              CT1.go=0.25,
                              false.go.CT1=TRUE,FGR.CT1=0.25,
                              CT1.nogo=0.25,
                              false.nogo.CT1=TRUE,FNGR.CT1=0.25,
                              CT2.go=0.3,
                              false.go.CT2=TRUE, FGR.CT2=0.5,
                              CT2.nogo=0.3,
                              false.nogo.CT2=TRUE, FNGR.CT2=0.5,
                              method='Frequentist',direction='Greater',
                              seed.num=369,para.exp=TRUE,
                              logic.go='and',logic.nogo='or'){

  # n=ifelse(eventinput,nevents,npatients*maturity)
  n = npatients

  est1.go=NA
  est2.go=NA
  est1.nogo=NA
  est2.nogo=NA
  flag=rep(0,4)
  overlap.flag=0
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
  if(para.exp==TRUE&method=='Frequentist'){
    if(direction=='Greater'){
      if(false.go.CT1==TRUE){
        est1.go<-findSAScutoff(n=n,CT=CT1.go,rate=1-FGR.CT1,direction='Greater')
      }
      if(false.go.CT2==TRUE){
        est2.go<-findSAScutoff(n=n,CT=CT2.go,rate=1-FGR.CT2,direction='Greater')
      }
      if(false.nogo.CT1==TRUE){
        est1.nogo<-findSAScutoff(n=n,CT=CT1.nogo,rate=1-FNGR.CT1,direction='Less')
      }
      if(false.nogo.CT2==TRUE){
        est2.nogo<-findSAScutoff(n=n,CT=CT2.nogo,rate=1-FNGR.CT2,direction='Less')
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
      if(false.go.CT1==TRUE){
        est1.go<-findSAScutoff(n=n,CT=CT1.go,rate=1-FGR.CT1,direction='Less')
      }
      if(false.go.CT2==TRUE){
        est2.go<-findSAScutoff(n=n,CT=CT2.go,rate=1-FGR.CT2,direction='Less')
      }
      if(false.nogo.CT1==TRUE){
        est1.nogo<-findSAScutoff(n=n,CT=CT1.nogo,rate=1-FNGR.CT1,direction='Greater')
      }
      if(false.nogo.CT2==TRUE){
        est2.nogo<-findSAScutoff(n=n,CT=CT2.nogo,rate=1-FNGR.CT2,direction='Greater')
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
  if(para.exp==FALSE){
    # The same as SAB
      temp<-SAB_Bin_Cutoff(n=n,a=a,b=b,
                           CT1.go=CT1.go,
                           false.go.CT1=false.go.CT1,FGR.CT1=FGR.CT1,
                           CT1.nogo=CT1.nogo,
                           false.nogo.CT1=false.nogo.CT1,FNGR.CT1=FNGR.CT1,
                           CT2.go=CT2.go,
                           false.go.CT2=false.go.CT2, FGR.CT2=FGR.CT2,
                           CT2.nogo=CT2.nogo,
                           false.nogo.CT2=false.nogo.CT2, FNGR.CT2=FNGR.CT2,
                           method=method,direction=direction,
                           logic.go=logic.go,logic.nogo=logic.nogo)
      return(list(cutoff=temp$cutoff/n,flag=temp$flag,overlap=temp$overlap))
      
    # if(method=='Bayesian'){
    #   temp<-SAB_Bin_Cutoff(n=n,a=a,b=b,
    #                        CT1.go=CT1.go,
    #                        false.go.CT1=false.go.CT1,FGR.CT1=FGR.CT1,
    #                        CT1.nogo=CT1.nogo,
    #                        false.nogo.CT1=false.nogo.CT1,FNGR.CT1=FNGR.CT1,
    #                        CT2.go=CT2.go,
    #                        false.go.CT2=false.go.CT2, FGR.CT2=FGR.CT2,
    #                        CT2.nogo=CT2.nogo,
    #                        false.nogo.CT2=false.nogo.CT2, FNGR.CT2=FNGR.CT2,
    #                        method=method,direction=direction,
    #                        logic.go=logic.go,logic.nogo=logic.nogo)
    #   return(list(cutoff=temp$cutoff/n,flag=temp$flag,overlap=temp$overlap))
    # }
    # if(method=='Frequentist'){
    #   if(direction=='Greater'){
    # 
    #     if(false.go.CT1==TRUE){
    #       sd.n1=sqrt(CT1.go*(1-CT1.go)/n)
    #       est1.go<-qnorm(1-FGR.CT1,mean=CT1.go,sd=sd.n1)
    #     }
    #     if(false.go.CT2==TRUE){
    #       sd.n2=sqrt(CT2.go*(1-CT2.go)/n)
    #       est2.go<-qnorm(1-FGR.CT2,mean=CT2.go,sd=sd.n2)
    #     }
    #     if(false.nogo.CT1==TRUE){
    #       sd.n1=sqrt(CT1.nogo*(1-CT1.nogo)/n)
    #       est1.nogo<-qnorm(FNGR.CT1,mean=CT1.nogo,sd=sd.n1)
    #     }
    #     if(false.nogo.CT2==TRUE){
    #       sd.n2=sqrt(CT2.nogo*(1-CT2.nogo)/n)
    #       est2.nogo<-qnorm(FNGR.CT2,mean=CT2.nogo,sd=sd.n2)
    #     }
    # 
    # 
    #     if(any(is.na(c(est1.go,est2.go)))){logic.go='and'}
    #     if(any(is.na(c(est1.nogo,est2.nogo)))){logic.nogo='and'}
    #     if(logic.go=='and'){
    #       go_cutoff=max(est1.go,est2.go,na.rm=TRUE)
    #     }
    #     if(logic.go=='or'){
    #       go_cutoff=min(est1.go,est2.go,na.rm=TRUE)
    #     }
    #     if(logic.nogo=='and')
    #     {
    #       nogo_cutoff=min(est1.nogo,est2.nogo,na.rm=TRUE)
    #     }
    #     if(logic.nogo=='or')
    #     {
    #       nogo_cutoff=max(est1.nogo,est2.nogo,na.rm=TRUE)
    #     }
    # 
    # 
    #     if(go_cutoff>=nogo_cutoff){return(list(cutoff=c(go_cutoff,nogo_cutoff),flag=flag,overlap=overlap.flag))}
    #     else{
    #       overlap.flag=1
    #       return(list(cutoff=c(go_cutoff,nogo_cutoff),flag=flag,overlap=overlap.flag))}
    #     if(all(flag==0)==FALSE){return(list(cutoff=c(NA,NA),flag=flag,overlap=overlap.flag))}
    # 
    #   }
    #   if(direction=='Less'){
    #     if(false.go.CT1==TRUE){
    #       sd.n1=sqrt(CT1.go*(1-CT1.go)/n)
    #       est1.go<-qnorm(FGR.CT1,mean=CT1.go,sd=sd.n1)
    #     }
    #     if(false.go.CT2==TRUE){
    #       sd.n2=sqrt(CT2.go*(1-CT2.go)/n)
    #       est2.go<-qnorm(FGR.CT2,mean=CT2.go,sd=sd.n2)
    #     }
    #     if(false.nogo.CT1==TRUE){
    #       sd.n1=sqrt(CT1.nogo*(1-CT1.nogo)/n)
    #       est1.nogo<-qnorm(1-FNGR.CT1,mean=CT1.nogo,sd=sd.n1)
    #     }
    #     if(false.nogo.CT2==TRUE){
    #       sd.n2=sqrt(CT2.nogo*(1-CT2.nogo)/n)
    #       est2.nogo<-qnorm(1-FNGR.CT2,mean=CT2.nogo,sd=sd.n2)
    #     }
    # 
    #     if(any(is.na(c(est1.go,est2.go)))){logic.go='and'}
    #     if(any(is.na(c(est1.nogo,est2.nogo)))){logic.nogo='and'}
    #     if(logic.go=='and'){
    #       go_cutoff=min(est1.go,est2.go,na.rm=TRUE)
    #     }
    #     if(logic.go=='or'){
    #       go_cutoff=max(est1.go,est2.go,na.rm=TRUE)
    #     }
    #     if(logic.nogo=='and')
    #     {
    #       nogo_cutoff=max(est1.nogo,est2.nogo,na.rm=TRUE)
    #     }
    #     if(logic.nogo=='or')
    #     {
    #       nogo_cutoff=min(est1.nogo,est2.nogo,na.rm=TRUE)
    #     }
    # 
    #     if(go_cutoff<=nogo_cutoff){return(list(cutoff=c(go_cutoff,nogo_cutoff),flag=flag,overlap=overlap.flag))}
    #     else{
    #       overlap.flag=1
    #       return(list(cutoff=c(go_cutoff,nogo_cutoff),flag=flag,overlap=overlap.flag))}
    #     if(all(flag==0)==FALSE){return(list(cutoff=c(NA,NA),flag=flag,overlap=overlap.flag))}
    #   }
    # 
    # }

  }

}
####################################################################

Fix_SS_SAS_Survival_Prob<-function(npatients=126,
                                   a=1,b=1,
                                   mean=c(0.15,0.9),
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
                                   seed.num=369,para.exp=TRUE,
                                   logic.go='and',logic.nogo='or'){
  # n=ifelse(eventinput,nevents,npatients*maturity)
  n=npatients
  
  minmean=min(mean)
  maxmean=max(mean)
  meanseq=seq(minmean,maxmean,0.05)
  go_prob<-rep(NA,length(mean))
  nogo_prob<-rep(NA,length(mean))
  inconclusive_prob<-rep(NA,length(mean))
  go_prob_plot<-rep(NA,length(meanseq))
  nogo_prob_plot<-rep(NA,length(meanseq))
  inconclusive_prob_plot<-rep(NA,length(meanseq))
  index=1
  unsatisfied.flag=0
  overlap.flag=0
  if(para.exp==FALSE){
    sd.n=sqrt(mean*(1-mean)/n)
    sd.nseq=sqrt(meanseq*(1-meanseq)/n)
  }
  if(para.exp==TRUE){
    lambdat=-log(mean)
    lambdatseq=-log(meanseq)
    var=(lambdat)^2*mean^2/n
    varseq=(lambdatseq)^2*meanseq^2/n
    sd.n=sqrt(var)
    sd.nseq=sqrt(varseq)
  }
  temp=SAS_Survival_Cutoff(npatients=n,
                           a=a,b=b,
                           CT1.go=CT1.go,
                           false.go.CT1 = false.go.CT1, FGR.CT1=FGR.CT1,
                           CT1.nogo = CT1.nogo,
                           false.nogo.CT1 = false.nogo.CT1, FNGR.CT1 = FNGR.CT1,
                           CT2.go=CT2.go,
                           false.go.CT2 = false.go.CT2, FGR.CT2=FGR.CT2,
                           CT2.nogo=CT2.nogo,
                           false.nogo.CT2 = false.nogo.CT2, FNGR.CT2 = FNGR.CT2,
                           method=method,direction=direction,
                           seed.num=seed.num,para.exp = para.exp,
                           logic.go = logic.go,logic.nogo = logic.nogo)

  ###SAS####
  true_go_cutoff<-ifelse(temp$overlap==0,temp$cutoff[1],(overlap.option=='GO')*temp$cutoff[1]+(overlap.option=='NOGO')*(temp$cutoff[2]))
  true_nogo_cutoff<-ifelse(temp$overlap==0,temp$cutoff[2],(overlap.option=='GO')*temp$cutoff[1]+(overlap.option=='NOGO')*(temp$cutoff[2]))
  
  if (para.exp==FALSE){
    if(direction=='Greater'){
      if(all(temp$flag==0)&temp$overlap==0){
        go_prob=1-pbinom(temp$cutoff[1]*n-1,n,mean)
        nogo_prob=pbinom(temp$cutoff[2]*n-1,n,mean)
        inconclusive_prob=1-go_prob-nogo_prob
        go_prob_plot=1-pbinom(temp$cutoff[1]*n-1,n,meanseq)
        nogo_prob_plot=pbinom(temp$cutoff[2]*n-1,n,meanseq)
        inconclusive_prob_plot=1-go_prob_plot-nogo_prob_plot
      }
      if(all(temp$flag==0)&temp$overlap==1&overlap.option=="GO"){
        go_prob=1-pbinom(temp$cutoff[1]*n-1,n,mean)
        nogo_prob=1-go_prob
        inconclusive_prob=1-go_prob-nogo_prob

        go_prob_plot=1-pbinom(temp$cutoff[1]*n-1,n,meanseq)
        nogo_prob_plot=1-go_prob_plot
        inconclusive_prob_plot=1-go_prob_plot-nogo_prob_plot

        overlap.flag=1
      }
      if(all(temp$flag==0)&temp$overlap==1&overlap.option=="NOGO"){
        go_prob=1-pbinom(temp$cutoff[2]*n-1,n,mean)
        nogo_prob=1-go_prob
        inconclusive_prob=1-go_prob-nogo_prob

        nogo_prob_plot=pbinom(temp$cutoff[2]*n-1,n,meanseq)
        go_prob_plot=1-nogo_prob_plot
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
        go_prob=pbinom(temp$cutoff[1]*n,n,mean)
        nogo_prob=1-pbinom(temp$cutoff[2]*n,n,mean)
        inconclusive_prob=1-go_prob-nogo_prob

        go_prob_plot=pbinom(temp$cutoff[1]*n,n,meanseq)
        nogo_prob_plot=1-pbinom(temp$cutoff[2]*n,n,meanseq)
        inconclusive_prob_plot=1-go_prob_plot-nogo_prob_plot
      }
      if(all(temp$flag==0)&temp$overlap==1&overlap.option=="GO"){
        go_prob=pbinom(temp$cutoff[1]*n,n,mean)
        nogo_prob=1-go_prob
        inconclusive_prob=1-go_prob-nogo_prob

        go_prob_plot=pbinom(temp$cutoff[1]*n,n,meanseq)
        nogo_prob_plot=1-go_prob_plot
        inconclusive_prob_plot=1-go_prob_plot-nogo_prob_plot

        overlap.flag=1
      }
      if(all(temp$flag==0)&temp$overlap==1&overlap.option=="NOGO"){
        go_prob=pbinom(temp$cutoff[2]*n,n,mean)
        nogo_prob=1-go_prob
        inconclusive_prob=1-go_prob-nogo_prob

        go_prob_plot=pbinom(temp$cutoff[2]*n,n,meanseq)
        nogo_prob_plot=1-go_prob_plot
        inconclusive_prob_plot=1-go_prob-nogo_prob_plot

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

  if (para.exp==TRUE){
    if(direction=='Greater'){
      if(all(temp$flag==0)&temp$overlap==0){
        go_prob=1-pnorm(temp$cutoff[1],mean=mean,sd=sd.n)
        nogo_prob=pnorm(temp$cutoff[2],mean=mean,sd=sd.n)
        go_prob_plot=1-pnorm(temp$cutoff[1],mean=meanseq,sd=sd.nseq)
        nogo_prob_plot=pnorm(temp$cutoff[2],mean=meanseq,sd=sd.nseq)
        if(any(temp$cutoff[1]>=1)){
          go_prob[temp$cutoff[1]>=1]=0
          go_prob_plot[temp$cutoff[1]>=1]=0
        }
        if(any(temp$cutoff[1]<=0)){
          go_prob[temp$cutoff[1]<=0]=1
          go_prob_plot[temp$cutoff[1]<=0]=1
        }
        if(any(temp$cutoff[2]<=0)){
          nogo_prob[temp$cutoff[2]<=0]=0
          nogo_prob_plot[temp$cutoff[2]<=0]=0
        }
        if(any(temp$cutoff[2]>=1)){
          nogo_prob[temp$cutoff[2]>=1]=1
          nogo_prob_plot[temp$cutoff[2]>=1]=1
        }
        inconclusive_prob_plot=1-go_prob_plot-nogo_prob_plot
        inconclusive_prob=1-go_prob-nogo_prob
      }
      if(all(temp$flag==0)&temp$overlap==1&overlap.option=="GO"){
        go_prob=1-pnorm(temp$cutoff[1],mean=mean,sd=sd.n)
        nogo_prob=pnorm(temp$cutoff[1],mean=mean,sd=sd.n)
        
        go_prob_plot=1-pnorm(temp$cutoff[1],mean=meanseq,sd=sd.nseq)
        nogo_prob_plot=pnorm(temp$cutoff[1],mean=meanseq,sd=sd.nseq)
        
        if(any(temp$cutoff[1]>=1)){
          go_prob[temp$cutoff[1]>=1]=0
          go_prob_plot[temp$cutoff[1]>=1]=0
        }
        if(any(temp$cutoff[1]<=0)){
          go_prob[temp$cutoff[1]<=0]=1
          go_prob_plot[temp$cutoff[1]<=0]=1
        }
        if(any(temp$cutoff[1]<=0)){
          nogo_prob[temp$cutoff[1]<=0]=0
          nogo_prob_plot[temp$cutoff[1]<=0]=0
        }
        if(any(temp$cutoff[1]>=1)){
          nogo_prob[temp$cutoff[1]>=1]=1
          nogo_prob_plot[temp$cutoff[1]>=1]=1
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
        
        if(any(temp$cutoff[2]>=1)){
          go_prob[temp$cutoff[2]>=1]=0
          go_prob_plot[temp$cutoff[2]>=1]=0
        }
        if(any(temp$cutoff[2]<=0)){
          go_prob[temp$cutoff[2]<=0]=1
          go_prob_plot[temp$cutoff[2]<=0]=1
        }
        if(any(temp$cutoff[2]<=0)){
          nogo_prob[temp$cutoff[2]<=0]=0
          nogo_prob_plot[temp$cutoff[2]<=0]=0
        }
        if(any(temp$cutoff[2]>=1)){
          nogo_prob[temp$cutoff[2]>=1]=1
          nogo_prob_plot[temp$cutoff[2]>=1]=1
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
        
        if(any(temp$cutoff[1]>=1)){
          go_prob[temp$cutoff[1]>=1]=1
          go_prob_plot[temp$cutoff[1]>=1]=1
        }
        if(any(temp$cutoff[1]<=0)){
          go_prob[temp$cutoff[1]<=0]=0
          go_prob_plot[temp$cutoff[1]<=0]=0
        }
        if(any(temp$cutoff[2]<=0)){
          nogo_prob[temp$cutoff[2]<=0]=1
          nogo_prob_plot[temp$cutoff[2]<=0]=1
        }
        if(any(temp$cutoff[2]>=1)){
          nogo_prob[temp$cutoff[2]>=1]=0
          nogo_prob_plot[temp$cutoff[2]>=1]=0
        }
        inconclusive_prob_plot=1-go_prob_plot-nogo_prob_plot
        inconclusive_prob=1-go_prob-nogo_prob
        
      }
      if(all(temp$flag==0)&temp$overlap==1&overlap.option=="GO"){
        go_prob=pnorm(temp$cutoff[1],mean=mean,sd=sd.n)
        nogo_prob=1-pnorm(temp$cutoff[1],mean=mean,sd=sd.n)
        
        
        go_prob_plot=pnorm(temp$cutoff[1],mean=meanseq,sd=sd.nseq)
        nogo_prob_plot=1-pnorm(temp$cutoff[1],mean=meanseq,sd=sd.nseq)
        
        if(any(temp$cutoff[1]>=1)){
          go_prob[temp$cutoff[1]>=1]=1
          go_prob_plot[temp$cutoff[1]>=1]=1
        }
        if(any(temp$cutoff[1]<=0)){
          go_prob[temp$cutoff[1]<=0]=0
          go_prob_plot[temp$cutoff[1]<=0]=0
        }
        if(any(temp$cutoff[1]<=0)){
          nogo_prob[temp$cutoff[1]<=0]=1
          nogo_prob_plot[temp$cutoff[1]<=0]=1
        }
        if(any(temp$cutoff[1]>=1)){
          nogo_prob[temp$cutoff[1]>=1]=0
          nogo_prob_plot[temp$cutoff[1]>=1]=0
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
        
        if(any(temp$cutoff[2]>=1)){
          go_prob[temp$cutoff[2]>=1]=1
          go_prob_plot[temp$cutoff[2]>=1]=1
        }
        if(any(temp$cutoff[2]<=0)){
          go_prob[temp$cutoff[2]<=0]=0
          go_prob_plot[temp$cutoff[2]<=0]=0
        }
        if(any(temp$cutoff[2]<=0)){
          nogo_prob[temp$cutoff[2]<=0]=1
          nogo_prob_plot[temp$cutoff[2]<=0]=1
        }
        if(any(temp$cutoff[2]>=1)){
          nogo_prob[temp$cutoff[2]>=1]=0
          nogo_prob_plot[temp$cutoff[2]>=1]=0
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
  

########################################################################

  if(plot.figure==TRUE){

    delta=meanseq
    p_nogo=nogo_prob_plot*100
    p_grey=inconclusive_prob_plot*100
    p_go=go_prob_plot*100

    if(true_go_cutoff!=true_nogo_cutoff){

    par(mfrow=c(1,1),pty="m",bty="l",cex=1.4)


    plot(delta,p_go,xlab=expression(paste("Survival probablity",sep="")),
         ylab="Prob of GO/ NOGO/ Inconclusive (%)",xlim=range(delta),
         ylim=c(0,100),type="n",axes=F)

    axis(1, at=seq(minmean,maxmean,0.05), labels=T)
    axis(2, at=seq(0,100,10),labels=T)
    #text(rep(-0.6,6),seq(0,100,20),seq(0,100,20),adj=1,xpd=T)
    box()
    #abline(h=seq(0,100,10),lty=3,col=rgb(0.8,0.8,0.8),lwd=1)
    points(delta,p_nogo,type="b",pch=16,col=rgb(0.9,0,0),lwd=3,lty=1)
    points(delta,p_grey,type="b",pch=4,col=rgb(0.9,0.6,0),lwd=3,lty=6)
    points(delta,p_go,type="b",pch=2,col=rgb(0,0.7,0),lwd=3,lty=2)

    # if(eventinput==TRUE){
    #   text(min(mean),140,bquote(Number~of~events==~.(nevents)),xpd=T,adj=0,cex=0.8,col=rgb(0, 0,0))
    # }
    # if(eventinput==FALSE){
    #   text(min(mean),140,bquote(Number~of~patients==~.(npatients)~';'~Maturity==~.(maturity)),xpd=T,adj=0,cex=0.8,col=rgb(0, 0,0))
    # }
    text(min(mean),140,bquote(Number~of~patients==~.(npatients)),xpd=T,adj=0,cex=0.8,col=rgb(0, 0,0))
    
    if(overlap.flag==0&unsatisfied.flag==0){
      if(direction=='Greater'){
        text(minmean,130,bquote(GO~symbol("\336")~Observed~survival~probability~symbol("\263")~.(round(temp$cutoff[1],digits=3))),xpd=T,adj=0,cex=0.8,col=rgb(0, 0.7,0))
        text(minmean,120,bquote(NOGO~symbol("\336")~Observed~survival~probability~symbol("\243")~.(round(temp$cutoff[2],digits=3))),xpd=T,adj=0,cex=0.8,col=rgb(0.9, 0,0))
        text(minmean,110,bquote(Inconclusive~symbol("\336")~anything~between~GO~and~NOGO),xpd=T,adj=0,cex=0.8,col=rgb(0.9,0.6,0))
      }
      if(direction=='Less'){
        text(minmean,130,bquote(GO~symbol("\336")~Observed~survival~probability~symbol("\243")~.(round(temp$cutoff[1],digits=3))),xpd=T,adj=0,cex=0.8,col=rgb(0, 0.7,0))
        text(minmean,120,bquote(NOGO~symbol("\336")~Observed~survival~probability~symbol("\263")~.(round(temp$cutoff[2],digits=3))),xpd=T,adj=0,cex=0.8,col=rgb(0.9, 0,0))
        text(minmean,110,bquote(Inconclusive~symbol("\336")~anything~between~GO~and~NOGO),xpd=T,adj=0,cex=0.8,col=rgb(0.9,0.6,0))
      }
    }
    if(overlap.flag==1&unsatisfied.flag==0&overlap.option=='GO'){
      if(direction=="Greater"){
        text(minmean,130,bquote(GO~symbol("\336")~Observed~survival~probability~symbol("\263")~.(round(temp$cutoff[1],digits=3))),xpd=T,adj=0,cex=0.8,col=rgb(0, 0.7,0))
        text(minmean,120,bquote(NOGO~symbol("\336")~Observed~survival~probability~symbol("\243")~.(round(temp$cutoff[1],digits=3))),xpd=T,adj=0,cex=0.8,col=rgb(0.9, 0,0))
        #text(minmean,115,paste0('Warning: GO and NOGO zones are overlaped at some sample sizes'),xpd=T,adj=0,cex=0.5,col=rgb(0.9,0.6,0))
        #text(minmean,110,paste0('The zones are classfied by cutoff of ', overlap.option),xpd=T,adj=0,cex=0.5,col=rgb(0.9,0.6,0))

      }
      if(direction=="Less"){
        text(minmean,130,bquote(GO~symbol("\336")~Observed~survival~probability~symbol("\243")~.(round(temp$cutoff[1],digits=3))),xpd=T,adj=0,cex=0.8,col=rgb(0, 0.7,0))
        text(minmean,120,bquote(NOGO~symbol("\336")~Observed~survival~probability~symbol("\263")~.(round(temp$cutoff[1],digits=3))),xpd=T,adj=0,cex=0.8,col=rgb(0.9, 0,0))
        #text(minmean,115,paste0('Warning: GO and NOGO zones are overlaped at some sample sizes'),xpd=T,adj=0,cex=0.5,col=rgb(0.9,0.6,0))
        #text(minmean,110,paste0('The zones are classfied by cutoff of ', overlap.option),xpd=T,adj=0,cex=0.5,col=rgb(0.9,0.6,0))

      }
    }

    if(overlap.flag==1&unsatisfied.flag==0&overlap.option=='NOGO'){
      if(direction=="Greater"){
        text(minmean,130,bquote(GO~symbol("\336")~Observed~survival~probability~symbol("\263")~.(round(temp$cutoff[2],digits=3))),xpd=T,adj=0,cex=0.8,col=rgb(0, 0.7,0))
        text(minmean,120,bquote(NOGO~symbol("\336")~Observed~survival~probability~symbol("\243")~.(round(temp$cutoff[2],digits=3))),xpd=T,adj=0,cex=0.8,col=rgb(0.9, 0,0))
        #text(minmean,115,paste0('Warning: GO and NOGO zones are overlaped at some sample sizes'),xpd=T,adj=0,cex=0.5,col=rgb(0.9,0.6,0))
        #text(minmean,110,paste0('The zones are classfied by cutoff of ', overlap.option),xpd=T,adj=0,cex=0.5,col=rgb(0.9,0.6,0))
      }
      if(direction=="Less"){
        text(minmean,130,bquote(GO~symbol("\336")~Observed~survival~probability~symbol("\243")~.(round(temp$cutoff[2],digits=3))),xpd=T,adj=0,cex=0.8,col=rgb(0, 0.7,0))
        text(minmean,120,bquote(NOGO~symbol("\336")~Observed~survival~probability~symbol("\263")~.(round(temp$cutoff[2],digits=3))),xpd=T,adj=0,cex=0.8,col=rgb(0.9, 0,0))
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

#######################################################################
  if(plot.figure==TRUE){

    delta=meanseq
    p_nogo=nogo_prob_plot*100
    p_grey=inconclusive_prob_plot*100
    p_go=go_prob_plot*100

    if(true_go_cutoff==true_nogo_cutoff){

    par(mfrow=c(1,1),pty="m",bty="l",cex=1.4)


    plot(delta,p_go,xlab=expression(paste("Survival probablity",sep="")),
         ylab="Prob of GO/ NOGO/ Inconclusive (%)",xlim=range(delta),
         ylim=c(0,100),type="n",axes=F)

    axis(1, at=seq(minmean,maxmean,0.05), labels=T)
    axis(2, at=seq(0,100,10),labels=T)
    #text(rep(-0.6,6),seq(0,100,20),seq(0,100,20),adj=1,xpd=T)
    box()
    #abline(h=seq(0,100,10),lty=3,col=rgb(0.8,0.8,0.8),lwd=1)
    points(delta,p_nogo,type="b",pch=16,col=rgb(0.9,0,0),lwd=3,lty=1)
    #points(delta,p_grey,type="b",pch=4,col=rgb(0.9,0.6,0),lwd=3,lty=6)
    points(delta,p_go,type="b",pch=2,col=rgb(0,0.7,0),lwd=3,lty=2)

    # if(eventinput==TRUE){
    #   text(min(mean),140,bquote(Number~of~events==~.(nevents)),xpd=T,adj=0,cex=0.8,col=rgb(0, 0,0))
    # }
    # if(eventinput==FALSE){
    #   text(min(mean),140,bquote(Number~of~patients==~.(npatients)~';'~Maturity==~.(maturity)),xpd=T,adj=0,cex=0.8,col=rgb(0, 0,0))
    # }
    text(min(mean),140,bquote(Number~of~patients==~.(npatients)),xpd=T,adj=0,cex=0.8,col=rgb(0, 0,0))
    
    if(overlap.flag==0&unsatisfied.flag==0){
      if(direction=='Greater'){
        text(minmean,130,bquote(GO~symbol("\336")~Observed~survival~probability~symbol("\263")~.(round(temp$cutoff[1],digits=3))),xpd=T,adj=0,cex=0.8,col=rgb(0, 0.7,0))
        text(minmean,120,bquote(NOGO~symbol("\336")~Observed~survival~probability~symbol("\243")~.(round(temp$cutoff[2],digits=3))),xpd=T,adj=0,cex=0.8,col=rgb(0.9, 0,0))
        #text(minmean,110,bquote(Inconclusive~symbol("\336")~anything~between~GO~and~NOGO),xpd=T,adj=0,cex=0.8,col=rgb(0.9,0.6,0))
      }
      if(direction=='Less'){
        text(minmean,130,bquote(GO~symbol("\336")~Observed~survival~probability~symbol("\243")~.(round(temp$cutoff[1],digits=3))),xpd=T,adj=0,cex=0.8,col=rgb(0, 0.7,0))
        text(minmean,120,bquote(NOGO~symbol("\336")~Observed~survival~probability~symbol("\263")~.(round(temp$cutoff[2],digits=3))),xpd=T,adj=0,cex=0.8,col=rgb(0.9, 0,0))
        #text(minmean,110,bquote(Inconclusive~symbol("\336")~anything~between~GO~and~NOGO),xpd=T,adj=0,cex=0.8,col=rgb(0.9,0.6,0))
      }
    }
    if(overlap.flag==1&unsatisfied.flag==0&overlap.option=='GO'){
      if(direction=="Greater"){
        text(minmean,130,bquote(GO~symbol("\336")~Observed~survival~probability~symbol("\263")~.(round(temp$cutoff[1],digits=3))),xpd=T,adj=0,cex=0.8,col=rgb(0, 0.7,0))
        text(minmean,120,bquote(NOGO~symbol("\336")~Observed~survival~probability~symbol("\243")~.(round(temp$cutoff[1],digits=3))),xpd=T,adj=0,cex=0.8,col=rgb(0.9, 0,0))
        #text(minmean,115,paste0('Warning: GO and NOGO zones are overlaped at some sample sizes'),xpd=T,adj=0,cex=0.5,col=rgb(0.9,0.6,0))
        #text(minmean,110,paste0('The zones are classfied by cutoff of ', overlap.option),xpd=T,adj=0,cex=0.5,col=rgb(0.9,0.6,0))

      }
      if(direction=="Less"){
        text(minmean,130,bquote(GO~symbol("\336")~Observed~survival~probability~symbol("\243")~.(round(temp$cutoff[1],digits=3))),xpd=T,adj=0,cex=0.8,col=rgb(0, 0.7,0))
        text(minmean,120,bquote(NOGO~symbol("\336")~Observed~survival~probability~symbol("\263")~.(round(temp$cutoff[1],digits=3))),xpd=T,adj=0,cex=0.8,col=rgb(0.9, 0,0))
        #text(minmean,115,paste0('Warning: GO and NOGO zones are overlaped at some sample sizes'),xpd=T,adj=0,cex=0.5,col=rgb(0.9,0.6,0))
        #text(minmean,110,paste0('The zones are classfied by cutoff of ', overlap.option),xpd=T,adj=0,cex=0.5,col=rgb(0.9,0.6,0))

      }
    }

    if(overlap.flag==1&unsatisfied.flag==0&overlap.option=='NOGO'){
      if(direction=="Greater"){
        text(minmean,130,bquote(GO~symbol("\336")~Observed~survival~probability~symbol("\263")~.(round(temp$cutoff[2],digits=3))),xpd=T,adj=0,cex=0.8,col=rgb(0, 0.7,0))
        text(minmean,120,bquote(NOGO~symbol("\336")~Observed~survival~probability~symbol("\243")~.(round(temp$cutoff[2],digits=3))),xpd=T,adj=0,cex=0.8,col=rgb(0.9, 0,0))
        #text(minmean,115,paste0('Warning: GO and NOGO zones are overlaped at some sample sizes'),xpd=T,adj=0,cex=0.5,col=rgb(0.9,0.6,0))
        #text(minmean,110,paste0('The zones are classfied by cutoff of ', overlap.option),xpd=T,adj=0,cex=0.5,col=rgb(0.9,0.6,0))
      }
      if(direction=="Less"){
        text(minmean,130,bquote(GO~symbol("\336")~Observed~survival~probability~symbol("\243")~.(round(temp$cutoff[2],digits=3))),xpd=T,adj=0,cex=0.8,col=rgb(0, 0.7,0))
        text(minmean,120,bquote(NOGO~symbol("\336")~Observed~survival~probability~symbol("\263")~.(round(temp$cutoff[2],digits=3))),xpd=T,adj=0,cex=0.8,col=rgb(0.9, 0,0))
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


#######################################################################

  return(list(go_prob=go_prob,nogo_prob=nogo_prob,
              inconclusive_prob=inconclusive_prob,
              overlap.flag=overlap.flag,overlap.option=overlap.option,
              unsatisfied.flag=unsatisfied.flag,cutoff=temp$cutoff,true_cutoff=c(true_go_cutoff,true_nogo_cutoff))
  )

}

Vary_SS_SAS_Survival_Prob<-function(npatientsmin=100,npatientsmax=500,
                                    a=1/3,b=1,
                                    mean=0.3,
                                    CT1.go=0.25,
                                    false.go.CT1=TRUE,FGR.CT1=0.25,
                                    CT1.nogo=0.3,
                                    false.nogo.CT1=TRUE,FNGR.CT1=0.25,
                                    CT2.go=0.3,
                                    false.go.CT2=TRUE, FGR.CT2=0.5,
                                    CT2.nogo=0.3,
                                    false.nogo.CT2=TRUE, FNGR.CT2=0.5,
                                    overlap.option='GO',plot.cutoff=TRUE,plot.prob=TRUE,
                                    method="Bayesian",direction="Greater",
                                    seed.num=369,para.exp=TRUE,
                                    logic.go='and',logic.nogo='or'){

  # nmin=ifelse(eventinput,neventsmin,npatientsmin*maturity)
  # nmax=ifelse(eventinput,neventsmax,npatientsmax*maturity)
  nmin=npatientsmin
  nmax=npatientsmax
  
  nseq=unique(round(c(seq(nmin,(nmax+nmin)/2,length=6)[-6],seq((nmax+nmin)/2,nmax,length=6))))
  go_prob<-matrix(NA,ncol=length(nseq),nrow=length(mean))
  nogo_prob<-matrix(NA,ncol=length(nseq),nrow=length(mean))
  inconclusive_prob<-matrix(NA,ncol=length(nseq),nrow=length(mean))
  go_cutoff<-rep(NA,length(nseq))
  nogo_cutoff<-rep(NA,length(nseq))

  index=1
  n_unsatisfied=NA
  n_overlap=NA


  if (para.exp==FALSE){
    for(i in nseq){
      temp=SAS_Survival_Cutoff(npatients=i,
                               a=a,b=b,
                               CT1.go=CT1.go,
                               false.go.CT1 = false.go.CT1, FGR.CT1=FGR.CT1,
                               CT1.nogo = CT1.nogo,
                               false.nogo.CT1 = false.nogo.CT1, FNGR.CT1 = FNGR.CT1,
                               CT2.go=CT2.go,
                               false.go.CT2 = false.go.CT2, FGR.CT2=FGR.CT2,
                               CT2.nogo=CT2.nogo,
                               false.nogo.CT2 = false.nogo.CT2, FNGR.CT2 = FNGR.CT2,
                               method=method,direction=direction,
                               seed.num=seed.num,para.exp = para.exp,
                               logic.go=logic.go,logic.nogo=logic.nogo)
      overlap = temp$overlap
      if(direction=='Greater'){
        if(all(temp$flag==0)&temp$overlap==0){
          go_prob[,index]=1-pbinom(temp$cutoff[1]*i-1,i,mean)
          nogo_prob[,index]=pbinom(temp$cutoff[2]*i-1,i,mean)
          inconclusive_prob[,index]=1-go_prob[,index]-nogo_prob[,index]
          go_cutoff[index]=temp$cutoff[1]
          nogo_cutoff[index]=temp$cutoff[2]
        }
        if(all(temp$flag==0)&temp$overlap==1&overlap.option=="GO"){
          go_prob[,index]=1-pbinom(temp$cutoff[1]*i-1,i,mean)
          nogo_prob[,index]=pbinom(temp$cutoff[1]*i-1,i,mean)
          inconclusive_prob[,index]=0
          go_cutoff[index]=temp$cutoff[1]
          nogo_cutoff[index]=temp$cutoff[2]
          n_overlap=c(n_overlap,i)
        }
        if(all(temp$flag==0)&temp$overlap==1&overlap.option=="NOGO"){
          go_prob[,index]=1-pbinom(temp$cutoff[2]*i-1,i,mean)
          nogo_prob[,index]=pbinom(temp$cutoff[2]*i-1,i,mean)
          inconclusive_prob[,index]=0
          go_cutoff[index]=temp$cutoff[1]
          nogo_cutoff[index]=temp$cutoff[2]
          n_overlap=c(n_overlap,i)
        }
        if(all(temp$flag==0)==FALSE){
          n_unsatisfied=c(n_unsatisfied,i)
        }
        index=index+1
      }

      if(direction=='Less'){
        if(all(temp$flag==0)&temp$overlap==0){
          go_prob[,index]=pbinom(temp$cutoff[1]*i,i,mean)
          nogo_prob[,index]=1-pbinom(temp$cutoff[2]*i,i,mean)
          inconclusive_prob[,index]=1-go_prob[,index]-nogo_prob[,index]
          go_cutoff[index]=temp$cutoff[1]
          nogo_cutoff[index]=temp$cutoff[2]
        }
        if(all(temp$flag==0)&temp$overlap==1&overlap.option=="GO"){
          go_prob[,index]=pbinom(temp$cutoff[1]*i,i,mean)
          nogo_prob[,index]=1-pbinom(temp$cutoff[1]*i,i,mean)
          inconclusive_prob[,index]=1-go_prob[,index]-nogo_prob[,index]
          go_cutoff[index]=temp$cutoff[1]
          nogo_cutoff[index]=temp$cutoff[2]
          n_overlap=c(n_overlap,i)
        }
        if(all(temp$flag==0)&temp$overlap==1&overlap.option=="NOGO"){
          go_prob[,index]=pbinom(temp$cutoff[2]*i,i,mean)
          nogo_prob[,index]=1-pbinom(temp$cutoff[2]*i,i,mean)
          inconclusive_prob[,index]=1-go_prob[,index]-nogo_prob[,index]
          go_cutoff[index]=temp$cutoff[1]
          nogo_cutoff[index]=temp$cutoff[2]
          n_overlap=c(n_overlap,i)
        }
        if(all(temp$flag==0)==FALSE){
          n_unsatisfied=c(n_unsatisfied,i)
        }
        index=index+1
      }

    }
  }
  
  if (para.exp==TRUE){
    
    ncore <- detectCores()
    cl<-makeCluster(ncore)
    registerDoParallel(cl)
    results<-foreach(i = nseq,.export=c('SAS_Survival_Cutoff',"SAB_Bin_Cutoff"),.combine=rbind) %dopar% {
      temp=SAS_Survival_Cutoff(npatients=i,
                               a=a,b=b,
                               CT1.go=CT1.go,
                               false.go.CT1 = false.go.CT1, FGR.CT1=FGR.CT1,
                               CT1.nogo = CT1.nogo,
                               false.nogo.CT1 = false.nogo.CT1, FNGR.CT1 = FNGR.CT1,
                               CT2.go=CT2.go,
                               false.go.CT2 = false.go.CT2, FGR.CT2=FGR.CT2,
                               CT2.nogo=CT2.nogo,
                               false.nogo.CT2 = false.nogo.CT2, FNGR.CT2 = FNGR.CT2,
                               method=method,direction=direction,
                               seed.num=seed.num,para.exp = para.exp,
                               logic.go=logic.go,logic.nogo=logic.nogo)
      c(temp$overlap,temp$cutoff)
    }
    stopCluster(cl)
    
    go_cutoff<-results[,2]
    nogo_cutoff<-results[,3]
    overlap<-results[,1]
    for(mean.index in 1:length(mean)){
      mean_temp=mean[mean.index]
      sd.i=sqrt(mean_temp*(1-mean_temp)/nseq)
      if(para.exp==FALSE){
        sd.i=sqrt(mean_temp*(1-mean_temp)/nseq)
      }
      if(para.exp==TRUE){
        lambdat=-log(mean_temp)
        var=(lambdat)^2*mean_temp^2/nseq
        sd.i=sqrt(var)
        
      }
      
      
      true_go_cutoff<-go_cutoff
      true_nogo_cutoff<-nogo_cutoff
      
      
      if(direction=='Greater'){
        true_go_cutoff[overlap==1]<-(overlap.option=='GO')*go_cutoff[overlap==1]+(overlap.option=='NOGO')*nogo_cutoff[overlap==1]
        true_nogo_cutoff[overlap==1]<-(overlap.option=='GO')*go_cutoff[overlap==1]+(overlap.option=='NOGO')*nogo_cutoff[overlap==1]
        go_prob[mean.index,]=1-pnorm(true_go_cutoff,mean=mean_temp,sd=sd.i)
        nogo_prob[mean.index,]=pnorm(true_nogo_cutoff,mean=mean_temp,sd=sd.i)
        if(any(true_go_cutoff>=1)){
          go_prob[mean.index,true_go_cutoff>=1]=0
        }
        if(any(true_go_cutoff<=0)){
          go_prob[mean.index,cutoff<=0]=1
        }
        if(any(true_nogo_cutoff<=0)){
          nogo_prob[mean.index,true_nogo_cutoff<=0]=0
        }
        if(any(true_nogo_cutoff>=1)){
          nogo_prob[mean.index,true_nogo_cutoff>=1]=1
        }
        inconclusive_prob[mean.index,]=1-go_prob[mean.index,]-nogo_prob[mean.index,]
        
      }
      
      if(direction=='Less'){
        
        true_go_cutoff[overlap==1]<-(overlap.option=='GO')*go_cutoff[overlap==1]+(overlap.option=='NOGO')*nogo_cutoff[overlap==1]
        true_nogo_cutoff[overlap==1]<-(overlap.option=='GO')*go_cutoff[overlap==1]+(overlap.option=='NOGO')*nogo_cutoff[overlap==1]
        
        go_prob[mean.index,]=pnorm(true_go_cutoff,mean=mean_temp,sd=sd.i)
        nogo_prob[mean.index,]=1-pnorm(true_nogo_cutoff,mean=mean_temp,sd=sd.i)
        if(any(true_go_cutoff>=1)){
          go_prob[mean.index,true_go_cutoff>=1]=1
        }
        if(any(true_go_cutoff<=0)){
          go_prob[mean.index,cutoff<=0]=0
        }
        if(any(true_nogo_cutoff<=0)){
          nogo_prob[mean.index,true_nogo_cutoff<=0]=1
        }
        if(any(true_nogo_cutoff>=1)){
          nogo_prob[mean.index,true_nogo_cutoff>=1]=0
        }
        inconclusive_prob[mean.index,]=1-go_prob[mean.index,]-nogo_prob[mean.index,]
        
      }
    }
    
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
      axis(1, at=nseq, labels=T)
      axis(2, at=seq(0,100,10),labels=T)
      #text(rep(-0.6,6),seq(0,100,20),seq(0,100,20),adj=1,xpd=T)
      box()
      #abline(h=seq(0,100,10),lty=3,col=rgb(0.8,0.8,0.8),lwd=1)
      polygon(c(delta,rev(delta)),c(cum_p_nogo,rev(rep(0,length(delta)))),col=rgb(0.9,0,0),border=NA)
      polygon(c(delta,rev(delta)),c(cum_p_grey,rev(cum_p_nogo)),col=rgb(0.9,0.6,0),border=NA)
      polygon(c(delta,rev(delta)),c(cum_p_go,rev(cum_p_grey)),col=rgb(0,0.7,0),border=NA)
      text((nmin+nmax)/2,120,paste0("True survival probability=",round(mean[j],3)),xpd=T,adj=0.5,cex=0.8)
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
    ylim_max=max(c(go_cutoff,nogo_cutoff,round(CT1.go,3),round(CT2.go,3),round(CT1.nogo,3),round(CT2.nogo,3)),na.rm=TRUE)+0.1
    ylim_min=min(c(go_cutoff,nogo_cutoff,round(CT1.go,3),round(CT2.go,3),round(CT1.nogo,3),round(CT2.nogo,3)),na.rm=TRUE)-0.1
    plot(NA,NA,xlab='Sample size',ylab="Observed survival probability",xlim=c(max(0,range(nseq)[1]-10),range(nseq)[2]+10),ylim=c(ylim_min,ylim_max),type="n",axes=F,col=rgb(1,0,0),lty=1,lwd=2)
    axis(1, at=nseq, labels=T)
    axis(2, at=round(c(seq(ylim_min,ylim_max,round((ylim_max-ylim_min)/10,digits=2)),CT1.go,CT2.go,CT1.nogo,CT2.nogo),digits=2),labels=T)
    #text(rep(-0.6,6),seq(0,100,20),seq(0,300,20),adj=1,xpd=T)
    box()
    lines(nseq,nogo_cutoff,col=rgb(0.9,0,0),lwd=2)
    lines(nseq,go_cutoff,col=rgb(0,0.7,0),lwd=2)
    legend('bottomright',legend=c("Cut off of GO","Cut off of NOGO"),
           col=c(rgb(0,0.7,0),rgb(0.9,0,0)),
           lwd=c(2,2),
           lty=c(1,1),cex=0.5)

    if(any(overlap!=0)){
      text(nmin,ylim_max+abs(ylim_max-ylim_min)/10,paste0('Warning: GO and NOGO zones are overlaped at some sample sizes'),xpd=T,adj=0,cex=0.5,col=rgb(0.9,0.6,0))
      #text(nmin,110,paste(n_overlap),xpd=T,adj=0,cex=0.5,col=rgb(0.9,0.6,0))
      #print(n_overlap)
      text(nmin,ylim_max+abs(ylim_max-ylim_min)/20,paste0('The zones are classfied by cutoff of ', overlap.option),xpd=T,adj=0,cex=0.5,col=rgb(0.9,0.6,0))
    }



  }
  return(list(overlap))


}



Interim_SAS<-function(num_interim=3,interim_n=c(50,100,150),
                      CT1.go=c(0.25,0.25,0.25),
                      false.go.CT1=c(TRUE,TRUE,TRUE),FGR.CT1=c(0.25,0.25,0.25),
                      CT1.nogo=c(0.25,0.25,0.25),
                      false.nogo.CT1=c(TRUE,TRUE,TRUE),FNGR.CT1=c(0.25,0.25,0.25),
                      CT2.go=c(0.25,0.25,0.25),
                      false.go.CT2=c(TRUE,TRUE,TRUE),FGR.CT2=c(0.25,0.25,0.25),
                      CT2.nogo=c(0.25,0.25,0.25),
                      false.nogo.CT2=c(TRUE,TRUE,TRUE),FNGR.CT2=c(0.25,0.25,0.25),
                      overlap.option=c('GO','GO','GO'),
                      method='Bayesian',direction=c("Greater",'Greater','Greater'),nsim_IA=10000,seed.num=369,
                      task=c('Futility','Superiority','Futility and superiority'),
                      logic.go=c('and','or','and'),
                      logic.nogo=c('and','or','or'),
                      a=1/3,b=1,mean=c(0.25),para.exp=TRUE){
  # if(eventinput==FALSE){
  #   interim_n=floor(interim_n*maturity)
  # }
  interim_n=sort(interim_n)
  go_cutoff<-rep(NA,num_interim)
  nogo_cutoff<-rep(NA,num_interim)
  true_go_cutoff<-rep(NA,num_interim)
  true_nogo_cutoff<-rep(NA,num_interim)
  overlap<-rep(NA,num_interim)
  for(i in 1:num_interim){
    temp<-SAS_Survival_Cutoff(npatients=interim_n[i],a=a,b=b,
                              CT1.go=CT1.go[i],
                              false.go.CT1=false.go.CT1[i],FGR.CT1=FGR.CT1[i],
                              CT1.nogo=CT1.nogo[i],
                              false.nogo.CT1=false.nogo.CT1[i],FNGR.CT1=FNGR.CT1[i],
                              CT2.go=CT2.go[i],
                              false.go.CT2=false.go.CT2[i], FGR.CT2=FGR.CT2[i],
                              CT2.nogo=CT2.nogo[i],
                              false.nogo.CT2=false.nogo.CT2[i], FNGR.CT2=FNGR.CT2[i],
                              method=method,direction=direction[i],
                              seed.num=seed.num,para.exp=para.exp,
                              logic.go=logic.go[i],logic.nogo=logic.nogo[i])
    go_cutoff[i]<-temp$cutoff[1]
    nogo_cutoff[i]<-temp$cutoff[2]
    overlap[i]<-temp$overlap
    ###SAS####
    true_go_cutoff[i]<-ifelse(overlap[i]==0,go_cutoff[i],(overlap.option[i]=='GO')*go_cutoff[i]+(overlap.option[i]=='NOGO')*(nogo_cutoff[i]))
    true_nogo_cutoff[i]<-ifelse(overlap[i]==0,nogo_cutoff[i],(overlap.option[i]=='GO')*(go_cutoff[i])+(overlap.option[i]=='NOGO')*(nogo_cutoff[i]))

    ####
  }




  #####
  temptable=c()
  for(meanindex in 1:length(mean)){

    set.seed(seed.num)
    
    ###SAS####
    
    sim_data<-matrix(NA,nrow=nsim_IA,ncol=max(interim_n))
    
    # if(para.exp==FALSE){
    #   sd=sqrt(mean[meanindex]*(1-mean[meanindex]))
    # }
    if(para.exp==FALSE){
      for(k in 1:max(interim_n)){
        sim_data[,k]<- rbinom(nsim_IA,1,mean[meanindex])
      }
    }
    if(para.exp==TRUE){
      lambdat=-log(mean[meanindex])
      var=(lambdat)^2*mean[meanindex]^2
      sd=sqrt(var)
      for(k in 1:max(interim_n)){
        sim_data[,k]<- rnorm(nsim_IA,mean=mean[meanindex],sd=sd)
      } 
    }
 

    cum_sim_data_temp<-apply(sim_data,1,cumsum)
    cum_sim_data<-t(cum_sim_data_temp[interim_n,]/interim_n)
    #####
    go_matrix<-matrix(NA,nrow=nsim_IA,ncol=num_interim)
    nogo_matrix<-matrix(NA,nrow=nsim_IA,ncol=num_interim)
    inconclusive_matrix<-matrix(NA,nrow=nsim_IA,ncol=num_interim)
    table<-matrix(NA,ncol=num_interim+1,nrow=10)
    IA_go_matrix<-matrix(NA,nrow=nsim_IA,ncol=num_interim) ###whether continue to next stage
    for(j in 1:num_interim ){
      if(direction[j]=='Greater'){
        go_matrix[,j]<-cum_sim_data[,j]>=true_go_cutoff[j]
        nogo_matrix[,j]<-cum_sim_data[,j]<true_nogo_cutoff[j]
        inconclusive_matrix[,j]<-rep(1,nsim_IA)-go_matrix[,j]-nogo_matrix[,j]
      }
    }

    for(j in 1:num_interim ){
      if(direction[j]=='Less'){
        go_matrix[,j]<-cum_sim_data[,j]<=true_go_cutoff[j]
        nogo_matrix[,j]<-cum_sim_data[,j]>true_nogo_cutoff[j]
        inconclusive_matrix[,j]<-rep(1,nsim_IA)-go_matrix[,j]-nogo_matrix[,j]
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
      if(task[j]=='Futility'){
        table[7,j]=ifelse(direction[j]=='Greater',HTML(paste0('-/<',round(true_nogo_cutoff[j],3))),HTML(paste0('-/>',round(true_nogo_cutoff[j],3))))
      }
      if(task[j]=='Superiority'){
        table[7,j]=ifelse(direction[j]=='Greater',HTML(paste0('>=',round(true_go_cutoff[j],3),'/-')),HTML(paste0('<=',round(true_go_cutoff[j],3),'/-')))
      }
      if(task[j]=='Futility and superiority'){
        table[7,j]=ifelse(direction[j]=='Greater',HTML(paste0('>=',round(true_go_cutoff[j],3),' /','<',round(true_nogo_cutoff[j],3))),HTML(paste0('<=',round(true_go_cutoff[j],3),' / ','>',round(true_nogo_cutoff[j],3))))
      }
      table[10,j]=ifelse(direction[j]=='Greater',HTML(paste0('>=',round(go_cutoff[j],3))),HTML(paste0('<=',round(go_cutoff[j],3))))
      table[9,j]=ifelse(direction[j]=='Greater',HTML(paste0('<',round(nogo_cutoff[j],3))),HTML(paste0('>',round(nogo_cutoff[j],3))))
      table[8,j]<-ifelse(overlap[j]==1,paste0('GO/NOGO zones overlapped, classified by criterion of ',overlap.option[j]),'None')

    }
    expectss<-round(sum(as.numeric(table[1,1:num_interim])*(c(as.numeric(table[4,1:num_interim-1])+as.numeric(table[6,1:num_interim-1]),as.numeric(table[5,num_interim-1])))),3)
    table[1,num_interim+1]=HTML(paste0(expectss,' (expected)'))
    table[2,num_interim+1]=mean[meanindex]
    table[3,num_interim+1]=''
    table[4,num_interim+1]=round(sum(as.numeric(table[4,1:num_interim])),3)
    table[5,num_interim+1]=round(as.numeric(table[5,num_interim]),3)
    table[6,num_interim+1]=round(sum(as.numeric(table[6,1:num_interim])),3)
    table[7,num_interim+1]=''
    table[8,num_interim+1]=''
    table[9,num_interim+1]=''
    table[10,num_interim+1]=''
    table<-as.table(table)
    tablecolname<-c(paste0('Interim analysis ',1:(num_interim-1)),'Final analysis',"Summary")
    tablerowname<-c('Sample size','True survival probability','Task','Success','To next interim/final or inconclusive',
                    'Stop','Superority/Futility zone','Warning',
                    'Cut off for NOGO rule',
                    'Cut off for GO rule')
    table<-cbind(tablerowname,rep(meanindex,10),table)
    colnames(table)<-c(" ",'Setting',tablecolname)
    temptable=rbind(temptable,table)

  }
  return(temptable)
}

# #Interim_SAS(interim_n = c(66,131,197),CT1.go=c(0.5,0.9,0.7),FGR.CT1=c(1-0.71,0.653,0.0238),FNGR.CT1=c(0.707,1-0.653,1-0.0238),CT2.go=c(NA,NA,NA,NA),method='Frequentist',direction=c('Greater','Less','Greater'),task=c('Futility','Superiority','Superiority'),mean=c(0.9,1.8))
#
# Fix_SS_SAS_Survival_Prob(eventinput=TRUE,nevents=100,
#                                    npatients=180,maturity=0.7,
#                                    a=1,b=1,
#                                    mean=c(0.15,0.9),
#                                    CT1.go=0.38,
#                                    false.go.CT1=TRUE,FGR.CT1=0.7,
#                                    CT1.nogo=0.38,
#                                    false.nogo.CT1=TRUE,FNGR.CT1=0.2,
#                                    CT2.go=0.3,
#                                    false.go.CT2=TRUE, FGR.CT2=0.5,
#                                    CT2.nogo=0.3,
#                                    false.nogo.CT2=TRUE, FNGR.CT2=0.5,
#                                    overlap.option='GO',plot.figure=TRUE,
#                                    method='Bayesian',direction='Greater',
#                                    seed.num=369,para.exp=FALSE,
#                                    logic.go='and',logic.nogo='or')


