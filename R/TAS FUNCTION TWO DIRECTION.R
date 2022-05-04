library(shiny)
TAS_Survival_Cutoff<-function(eventinput=TRUE,nevents=126,
                              npatients=180,maturity=0.7,
                              alloc.ratio=1,
                              CT1.go=0.25,
                              false.go.CT1=TRUE,FGR.CT1=0.25,
                              CT1.nogo=0.25,
                              false.nogo.CT1=TRUE,FNGR.CT1=0.25,
                              CT2.go=0.3,
                              false.go.CT2=TRUE, FGR.CT2=0.5,
                              CT2.nogo=0.3,
                              false.nogo.CT2=TRUE, FNGR.CT2=0.5,
                              method='Bayesian',direction='Greater',
                              prior.mean=1,prior.sd=1,
                              logic.go='and',
                              logic.nogo='or'){
  if(eventinput==TRUE){n=nevents}else{n=npatients*maturity}
  n1=n*alloc.ratio/(1+alloc.ratio)
  n2=n-n1
  sd=sqrt(1/n1+1/n2)
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
  if(method=='Frequentist'){
    if(direction=='Less'){
      if(false.go.CT1==TRUE){
        est1.go=exp(qnorm(FGR.CT1,mean=log(CT1.go),sd=sd))

      }
      if(false.go.CT2==TRUE){
        est2.go=exp(qnorm(FGR.CT2,mean=log(CT2.go),sd=sd))

      }

      if(false.nogo.CT1==TRUE){
        est1.nogo=exp(qnorm(1-FNGR.CT1,mean=log(CT1.nogo),sd=sd))

      }
      if(false.nogo.CT2==TRUE){
        est2.nogo=exp(qnorm(1-FNGR.CT2,mean=log(CT2.nogo),sd=sd))
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

    if(direction=='Greater'){
      if(false.go.CT1==TRUE){
        est1.go=exp(qnorm(1-FGR.CT1,mean=log(CT1.go),sd=sd))

      }
      if(false.go.CT2==TRUE){
        est2.go=exp(qnorm(1-FGR.CT2,mean=log(CT2.go),sd=sd))

      }

      if(false.nogo.CT1==TRUE){
        est1.nogo=exp(qnorm(FNGR.CT1,mean=log(CT1.nogo),sd=sd))

      }
      if(false.nogo.CT2==TRUE){
        est2.nogo=exp(qnorm(FNGR.CT2,mean=log(CT2.nogo),sd=sd))
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
  }
  if(method=='Bayesian'){
    temp<-SAN_Normal_Cutoff(n=n,prior.mean = prior.mean,prior.sd=prior.sd,sd=sd*sqrt(n),
                            CT1.go=log(CT1.go),
                            false.go.CT1 = false.go.CT1, FGR.CT1=FGR.CT1,
                            CT1.nogo=log(CT1.nogo),
                            false.nogo.CT1 = false.nogo.CT1, FNGR.CT1 = FNGR.CT1,
                            CT2.go=log(CT2.go),
                            false.go.CT2 = false.go.CT2, FGR.CT2=FGR.CT2,
                            CT2.nogo=log(CT2.nogo),
                            false.nogo.CT2 = false.nogo.CT2, FNGR.CT2 = FNGR.CT2,
                            method='Bayesian',direction=direction,
                            fix.var=TRUE,noninfo=TRUE,
                            logic.go=logic.go,logic.nogo=logic.nogo)

    return(list(cutoff=exp(temp$cutoff),flag=temp$flag,overlap=temp$overlap))


  }
}

#####################################################
Fix_SS_TAS_Survival_Prob<-function(eventinput=TRUE,nevents=126,
                                   npatients=180,maturity=0.7,
                                   alloc.ratio=1,
                                   CT1.go=10/14.6,
                                   false.go.CT1=TRUE,FGR.CT1=0.1,
                                   CT1.nogo=10/14.6,
                                   false.nogo.CT1=TRUE,FNGR.CT1=0.1,
                                   CT2.go=10/12,
                                   false.go.CT2=TRUE, FGR.CT2=0.2,
                                   CT2.nogo=10/12,
                                   false.nogo.CT2=FALSE, FNGR.CT2=0.2,
                                   RR=c(0.933,0.483),overlap.option='GO',plot.figure=TRUE,
                                   method='Bayesian',direction='Less',prior.mean=1,prior.sd=1,
                                   logic.go='and',logic.nogo='or'){
  if(eventinput==TRUE){n=nevents}else{n=npatients*maturity}
  n1=n*alloc.ratio/(1+alloc.ratio)
  n2=n-n1
  sd=sqrt(1/n1+1/n2)
  go_prob<-rep(NA,length(RR))
  nogo_prob<-rep(NA,length(RR))
  inconclusive_prob<-rep(NA,length(RR))
  index=1
  unsatisfied.flag=0
  overlap.flag=0

  temp=TAS_Survival_Cutoff(eventinput=eventinput,nevents=nevents,
                           npatients=npatients,maturity=maturity,
                           alloc.ratio = alloc.ratio,
                           CT1.go=CT1.go,
                           false.go.CT1 = false.go.CT1, FGR.CT1=FGR.CT1,
                           CT1.nogo = CT1.nogo,
                           false.nogo.CT1 = false.nogo.CT1, FNGR.CT1 = FNGR.CT1,
                           CT2.go=CT2.go,
                           false.go.CT2 = false.go.CT2, FGR.CT2=FGR.CT2,
                           CT2.nogo=CT2.nogo,
                           false.nogo.CT2 = false.nogo.CT2, FNGR.CT2 = FNGR.CT2,
                           method=method,direction=direction,
                           prior.mean=prior.mean,prior.sd=prior.sd,
                           logic.go = logic.go,logic.nogo=logic.nogo)
  ###TAS####
  true_go_cutoff<-ifelse(temp$overlap==0,temp$cutoff[1],(overlap.option=='GO')*temp$cutoff[1]+(overlap.option=='NOGO')*(temp$cutoff[2]))
  true_nogo_cutoff<-ifelse(temp$overlap==0,temp$cutoff[2],(overlap.option=='GO')*temp$cutoff[1]+(overlap.option=='NOGO')*(temp$cutoff[2]))


  if(direction=='Less'){

    if(all(temp$flag==0)&temp$overlap==0){
      go_prob=pnorm(log(temp$cutoff[1]),mean=log(RR),sd=sd)
      nogo_prob=1-pnorm(log(temp$cutoff[2]),mean=log(RR),sd=sd)
      inconclusive_prob=1-go_prob-nogo_prob
    }
    if(all(temp$flag==0)&temp$overlap==1&overlap.option=="GO"){
      go_prob=pnorm(log(temp$cutoff[1]),mean=log(RR),sd=sd)
      nogo_prob=1-pnorm(log(temp$cutoff[1]),mean=log(RR),sd=sd)
      inconclusive_prob=1-go_prob-nogo_prob
      overlap.flag=1
    }
    if(all(temp$flag==0)&temp$overlap==1&overlap.option=="NOGO"){
      nogo_prob=1-pnorm(log(temp$cutoff[2]),mean=log(RR),sd=sd)
      go_prob=pnorm(log(temp$cutoff[2]),mean=log(RR),sd=sd)
      inconclusive_prob=1-go_prob-nogo_prob
      overlap.flag=1
    }
    if(all(temp$flag==0)==FALSE){
      unsatisfied.flag=1
      nogo_prob=rep(NA,length(RR))
      go_prob=rep(NA,length(RR))
      inconclusive_prob=rep(NA,length(RR))
    }

  }

  if(direction=='Greater'){

    if(all(temp$flag==0)&temp$overlap==0){
      go_prob=1-pnorm(log(temp$cutoff[1]),mean=log(RR),sd=sd)
      nogo_prob=pnorm(log(temp$cutoff[2]),mean=log(RR),sd=sd)
      inconclusive_prob=1-go_prob-nogo_prob
    }
    if(all(temp$flag==0)&temp$overlap==1&overlap.option=="GO"){
      go_prob=1-pnorm(log(temp$cutoff[1]),mean=log(RR),sd=sd)
      nogo_prob=pnorm(log(temp$cutoff[1]),mean=log(RR),sd=sd)
      inconclusive_prob=1-go_prob-nogo_prob
      overlap.flag=1
    }
    if(all(temp$flag==0)&temp$overlap==1&overlap.option=="NOGO"){
      go_prob=1-pnorm(log(temp$cutoff[2]),mean=log(RR),sd=sd)
      nogo_prob=pnorm(log(temp$cutoff[2]),mean=log(RR),sd=sd)
      inconclusive_prob=1-go_prob-nogo_prob
      overlap.flag=1
    }
    if(all(temp$flag==0)==FALSE){
      unsatisfied.flag=1
      nogo_prob=rep(NA,length(RR))
      go_prob=rep(NA,length(RR))
      inconclusive_prob=rep(NA,length(RR))
    }

  }

##########################################################################################

  if(plot.figure==TRUE){

    GOdirectsymbol=NA
    NOGOdirectsymbol=NA
    RRseq=seq(min(RR),max(RR),0.05)
    go_prob_plot<-rep(NA,length(RRseq))
    nogo_prob_plot<-rep(NA,length(RRseq))
    inconclusive_prob_plot<-rep(NA,length(RRseq))
    if(direction=='Greater'){
      if(all(temp$flag==0)&temp$overlap==0){
        go_prob_plot=1-pnorm(log(temp$cutoff[1]),mean=log(RRseq),sd=sd)
        nogo_prob_plot=pnorm(log(temp$cutoff[2]),mean=log(RRseq),sd=sd)
        inconclusive_prob_plot=1-go_prob_plot-nogo_prob_plot
      }
      if(all(temp$flag==0)&temp$overlap==1&overlap.option=="GO"){
        go_prob_plot=1-pnorm(log(temp$cutoff[1]),mean=log(RRseq),sd=sd)
        nogo_prob_plot=pnorm(log(temp$cutoff[1]),mean=log(RRseq),sd=sd)
        inconclusive_prob_plot=1-go_prob_plot-nogo_prob_plot
        overlap.flag=1
      }
      if(all(temp$flag==0)&temp$overlap==1&overlap.option=="NOGO"){
        go_prob_plot=1-pnorm(log(temp$cutoff[2]),mean=log(RRseq),sd=sd)
        nogo_prob_plot=pnorm(log(temp$cutoff[2]),mean=log(RRseq),sd=sd)
        inconclusive_prob_plot=1-go_prob_plot-nogo_prob_plot
        overlap.flag=1
      }
      if(all(temp$flag==0)==FALSE){
        unsatisfied.flag=1
      }

    }
    if(direction=='Less'){
      if(all(temp$flag==0)&temp$overlap==0){
        go_prob_plot=pnorm(log(temp$cutoff[1]),mean=log(RRseq),sd=sd)
        nogo_prob_plot=1-pnorm(log(temp$cutoff[2]),mean=log(RRseq),sd=sd)
        inconclusive_prob_plot=1-go_prob_plot-nogo_prob_plot
      }
      if(all(temp$flag==0)&temp$overlap==1&overlap.option=="GO"){
        go_prob_plot=pnorm(log(temp$cutoff[1]),mean=log(RRseq),sd=sd)
        nogo_prob_plot=1-pnorm(log(temp$cutoff[1]),mean=log(RRseq),sd=sd)
        inconclusive_prob_plot=1-go_prob_plot-nogo_prob_plot
        overlap.flag=1
      }
      if(all(temp$flag==0)&temp$overlap==1&overlap.option=="NOGO"){
        go_prob_plot=pnorm(log(temp$cutoff[2]),mean=log(RRseq),sd=sd)
        nogo_prob_plot=1-pnorm(log(temp$cutoff[2]),mean=log(RRseq),sd=sd)
        inconclusive_prob_plot=1-go_prob_plot-nogo_prob_plot
        overlap.flag=1
      }
      if(all(temp$flag==0)==FALSE){
        unsatisfied.flag=1
      }
    }


    delta=RRseq
    p_nogo=nogo_prob_plot*100
    p_grey=inconclusive_prob_plot*100
    p_go=go_prob_plot*100

    if(true_go_cutoff!=true_nogo_cutoff){

    par(mfrow=c(1,1),pty="m",bty="l",cex=1.4)

    plot(delta,p_go,xlab=expression(paste("True ",HR,sep="")),
         ylab="Prob of GO/ NOGO/ Inconclusive (%)",xlim=range(delta),
         ylim=c(0,100),type="n",axes=F)


    axis(1, at=seq(min(RR),max(RR),0.05), labels=T)
    axis(2, at=seq(0,100,10),labels=T)
    box()

    points(delta,p_nogo,type="b",pch=16,col=rgb(0.9,0,0),lwd=3,lty=1)
    points(delta,p_grey,type="b",pch=4,col=rgb(0.9,0.6,0),lwd=3,lty=6)
    points(delta,p_go,type="b",pch=2,col=rgb(0,0.7,0),lwd=3,lty=2)


    if(eventinput==TRUE){
      text(min(RR),140,bquote(Number~of~events==~.(nevents)),xpd=T,adj=0,cex=0.8,col=rgb(0, 0,0))
    }
    if(eventinput==FALSE){
      text(min(RR),140,bquote(Number~of~patients==~.(npatients)~';'~Maturity==~.(maturity)),xpd=T,adj=0,cex=0.8,col=rgb(0, 0,0))
    }
    if(direction=='Greater'){
      if(overlap.flag==0&unsatisfied.flag==0){
        text(min(RR),130,bquote(GO~symbol("\336")~Observed~hazard~ratio~symbol("\263")~.(round(temp$cutoff[1],digits=3))),xpd=T,adj=0,cex=0.8,col=rgb(0, 0.7,0))
        text(min(RR),120,bquote(NOGO~symbol("\336")~Observed~hazard~ratio~symbol("\243")~.(round(temp$cutoff[2],digits=3))),xpd=T,adj=0,cex=0.8,col=rgb(0.9, 0,0))
        text(min(RR),110,bquote(Inconclusive~symbol("\336")~anything~between~GO~and~NOGO),xpd=T,adj=0,cex=0.8,col=rgb(0.9,0.6,0))
      }

      if(overlap.flag==1&unsatisfied.flag==0&overlap.option=="GO"){
        text(min(RR),130,bquote(GO~symbol("\336")~Observed~hazard~ratio~symbol("\263")~.(round(temp$cutoff[1],digits=3))),xpd=T,adj=0,cex=0.8,col=rgb(0, 0.7,0))
        text(min(RR),120,bquote(NOGO~symbol("\336")~Observed~hazard~ratio~symbol("\243")~.(round(temp$cutoff[1],digits=3))),xpd=T,adj=0,cex=0.8,col=rgb(0.9, 0,0))
        #text(min(RR),115,paste0('Warning: GO and NOGO Zones are overlaped at some sample sizes'),xpd=T,adj=0,cex=0.5,col=rgb(0.9,0.6,0))
        #text(min(RR),110,paste0('The zones are classfied by cutoff of ', overlap.option),xpd=T,adj=0,cex=0.5,col=rgb(0.9,0.6,0))
      }
      if(overlap.flag==1&unsatisfied.flag==0&overlap.option=="NOGO"){
        text(min(RR),130,bquote(GO~symbol("\336")~Observed~hazard~ratio~symbol("\263")~.(round(temp$cutoff[2],digits=3))),xpd=T,adj=0,cex=0.8,col=rgb(0, 0.7,0))
        text(min(RR),120,bquote(NOGO~symbol("\336")~Observed~hazard~ratio~symbol("\243")~.(round(temp$cutoff[2],digits=3))),xpd=T,adj=0,cex=0.8,col=rgb(0.9, 0,0))
        #text(min(RR),115,paste0('Warning: GO and NOGO Zones are overlaped at some sample sizes'),xpd=T,adj=0,cex=0.5,col=rgb(0.9,0.6,0))
        #text(min(RR),110,paste0('The zones are classfied by cutoff of ', overlap.option),xpd=T,adj=0,cex=0.5,col=rgb(0.9,0.6,0))
      }
      if(unsatisfied.flag==1){
        text(min(RR),130,paste0('Warning: We could not find classification of zones to satisfy your both GO and NOGO criterions at some sample sizes'),xpd=T,adj=0,cex=0.5,col=rgb(0.9,0,0))
        text(min(RR),120,paste0('Please check and modify your desicion rule or the range of sample size!'),xpd=T,adj=0,cex=0.5,col=rgb(0.9,0,0))

      }

    }
    if(direction=='Less'){
      if(overlap.flag==0&unsatisfied.flag==0){
        text(min(RR),130,bquote(GO~symbol("\336")~Observed~hazard~ratio~symbol("\243")~.(round(temp$cutoff[1],digits=3))),xpd=T,adj=0,cex=0.8,col=rgb(0, 0.7,0))
        text(min(RR),120,bquote(NOGO~symbol("\336")~Observed~hazard~ratio~symbol("\263")~.(round(temp$cutoff[2],digits=3))),xpd=T,adj=0,cex=0.8,col=rgb(0.9, 0,0))
        text(min(RR),110,bquote(Inconclusive~symbol("\336")~anything~between~GO~and~NOGO),xpd=T,adj=0,cex=0.8,col=rgb(0.9,0.6,0))
      }

      if(overlap.flag==1&unsatisfied.flag==0&overlap.option=="GO"){
        text(min(RR),130,bquote(GO~symbol("\336")~Observed~hazard~ratio~symbol("\243")~.(round(temp$cutoff[1],digits=3))),xpd=T,adj=0,cex=0.8,col=rgb(0, 0.7,0))
        text(min(RR),120,bquote(NOGO~symbol("\336")~Observed~hazard~ratio~symbol("\263")~.(round(temp$cutoff[1],digits=3))),xpd=T,adj=0,cex=0.8,col=rgb(0.9, 0,0))
        #text(min(RR),115,paste0('Warning: GO and NOGO Zones are overlaped at some sample sizes'),xpd=T,adj=0,cex=0.5,col=rgb(0.9,0.6,0))
        #text(min(RR),110,paste0('The zones are classfied by cutoff of ', overlap.option),xpd=T,adj=0,cex=0.5,col=rgb(0.9,0.6,0))
      }
      if(overlap.flag==1&unsatisfied.flag==0&overlap.option=="NOGO"){
        text(min(RR),130,bquote(GO~symbol("\336")~Observed~hazard~ratio~symbol("\243")~.(round(temp$cutoff[2],digits=3))),xpd=T,adj=0,cex=0.8,col=rgb(0, 0.7,0))
        text(min(RR),120,bquote(NOGO~symbol("\336")~Observed~hazard~ratio~symbol("\263")~.(round(temp$cutoff[2],digits=3))),xpd=T,adj=0,cex=0.8,col=rgb(0.9, 0,0))
        #text(min(RR),115,paste0('Warning: GO and NOGO Zones are overlaped at some sample sizes'),xpd=T,adj=0,cex=0.5,col=rgb(0.9,0.6,0))
        #text(min(RR),110,paste0('The zones are classfied by cutoff of ', overlap.option),xpd=T,adj=0,cex=0.5,col=rgb(0.9,0.6,0))
      }
      if(unsatisfied.flag==1){
        text(min(RR),130,paste0('Warning: We could not find classification of zones to satisfy your both GO and NOGO criterions at some sample sizes'),xpd=T,adj=0,cex=0.5,col=rgb(0.9,0,0))
        text(min(RR),120,paste0('Please check and modify your desicion rule or the range of sample size!'),xpd=T,adj=0,cex=0.5,col=rgb(0.9,0,0))

        }

      }

     }
  }
###########################&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&

  if(plot.figure==TRUE){

    GOdirectsymbol=NA
    NOGOdirectsymbol=NA
    RRseq=seq(min(RR),max(RR),0.05)
    go_prob_plot<-rep(NA,length(RRseq))
    nogo_prob_plot<-rep(NA,length(RRseq))
    inconclusive_prob_plot<-rep(NA,length(RRseq))
    if(direction=='Greater'){
      if(all(temp$flag==0)&temp$overlap==0){
        go_prob_plot=1-pnorm(log(temp$cutoff[1]),mean=log(RRseq),sd=sd)
        nogo_prob_plot=pnorm(log(temp$cutoff[2]),mean=log(RRseq),sd=sd)
        inconclusive_prob_plot=1-go_prob_plot-nogo_prob_plot
      }
      if(all(temp$flag==0)&temp$overlap==1&overlap.option=="GO"){
        go_prob_plot=1-pnorm(log(temp$cutoff[1]),mean=log(RRseq),sd=sd)
        nogo_prob_plot=pnorm(log(temp$cutoff[1]),mean=log(RRseq),sd=sd)
        inconclusive_prob_plot=1-go_prob_plot-nogo_prob_plot
        overlap.flag=1
      }
      if(all(temp$flag==0)&temp$overlap==1&overlap.option=="NOGO"){
        go_prob_plot=1-pnorm(log(temp$cutoff[2]),mean=log(RRseq),sd=sd)
        nogo_prob_plot=pnorm(log(temp$cutoff[2]),mean=log(RRseq),sd=sd)
        inconclusive_prob_plot=1-go_prob_plot-nogo_prob_plot
        overlap.flag=1
      }
      if(all(temp$flag==0)==FALSE){
        unsatisfied.flag=1
      }

    }
    if(direction=='Less'){
      if(all(temp$flag==0)&temp$overlap==0){
        go_prob_plot=pnorm(log(temp$cutoff[1]),mean=log(RRseq),sd=sd)
        nogo_prob_plot=1-pnorm(log(temp$cutoff[2]),mean=log(RRseq),sd=sd)
        inconclusive_prob_plot=1-go_prob_plot-nogo_prob_plot
      }
      if(all(temp$flag==0)&temp$overlap==1&overlap.option=="GO"){
        go_prob_plot=pnorm(log(temp$cutoff[1]),mean=log(RRseq),sd=sd)
        nogo_prob_plot=1-pnorm(log(temp$cutoff[1]),mean=log(RRseq),sd=sd)
        inconclusive_prob_plot=1-go_prob_plot-nogo_prob_plot
        overlap.flag=1
      }
      if(all(temp$flag==0)&temp$overlap==1&overlap.option=="NOGO"){
        go_prob_plot=pnorm(log(temp$cutoff[2]),mean=log(RRseq),sd=sd)
        nogo_prob_plot=1-pnorm(log(temp$cutoff[2]),mean=log(RRseq),sd=sd)
        inconclusive_prob_plot=1-go_prob_plot-nogo_prob_plot
        overlap.flag=1
      }
      if(all(temp$flag==0)==FALSE){
        unsatisfied.flag=1
      }
    }


    delta=RRseq
    p_nogo=nogo_prob_plot*100
    p_grey=inconclusive_prob_plot*100
    p_go=go_prob_plot*100

    if(true_go_cutoff==true_nogo_cutoff){

      par(mfrow=c(1,1),pty="m",bty="l",cex=1.4)

      plot(delta,p_go,xlab=expression(paste("True ",HR,sep="")),
           ylab="Prob of GO/ NOGO/ Inconclusive (%)",xlim=range(delta),
           ylim=c(0,100),type="n",axes=F)


      axis(1, at=seq(min(RR),max(RR),0.05), labels=T)
      axis(2, at=seq(0,100,10),labels=T)
      box()

      points(delta,p_nogo,type="b",pch=16,col=rgb(0.9,0,0),lwd=3,lty=1)
      #points(delta,p_grey,type="b",pch=4,col=rgb(0.9,0.6,0),lwd=3,lty=6)
      points(delta,p_go,type="b",pch=2,col=rgb(0,0.7,0),lwd=3,lty=2)


      if(eventinput==TRUE){
        text(min(RR),140,bquote(Number~of~events==~.(nevents)),xpd=T,adj=0,cex=0.8,col=rgb(0, 0,0))
      }
      if(eventinput==FALSE){
        text(min(RR),140,bquote(Number~of~patients==~.(npatients)~';'~Maturity==~.(maturity)),xpd=T,adj=0,cex=0.8,col=rgb(0, 0,0))
      }
      if(direction=='Greater'){
        if(overlap.flag==0&unsatisfied.flag==0){
          text(min(RR),130,bquote(GO~symbol("\336")~Observed~hazard~ratio~symbol("\263")~.(round(temp$cutoff[1],digits=3))),xpd=T,adj=0,cex=0.8,col=rgb(0, 0.7,0))
          text(min(RR),120,bquote(NOGO~symbol("\336")~Observed~hazard~ratio~symbol("\243")~.(round(temp$cutoff[2],digits=3))),xpd=T,adj=0,cex=0.8,col=rgb(0.9, 0,0))
          #text(min(RR),110,bquote(Inconclusive~symbol("\336")~anything~between~GO~and~NOGO),xpd=T,adj=0,cex=0.8,col=rgb(0.9,0.6,0))
        }

        if(overlap.flag==1&unsatisfied.flag==0&overlap.option=="GO"){
          text(min(RR),130,bquote(GO~symbol("\336")~Observed~hazard~ratio~symbol("\263")~.(round(temp$cutoff[1],digits=3))),xpd=T,adj=0,cex=0.8,col=rgb(0, 0.7,0))
          text(min(RR),120,bquote(NOGO~symbol("\336")~Observed~hazard~ratio~symbol("\243")~.(round(temp$cutoff[1],digits=3))),xpd=T,adj=0,cex=0.8,col=rgb(0.9, 0,0))
          #text(min(RR),115,paste0('Warning: GO and NOGO Zones are overlaped at some sample sizes'),xpd=T,adj=0,cex=0.5,col=rgb(0.9,0.6,0))
          #text(min(RR),110,paste0('The zones are classfied by cutoff of ', overlap.option),xpd=T,adj=0,cex=0.5,col=rgb(0.9,0.6,0))
        }
        if(overlap.flag==1&unsatisfied.flag==0&overlap.option=="NOGO"){
          text(min(RR),130,bquote(GO~symbol("\336")~Observed~hazard~ratio~symbol("\263")~.(round(temp$cutoff[2],digits=3))),xpd=T,adj=0,cex=0.8,col=rgb(0, 0.7,0))
          text(min(RR),120,bquote(NOGO~symbol("\336")~Observed~hazard~ratio~symbol("\243")~.(round(temp$cutoff[2],digits=3))),xpd=T,adj=0,cex=0.8,col=rgb(0.9, 0,0))
          #text(min(RR),115,paste0('Warning: GO and NOGO Zones are overlaped at some sample sizes'),xpd=T,adj=0,cex=0.5,col=rgb(0.9,0.6,0))
          #text(min(RR),110,paste0('The zones are classfied by cutoff of ', overlap.option),xpd=T,adj=0,cex=0.5,col=rgb(0.9,0.6,0))
        }
        if(unsatisfied.flag==1){
          text(min(RR),130,paste0('Warning: We could not find classification of zones to satisfy your both GO and NOGO criterions at some sample sizes'),xpd=T,adj=0,cex=0.5,col=rgb(0.9,0,0))
          text(min(RR),120,paste0('Please check and modify your desicion rule or the range of sample size!'),xpd=T,adj=0,cex=0.5,col=rgb(0.9,0,0))

        }

      }
      if(direction=='Less'){
        if(overlap.flag==0&unsatisfied.flag==0){
          text(min(RR),130,bquote(GO~symbol("\336")~Observed~hazard~ratio~symbol("\243")~.(round(temp$cutoff[1],digits=3))),xpd=T,adj=0,cex=0.8,col=rgb(0, 0.7,0))
          text(min(RR),120,bquote(NOGO~symbol("\336")~Observed~hazard~ratio~symbol("\263")~.(round(temp$cutoff[2],digits=3))),xpd=T,adj=0,cex=0.8,col=rgb(0.9, 0,0))
          #text(min(RR),110,bquote(Inconclusive~symbol("\336")~anything~between~GO~and~NOGO),xpd=T,adj=0,cex=0.8,col=rgb(0.9,0.6,0))
        }

        if(overlap.flag==1&unsatisfied.flag==0&overlap.option=="GO"){
          text(min(RR),130,bquote(GO~symbol("\336")~Observed~hazard~ratio~symbol("\243")~.(round(temp$cutoff[1],digits=3))),xpd=T,adj=0,cex=0.8,col=rgb(0, 0.7,0))
          text(min(RR),120,bquote(NOGO~symbol("\336")~Observed~hazard~ratio~symbol("\263")~.(round(temp$cutoff[1],digits=3))),xpd=T,adj=0,cex=0.8,col=rgb(0.9, 0,0))
          #text(min(RR),115,paste0('Warning: GO and NOGO Zones are overlaped at some sample sizes'),xpd=T,adj=0,cex=0.5,col=rgb(0.9,0.6,0))
          #text(min(RR),110,paste0('The zones are classfied by cutoff of ', overlap.option),xpd=T,adj=0,cex=0.5,col=rgb(0.9,0.6,0))
        }
        if(overlap.flag==1&unsatisfied.flag==0&overlap.option=="NOGO"){
          text(min(RR),130,bquote(GO~symbol("\336")~Observed~hazard~ratio~symbol("\243")~.(round(temp$cutoff[2],digits=3))),xpd=T,adj=0,cex=0.8,col=rgb(0, 0.7,0))
          text(min(RR),120,bquote(NOGO~symbol("\336")~Observed~hazard~ratio~symbol("\263")~.(round(temp$cutoff[2],digits=3))),xpd=T,adj=0,cex=0.8,col=rgb(0.9, 0,0))
          #text(min(RR),115,paste0('Warning: GO and NOGO Zones are overlaped at some sample sizes'),xpd=T,adj=0,cex=0.5,col=rgb(0.9,0.6,0))
          #text(min(RR),110,paste0('The zones are classfied by cutoff of ', overlap.option),xpd=T,adj=0,cex=0.5,col=rgb(0.9,0.6,0))
        }
        if(unsatisfied.flag==1){
          text(min(RR),130,paste0('Warning: We could not find classification of zones to satisfy your both GO and NOGO criterions at some sample sizes'),xpd=T,adj=0,cex=0.5,col=rgb(0.9,0,0))
          text(min(RR),120,paste0('Please check and modify your desicion rule or the range of sample size!'),xpd=T,adj=0,cex=0.5,col=rgb(0.9,0,0))

        }

      }

    }
  }

######################################################################################

  return(list(go_prob=go_prob,nogo_prob=nogo_prob,
              inconclusive_prob=inconclusive_prob,
              overlap.flag=overlap.flag,overlap.option=overlap.option,
              unsatisfied.flag=unsatisfied.flag,cutoff=temp$cutoff,true_cutoff=c(true_go_cutoff,true_nogo_cutoff))
  )

}


Vary_SS_TAS_Survival_Prob<-function(eventinput=TRUE,
                                    neventsmin=70,neventsmax=350,
                                    npatientsmin=100,npatientsmax=500,
                                    maturity=0.7,
                                    alloc.ratio=1,
                                    CT1.go=0.25,
                                    false.go.CT1=TRUE,FGR.CT1=0.25,
                                    CT1.nogo=0.25,
                                    false.nogo.CT1=TRUE,FNGR.CT1=0.25,
                                    CT2.go=0.3,
                                    false.go.CT2=TRUE, FGR.CT2=0.5,
                                    CT2.nogo=0.3,
                                    false.nogo.CT2=TRUE, FNGR.CT2=0.5,
                                    RR=c(0.25),overlap.option='GO',
                                    plot.cutoff=TRUE,plot.prob=TRUE,
                                    method='Bayesian',direction='Greater',
                                    prior.mean=1,prior.sd=1,
                                    logic.go='and',logic.nogo='or'){
  if(eventinput==TRUE){nmin=neventsmin
  nmax=neventsmax
  nseq=unique(round(c(seq(nmin,(nmax+nmin)/2,length=6)[-6],seq((nmax+nmin)/2,nmax,length=6))))
  n1=nseq*alloc.ratio/(1+alloc.ratio)
  n2=nseq*1/(1+alloc.ratio)
  sdseq=sqrt(1/n1+1/n2)
  xlabtitle=paste0('Number of events')}else{
    nmin=npatientsmin
    nmax=npatientsmax
    nseq=seq(nmin,nmax,by=1)
    n1=nseq*maturity*alloc.ratio/(1+alloc.ratio)
    n2=nseq*maturity*1/(1+alloc.ratio)
    sdseq=sqrt(1/n1+1/n2)
    xlabtitle=paste0('Number of patients (Maturity= ',maturity,')')
  }
  go_prob<-matrix(NA,ncol=length(nseq),nrow=length(RR))
  nogo_prob<-matrix(NA,ncol=length(nseq),nrow=length(RR))
  inconclusive_prob<-matrix(NA,ncol=length(nseq),nrow=length(RR))
  go_cutoff<-rep(NA,length(nseq))
  nogo_cutoff<-rep(NA,length(nseq))
  index=1
  n_unsatisfied=NA
  n_overlap=NA
  for(i in nseq){
    temp=TAS_Survival_Cutoff(eventinput = eventinput,
                             npatients = i,maturity = maturity,
                             nevents=i,alloc.ratio = alloc.ratio,
                             CT1.go=CT1.go,
                             false.go.CT1 = false.go.CT1, FGR.CT1=FGR.CT1,
                             CT1.nogo=CT1.nogo,
                             false.nogo.CT1 = false.nogo.CT1, FNGR.CT1 = FNGR.CT1,
                             CT2.go=CT2.go,
                             false.go.CT2 = false.go.CT2,FGR.CT2=FGR.CT2,
                             CT2.nogo=CT2.nogo,
                             false.nogo.CT2 = false.nogo.CT2, FNGR.CT2 = FNGR.CT2,
                             method=method,direction=direction,
                             prior.mean=prior.mean,prior.sd=prior.sd,
                             logic.go=logic.go,logic.nogo=logic.nogo)
    if(direction=='Less'){

      if(all(temp$flag==0)&temp$overlap==0){
        go_prob[,index]=pnorm(log(temp$cutoff[1]),mean=log(RR),sd=sdseq[index])
        nogo_prob[,index]=1-pnorm(log(temp$cutoff[2]),mean=log(RR),sd=sdseq[index])
        inconclusive_prob[,index]=1-go_prob[,index]-nogo_prob[,index]
        go_cutoff[index]=temp$cutoff[1]
        nogo_cutoff[index]=temp$cutoff[2]
      }
      if(all(temp$flag==0)&temp$overlap==1&overlap.option=="GO"){
        go_prob[,index]=pnorm(log(temp$cutoff[1]),mean=log(RR),sd=sdseq[index])
        nogo_prob[,index]=1-pnorm(log(temp$cutoff[1]),mean=log(RR),sd=sdseq[index])
        inconclusive_prob[,index]=0
        go_cutoff[index]=temp$cutoff[1]
        nogo_cutoff[index]=temp$cutoff[2]
        n_overlap=c(n_overlap,i)
      }
      if(all(temp$flag==0)&temp$overlap==1&overlap.option=="NOGO"){
        go_prob[,index]=pnorm(log(temp$cutoff[2]),mean=log(RR),sd=sdseq[index])
        nogo_prob[,index]=1-pnorm(log(temp$cutoff[2]),mean=log(RR),sd=sdseq[index])
        inconclusive_prob[,index]=0
        go_cutoff[index]=temp$cutoff[1]
        nogo_cutoff[index]=temp$cutoff[2]
        n_overlap=c(n_overlap,i)
      }
      if(all(temp$flag==0)==FALSE){
        n_unsatisfied=c(n_unsatisfied,i)
      }

    }


    if(direction=='Greater'){

      if(all(temp$flag==0)&temp$overlap==0){
        go_prob[,index]=1-pnorm(log(temp$cutoff[1]),mean=log(RR),sd=sdseq[index])
        nogo_prob[,index]=pnorm(log(temp$cutoff[2]),mean=log(RR),sd=sdseq[index])
        inconclusive_prob[,index]=1-go_prob[,index]-nogo_prob[,index]
        go_cutoff[index]=temp$cutoff[1]
        nogo_cutoff[index]=temp$cutoff[2]
      }
      if(all(temp$flag==0)&temp$overlap==1&overlap.option=="GO"){
        go_prob[,index]=1-pnorm(log(temp$cutoff[1]),mean=log(RR),sd=sdseq[index])
        nogo_prob[,index]=pnorm(log(temp$cutoff[1]),mean=log(RR),sd=sdseq[index])
        inconclusive_prob[,index]=0
        go_cutoff[index]=temp$cutoff[1]
        nogo_cutoff[index]=temp$cutoff[2]
        n_overlap=c(n_overlap,i)
      }
      if(all(temp$flag==0)&temp$overlap==1&overlap.option=="NOGO"){
        go_prob[,index]=1-pnorm(log(temp$cutoff[2]),mean=log(RR),sd=sdseq[index])
        nogo_prob[,index]=pnorm(log(temp$cutoff[2]),mean=log(RR),sd=sdseq[index])
        inconclusive_prob[,index]=0
        go_cutoff[index]=temp$cutoff[1]
        nogo_cutoff[index]=temp$cutoff[2]
        n_overlap=c(n_overlap,i)
      }
      if(all(temp$flag==0)==FALSE){
        n_unsatisfied=c(n_unsatisfied,i)
      }

    }


    index=index+1
  }
  n_overlap=n_overlap[-1]
  n_unsatisfied=n_unsatisfied[-1]
  ####plot figure
  if(plot.prob==TRUE){
    for(j in 1:length(RR)){
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
      plot(delta,cum_p_go,xlab=xlabtitle,
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
      text((nmin+nmax)/2,120,paste0("True HR=",round(RR[j],digits=3)),xpd=T,adj=0.5,cex=0.8)
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
    p.obs_go=go_cutoff
    p.obs_nogo=nogo_cutoff
    P_go=p.obs_go
    P_nogo=p.obs_nogo
    par(mfrow=c(1,1),pty="m",bty="l",cex=1.4)

    ylim_max=max(c(P_go,P_nogo,CT1.go,CT2.go,CT1.nogo,CT2.nogo),na.rm=TRUE)+0.1
    ylim_min=min(c(P_go,P_nogo,CT1.go,CT2.go,CT1.nogo,CT2.nogo),na.rm=TRUE)-0.1
    plot(NA,NA,xlab=xlabtitle,ylab="Estimated hazard ratio",xlim=c(max(0,range(nseq)[1]-10),range(nseq)[2]+10),ylim=c(ylim_min,ylim_max),type="n",axes=F,col=rgb(1,0,0),lty=1,lwd=2)
    axis(1, at=nseq, labels=T)
    axis(2, at=round(c(seq(ylim_min,ylim_max,0.1),CT1.go,CT2.go,CT1.nogo,CT2.nogo),2),labels=T)
    box()
    lines(nseq,P_nogo,col=rgb(0.9,0,0),lwd=2)
    lines(nseq,P_go,col=rgb(0,0.7,0),lwd=2)
    legend('bottomright',legend=c("Cut off of GO","Cut off of NOGO"),
           col=c(rgb(0,0.7,0),rgb(0.9,0,0)),
           lwd=c(2,2),
           lty=c(1,1),cex=0.5)
    #text((nmin+nmax)/2,ylim_max+0.075,paste0("True HR=",round(RR[j],digits=3)),xpd=T,adj=0.5,cex=0.8)
    if(length(n_overlap)!=0){
      text(nmin,ylim_max+0.05,paste0('Warning: GO and NOGO zones are overlaped at some sample sizes'),xpd=T,adj=0,cex=0.5,col=rgb(0.9,0.6,0))
      #text(nmin,110,paste(n_overlap),xpd=T,adj=0,cex=0.5,col=rgb(0.9,0.6,0))
      #print(n_overlap)
      text(nmin,ylim_max+0.025,paste0('The zones are classfied by cutoff of ', overlap.option),xpd=T,adj=0,cex=0.5,col=rgb(0.9,0.6,0))
    }
    if(length(n_unsatisfied)!=0){
      text(nmin,ylim_max+0.05,paste0('Warning: We could not find classification of zones to satisfy your both GO and NOGO criterions at some sample sizes'),xpd=T,adj=0,cex=0.5,col=rgb(0.9,0,0))
      text(nmin,ylim_max+0.025,paste0('Please check and modify your desicion rule or the range of sample size!'),xpd=T,adj=0,cex=0.5,col=rgb(0.9,0,0))

    }


  }
  return(list(n_unsatisfied,n_overlap))


}



Interim_TAS<-function(interim_n=c(50,100,150),num_interim=3,
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
                      logic.go=c('and','or','and'),
                      logic.nogo=c('and','or','or'),
                      seed.num=369,nsim_IA=10000,stop.criterion=10^-3,
                      eventinput=TRUE,
                      maturity=0.7,
                      alloc.ratio = 1,
                      prior.mean = 1,prior.sd=1,
                      mean=c(0.25),
                      nsim=10000){

  if(eventinput==FALSE){
    interim_n=floor(interim_n*maturity)
  }
  interim_n=sort(interim_n)
  go_cutoff<-rep(NA,num_interim)
  nogo_cutoff<-rep(NA,num_interim)
  true_go_cutoff<-rep(NA,num_interim)
  true_nogo_cutoff<-rep(NA,num_interim)
  overlap<-rep(NA,num_interim)
  for(i in 1:num_interim){
    temp<-TAS_Survival_Cutoff(eventinput=TRUE,nevents=interim_n[i],npatients=interim_n[i],
                              maturity=maturity,
                              alloc.ratio = alloc.ratio,
                              prior.mean = prior.mean,prior.sd = prior.sd,
                              CT1.go=CT1.go[i],
                              false.go.CT1=false.go.CT1[i],FGR.CT1=FGR.CT1[i],
                              CT1.nogo=CT1.nogo[i],
                              false.nogo.CT1=false.nogo.CT1[i],FNGR.CT1=FNGR.CT1[i],
                              CT2.go=CT2.go[i],
                              false.go.CT2=false.go.CT2[i], FGR.CT2=FGR.CT2[i],
                              CT2.nogo=CT1.nogo[i],
                              false.nogo.CT2=false.nogo.CT2[i], FNGR.CT2=FNGR.CT2[i],
                              method=method,direction=direction[i],
                              logic.go=logic.go[i],logic.nogo=logic.nogo[i])

    go_cutoff[i]<-temp$cutoff[1]
    nogo_cutoff[i]<-temp$cutoff[2]
    overlap[i]<-temp$overlap
    ###TAS####
    true_go_cutoff[i]<-ifelse(overlap[i]==0,go_cutoff[i],(overlap.option[i]=='GO')*go_cutoff[i]+(overlap.option[i]=='NOGO')*(nogo_cutoff[i]))
    true_nogo_cutoff[i]<-ifelse(overlap[i]==0,nogo_cutoff[i],(overlap.option[i]=='GO')*(go_cutoff[i])+(overlap.option[i]=='NOGO')*(nogo_cutoff[i]))
    ####
  }

  temptable=c()
  for(meanindex in 1:length(mean)){
    set.seed(seed.num)
    ###TAS####
    precision1=alloc.ratio/(1+alloc.ratio)
    precision2=1-precision1
    sd=sqrt(1/precision1+1/precision2)
    sim_data<-matrix(NA,nrow=nsim_IA,ncol=max(interim_n))
    for(k in 1:max(interim_n)){
      sim_data[,k]<- rnorm(nsim_IA,mean=log(mean[meanindex]),sd=sd)
    }

    cum_sim_data_temp<-apply(sim_data,1,cumsum)
    cum_sim_data<-t(exp(cum_sim_data_temp[interim_n,]/interim_n))
    #####
    go_matrix<-matrix(NA,nrow=nsim_IA,ncol=num_interim)
    nogo_matrix<-matrix(NA,nrow=nsim_IA,ncol=num_interim)
    inconclusive_matrix<-matrix(NA,nrow=nsim_IA,ncol=num_interim)
    table<-matrix(NA,ncol=num_interim+1,nrow=10)
    IA_go_matrix<-matrix(NA,nrow=nsim_IA,ncol=num_interim) ###whether continue to next stage
    for(j in 1:num_interim ){
      if(direction[j]=='Greater'){
        go_matrix[,j]<-cum_sim_data[,j]>=true_go_cutoff[j]
        nogo_matrix[,j]<-cum_sim_data[,j]<=true_nogo_cutoff[j]
        inconclusive_matrix[,j]<-rep(1,nsim_IA)-go_matrix[,j]-nogo_matrix[,j]
      }
    }

    for(j in 1:num_interim ){
      if(direction[j]=='Less'){
        go_matrix[,j]<-cum_sim_data[,j]<=true_go_cutoff[j]
        nogo_matrix[,j]<-cum_sim_data[,j]>=true_nogo_cutoff[j]
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
    tablerowname<-c('Sample size','True HR','Task','Success','To next interim/final or inconclusive',
                    'Stop','Superority/Futility zone','Warning',
                    'Cut off for NOGO rule',
                    'Cut off for GO rule')
    table<-cbind(tablerowname,rep(meanindex,10),table)
    colnames(table)<-c(" ",'Setting',tablecolname)
    temptable=rbind(temptable,table)
  }
  return(temptable)
}
#Interim_TAS(interim_n = c(66,131,197),CT1.go=c(1,0.6,1),FGR.CT1=c(1-0.71,0.653,0.0238),FNGR.CT1=c(0.707,1-0.653,1-0.0238),CT2.go=c(NA,NA,NA),method='Frequentist',direction=c(rep('Less',2),'Greater'),task=c('Futility','Superiority','Superiority'),mean=c(0.25,2))
#
# Fix_SS_TAS_Survival_Prob(eventinput=FALSE,nevents=170,
#                                    npatients=180,maturity=0.4,
#                                    alloc.ratio=1,
#                                    CT1.go=0.65,
#                                    false.go.CT1=TRUE,FGR.CT1=0.7,
#                                    CT1.nogo=0.65,
#                                    false.nogo.CT1=TRUE,FNGR.CT1=0.1,
#                                    CT2.go=10/12,
#                                    false.go.CT2=FALSE, FGR.CT2=0.2,
#                                    CT2.nogo=10/12,
#                                    false.nogo.CT2=FALSE, FNGR.CT2=0.2,
#                                    RR=c(0.933,0.483),overlap.option='GO',plot.figure=TRUE,
#                                    method='Frequentist',direction='Less',prior.mean=1,prior.sd=1,
#                                    logic.go='and',logic.nogo='or')

