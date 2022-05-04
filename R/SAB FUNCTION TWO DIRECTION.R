library(shiny)
####Single arm binomial functions####
SAB_Bin_Cutoff<-function(n=30,a=1/3,b=1,
                         CT1.go=0.25,
                         false.go.CT1=TRUE,FGR.CT1=0.25,
                         CT1.nogo=0.25,
                         false.nogo.CT1=TRUE,FNGR.CT1=0.25,
                         CT2.go=0.3,
                         false.go.CT2=TRUE, FGR.CT2=0.5,
                         CT2.nogo=0.3,
                         false.nogo.CT2=TRUE, FNGR.CT2=0.5,
                         method='Bayesian',direction='Greater',logic.go='and',logic.nogo='and'){

  r=0:n
  index1.go=NA
  index2.go=NA
  index1.nogo=NA
  index2.nogo=NA
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



  if(direction=="Greater"){
    if(method=='Bayesian'){

      if(false.go.CT1==TRUE){
        est1.go=pbeta(CT1.go,a+r,b+n-r)
        if(all(est1.go>FGR.CT1)){
          flag[1]=1
        }else{

          index1.go=which.max(est1.go<=FGR.CT1)
        }
      }
      if(false.go.CT2==TRUE){
        est2.go=pbeta(CT2.go,a+r,b+n-r)
        if(all(est2.go>FGR.CT2)){
          flag[2]=1
        }else{
          index2.go=which.max(est2.go<=FGR.CT2)
        }
      }

      if(false.nogo.CT1==TRUE){
        est1.nogo=1-pbeta(CT1.nogo,a+r,b+n-r)
        if(all(est1.nogo>FNGR.CT1)){
          flag[3]=1
        }else{
          if(all(est1.nogo<=FNGR.CT1)){index1.nogo=n+1}else{
            index1.nogo=which.max(est1.nogo>FNGR.CT1)}
        }

      }
      if(false.nogo.CT2==TRUE){
        est2.nogo=1-pbeta(CT2.nogo,a+r,b+n-r)
        if(all(est2.nogo>FNGR.CT2)){
          flag[4]=1
        }else{
          if(all(est2.nogo<=FNGR.CT2)){index2.nogo=n+1}else{
            index2.nogo=which.max(est2.nogo>FNGR.CT2)}
        }

      }

    }
    if(method=='Frequentist'){
      if(false.go.CT1==TRUE){
        est1.go<-c()
        for(i in 1:length(r)){
          est1.go=c(est1.go,binom.test(r[i],n,p=CT1.go,alternative='greater',conf.level = 1-FGR.CT1)$conf.int[1])
        }
        if(all(est1.go<= CT1.go)){
          flag[1]=1
        }else{
          index1.go=which.max(est1.go>CT1.go)
        }
      }
      if(false.go.CT2==TRUE){
        est2.go<-c()
        for(i in 1:length(r)){
          est2.go=c(est2.go,binom.test(r[i],n,p=CT2.go,alternative='greater',conf.level = 1-FGR.CT2)$conf.int[1])
        }
        if(all(est2.go<= CT2.go)){
          flag[2]=1
        }else{
          index2.go=which.max(est2.go>CT2.go)
        }
      }

      if(false.nogo.CT1==TRUE){
        est1.nogo<-c()
        for(i in 1:length(r)){
          est1.nogo=c(est1.nogo,binom.test(r[i],n,p=CT1.nogo,alternative='greater',conf.level = FNGR.CT1)$conf.int[1])
        }
        if(all(est1.nogo>CT1.nogo)){
          flag[3]=1
        }else{
          if(all(est1.nogo<CT1.nogo)){index1.nogo=n+1}else{
            index1.nogo=which.max(est1.nogo>=CT1.nogo)}
        }

      }
      if(false.nogo.CT2==TRUE){
        est2.nogo<-c()
        for(i in 1:length(r)){
          est2.nogo=c(est2.nogo,binom.test(r[i],n,p=CT2.nogo,alternative='greater',conf.level = FNGR.CT2)$conf.int[1])
        }
        if(all(est2.nogo>CT2.nogo)){
          flag[4]=1
        }else{
          if(all(est2.nogo<CT2.nogo)){index1.nogo=n+1}else{
            index2.nogo=which.max(est2.nogo>=CT2.nogo)}
        }

      }

    }

    if(all(flag==0)==TRUE){
      if(any(is.na(c(index1.go,index2.go)))){logic.go='and'}
      if(any(is.na(c(index1.nogo,index2.nogo)))){logic.nogo='and'}
      if(logic.go=='and'){
        go_cutoff=r[max(c(index1.go,index2.go),na.rm=TRUE)]
      }
      if(logic.go=='or'){
        go_cutoff=r[min(c(index1.go,index2.go),na.rm=TRUE)]
      }
      if(logic.nogo=='and')
      {
        nogo_cutoff=r[min(c(index1.nogo,index2.nogo),na.rm=TRUE)]
      }
      if(logic.nogo=='or')
      {
        nogo_cutoff=r[max(c(index1.nogo,index2.nogo),na.rm=TRUE)]
      }
      if(go_cutoff>=nogo_cutoff){return(list(cutoff=c(go_cutoff,nogo_cutoff),flag=flag,overlap=overlap.flag))}
      else{
        overlap.flag=1
        return(list(cutoff=c(go_cutoff,nogo_cutoff),flag=flag,overlap=overlap.flag))}
    }
    if(all(flag==0)==FALSE){return(list(cutoff=c(NA,NA),flag=flag,overlap=overlap.flag))}
  }


  if(direction=="Less"){
    if(method=='Bayesian'){

      if(false.go.CT1==TRUE){
        est1.go=pbeta(CT1.go,a+r,b+n-r,lower.tail = FALSE)
        if(all(est1.go>FGR.CT1)){
          flag[1]=1
        }else{
          if(all(est1.go<=FGR.CT1)){index1.go=n+1}else{
            index1.go=which.max(est1.go>FGR.CT1)-1}
        }
      }
      if(false.go.CT2==TRUE){
        est2.go=pbeta(CT2.go,a+r,b+n-r,lower.tail = FALSE)
        if(all(est2.go>FGR.CT2)){
          flag[2]=1
        }else{
          if(all(est2.go<=FGR.CT2)){index2.go=n+1}else{
            index2.go=which.max(est2.go>FGR.CT2)-1}
        }
      }

      if(false.nogo.CT1==TRUE){
        est1.nogo=1-pbeta(CT1.nogo,a+r,b+n-r,lower.tail = FALSE)
        if(all(est1.nogo>FNGR.CT1)){
          flag[3]=1
        }else{
          index1.nogo=which.max(est1.nogo<=FNGR.CT1)-1
        }

      }

      if(false.nogo.CT2==TRUE){
        est2.nogo=1-pbeta(CT2.nogo,a+r,b+n-r,lower.tail = FALSE)
        if(all(est2.nogo>FNGR.CT2)){
          flag[4]=1
        }else{
          index2.nogo=which.max(est2.nogo<=FNGR.CT2)-1
        }

      }

    }
    if(method=='Frequentist'){

      if(false.go.CT1==TRUE){
        est1.go<-c()
        for(i in 1:length(r)){
          est1.go=c(est1.go,binom.test(r[i],n,p=CT1.go,alternative='less',conf.level = 1-FGR.CT1)$conf.int[2])
        }
        if(all(est1.go>=CT1.go)){
          flag[1]=1
        }else{
          if(all(est1.go<CT1.go)){index1.go=n+1}else{
            index1.go=which.max(est1.go>=CT1.go)-1}
        }

      }
      if(false.go.CT2==TRUE){
        est2.go<-c()
        for(i in 1:length(r)){
          est2.go=c(est2.go,binom.test(r[i],n,p=CT2.go,alternative='less',conf.level = 1-FGR.CT2)$conf.int[2])
        }
        if(all(est2.go>=CT2.go)){
          flag[2]=1
        }else{
          if(all(est2.go<CT2.go)){index1.go=n+1}else{
            index2.go=which.max(est2.go>=CT2.go)-1}
        }

      }
      if(false.nogo.CT1==TRUE){
        est1.nogo<-c()
        for(i in 1:length(r)){
          est1.nogo=c(est1.nogo,binom.test(r[i],n,p=CT1.nogo,alternative='less',conf.level = FNGR.CT1)$conf.int[2])
        }
        if(all(est1.nogo> CT1.nogo)){
          flag[3]=1
        }else{
          index1.nogo=which.max(est1.nogo>=CT1.nogo)-1
        }
      }
      if(false.nogo.CT2==TRUE){
        est2.nogo<-c()
        for(i in 1:length(r)){
          est2.nogo=c(est2.nogo,binom.test(r[i],n,p=CT2.nogo,alternative='less',conf.level = FNGR.CT2)$conf.int[2])
        }
        if(all(est2.nogo> CT2.nogo)){
          flag[4]=1
        }else{
          index2.nogo=which.max(est2.nogo>=CT2.nogo)-1
        }
      }
    }

    if(all(flag==0)==TRUE){
      if(any(is.na(c(index1.go,index2.go)))){logic.go='and'}
      if(any(is.na(c(index1.nogo,index2.nogo)))){logic.nogo='and'}
      if(logic.go=='and'){
        go_cutoff=r[min(c(index1.go,index2.go),na.rm=TRUE)]
      }
      if(logic.go=='or'){
        go_cutoff=r[max(c(index1.go,index2.go),na.rm=TRUE)]
      }
      if(logic.nogo=='and')
      {
        nogo_cutoff=r[max(c(index1.nogo,index2.nogo),na.rm=TRUE)]
      }
      if(logic.nogo=='or')
      {
        nogo_cutoff=r[min(c(index1.nogo,index2.nogo),na.rm=TRUE)]
      }

      if(go_cutoff<=nogo_cutoff){
        return(list(cutoff=c(go_cutoff,nogo_cutoff),flag=flag,overlap=overlap.flag))}
      else{
        overlap.flag=1
        return(list(cutoff=c(go_cutoff,nogo_cutoff),flag=flag,overlap=overlap.flag))}
    }
    if(all(flag==0)==FALSE){return(list(cutoff=c(NA,NA),flag=flag,overlap=overlap.flag))}
  }
}

########################################
Fix_SS_SAB_Bin_Prob<-function(n=100,a=1/3,b=1,
                              CT1.go=0.25,
                              false.go.CT1=TRUE,FGR.CT1=0.25,
                              CT1.nogo=0.25,
                              false.nogo.CT1=TRUE,FNGR.CT1=0.25,
                              CT2.go=0.3,
                              false.go.CT2=TRUE, FGR.CT2=0.5,
                              CT2.nogo=0.3,
                              false.nogo.CT2=TRUE, FNGR.CT2=0.5,
                              RR=c(0.0,0.9),overlap.option='GO',
                              plot.figure=TRUE,
                              method='Bayesian',direction='Greater',logic.go='and',logic.nogo='and'){
  RRseq=seq(min(RR),max(RR),0.05)
  go_prob<-rep(NA,length(RR))
  nogo_prob<-rep(NA,length(RR))
  inconclusive_prob<-rep(NA,length(RR))
  go_prob_plot<-rep(NA,length(RRseq))
  nogo_prob_plot<-rep(NA,length(RRseq))
  inconclusive_prob_plot<-rep(NA,length(RRseq))
  index=1
  unsatisfied.flag=0
  overlap.flag=0

  temp=SAB_Bin_Cutoff(n=n,a=a,b=b,
                      CT1.go=CT1.go,
                      false.go.CT1 = false.go.CT1, FGR.CT1=FGR.CT1,
                      CT1.nogo=CT1.nogo,
                      false.nogo.CT1 = false.nogo.CT1, FNGR.CT1 = FNGR.CT1,
                      CT2.go=CT2.go,
                      false.go.CT2 = false.go.CT2, FGR.CT2=FGR.CT2,
                      CT2.nogo=CT2.nogo,
                      false.nogo.CT2 = false.nogo.CT2, FNGR.CT2 = FNGR.CT2,
                      method=method,direction=direction,logic.go = logic.go,logic.nogo=logic.nogo)

  ###SAB####
  if(direction=='Greater'){
    true_go_cutoff<-ifelse(temp$overlap==0,temp$cutoff[1],(overlap.option=='GO')*temp$cutoff[1]+(overlap.option=='NOGO')*(temp$cutoff[2]))
    true_nogo_cutoff<-ifelse(temp$overlap==0,temp$cutoff[2],(overlap.option=='GO')*(temp$cutoff[1])+(overlap.option=='NOGO')*(temp$cutoff[2]))
  }
  if(direction=='Less'){
    true_go_cutoff<-ifelse(temp$overlap==0,temp$cutoff[1],(overlap.option=='GO')*temp$cutoff[1]+(overlap.option=='NOGO')*(temp$cutoff[2]))
    true_nogo_cutoff<-ifelse(temp$overlap==0,temp$cutoff[2],(overlap.option=='GO')*(temp$cutoff[1])+(overlap.option=='NOGO')*(temp$cutoff[2]))
  }


  if(direction=='Greater'){
    if(all(temp$flag==0)&temp$overlap==0){
      go_prob=1-pbinom(temp$cutoff[1]-1,n,RR)
      nogo_prob=pbinom(temp$cutoff[2]-1,n,RR)
      inconclusive_prob=1-go_prob-nogo_prob
      go_prob_plot=1-pbinom(temp$cutoff[1]-1,n,RRseq)
      nogo_prob_plot=pbinom(temp$cutoff[2]-1,n,RRseq)
      inconclusive_prob_plot=1-go_prob_plot-nogo_prob_plot
    }
    if(all(temp$flag==0)&temp$overlap==1&overlap.option=="GO"){
      go_prob=1-pbinom(temp$cutoff[1]-1,n,RR)
      nogo_prob=1-go_prob
      inconclusive_prob=1-go_prob-nogo_prob

      go_prob_plot=1-pbinom(temp$cutoff[1]-1,n,RRseq)
      nogo_prob_plot=1-go_prob_plot
      inconclusive_prob_plot=1-go_prob_plot-nogo_prob_plot

      overlap.flag=1
    }
    if(all(temp$flag==0)&temp$overlap==1&overlap.option=="NOGO"){
      go_prob=1-pbinom(temp$cutoff[2]-1,n,RR)
      nogo_prob=1-go_prob
      inconclusive_prob=1-go_prob-nogo_prob

      nogo_prob_plot=pbinom(temp$cutoff[2]-1,n,RRseq)
      go_prob_plot=1-nogo_prob_plot
      inconclusive_prob_plot=1-go_prob_plot-nogo_prob_plot

      overlap.flag=1
    }
    if(all(temp$flag==0)==FALSE){
      unsatisfied.flag=1
      nogo_prob=rep(NA,length(RR))
      go_prob=rep(NA,length(RR))
      inconclusive_prob=rep(NA,length(RR))
    }

  }

  if(direction=='Less'){
    if(all(temp$flag==0)&temp$overlap==0){
      go_prob=pbinom(temp$cutoff[1],n,RR)
      nogo_prob=1-pbinom(temp$cutoff[2],n,RR)
      inconclusive_prob=1-go_prob-nogo_prob

      go_prob_plot=pbinom(temp$cutoff[1],n,RRseq)
      nogo_prob_plot=1-pbinom(temp$cutoff[2],n,RRseq)
      inconclusive_prob_plot=1-go_prob_plot-nogo_prob_plot
    }
    if(all(temp$flag==0)&temp$overlap==1&overlap.option=="GO"){
      go_prob=pbinom(temp$cutoff[1],n,RR)
      nogo_prob=1-go_prob
      inconclusive_prob=1-go_prob-nogo_prob

      go_prob_plot=pbinom(temp$cutoff[1],n,RRseq)
      nogo_prob_plot=1-go_prob_plot
      inconclusive_prob_plot=1-go_prob_plot-nogo_prob_plot

      overlap.flag=1
    }
    if(all(temp$flag==0)&temp$overlap==1&overlap.option=="NOGO"){
      go_prob=pbinom(temp$cutoff[2],n,RR)
      nogo_prob=1-go_prob
      inconclusive_prob=1-go_prob-nogo_prob

      go_prob_plot=pbinom(temp$cutoff[2],n,RRseq)
      nogo_prob_plot=1-go_prob_plot
      inconclusive_prob_plot=1-go_prob-nogo_prob_plot

      overlap.flag=1
    }
    if(all(temp$flag==0)==FALSE){
      unsatisfied.flag=1
      nogo_prob=rep(NA,length(RR))
      go_prob=rep(NA,length(RR))
      inconclusive_prob=rep(NA,length(RR))
    }

  }

###########################@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

##########################################################################################
  if(plot.figure==TRUE){

    delta=RRseq
    p_nogo=nogo_prob_plot*100
    p_grey=inconclusive_prob_plot*100
    p_go=go_prob_plot*100

    if(true_go_cutoff!=true_nogo_cutoff){

      par(mfrow=c(1,1),pty="m",bty="l",cex=1.4)

      plot(delta,p_nogo,xlab=expression(paste("True ",RR,sep="")),
           ylab="Prob of GO/ NOGO/ Inconclusive (%)",xlim=range(delta),
           ylim=c(0,100),type="n",axes=F)

      axis(1, at=RRseq, labels=T)
      axis(2, at=seq(0,100,10),labels=T)
      #text(rep(-0.6,6),seq(0,100,20),seq(0,100,20),adj=1,xpd=T)
      box()

      points(delta,p_nogo,type="b",pch=16,col=rgb(0.9,0,0),lwd=3,lty=1)
      points(delta,p_grey,type="b",pch=4,col=rgb(0.9,0.6,0),lwd=3,lty=6)
      points(delta,p_go,type="b",pch=2,col=rgb(0,0.7,0),lwd=3,lty=2)

      text(min(RR),140,bquote(n==~.(n)),xpd=T,adj=0,cex=0.8,col=rgb(0, 0,0))
      if(overlap.flag==0&unsatisfied.flag==0){
        if(direction=='Greater'){
          text(min(RR),130,bquote(GO~symbol("\336")~Number~of~responders~symbol("\263")~.(temp$cutoff[1])),xpd=T,adj=0,cex=0.8,col=rgb(0, 0.7,0))
          text(min(RR),120,bquote(NOGO~symbol("\336")~Number~of~responders~symbol("\74")~.(temp$cutoff[2])),xpd=T,adj=0,cex=0.8,col=rgb(0.9, 0,0))
          text(min(RR),110,bquote(Inconlusive~symbol("\336")~anything~between~GO~and~NOGO),xpd=T,adj=0,cex=0.8,col=rgb(0.9,0.6,0))
        }
        if(direction=='Less'){
          text(min(RR),130,bquote(GO~symbol("\336")~Number~of~responders~symbol("\243")~.(temp$cutoff[1])),xpd=T,adj=0,cex=0.8,col=rgb(0, 0.7,0))
          text(min(RR),120,bquote(NOGO~symbol("\336")~Number~of~responders~symbol("\76")~.(temp$cutoff[2])),xpd=T,adj=0,cex=0.8,col=rgb(0.9, 0,0))
          text(min(RR),110,bquote(Inconclusive~symbol("\336")~anything~between~GO~and~NOGO),xpd=T,adj=0,cex=0.8,col=rgb(0.9,0.6,0))
        }
      }
      if(overlap.flag==1&unsatisfied.flag==0&overlap.option=='GO'){
        if(direction=="Greater"){
          text(min(RR),130,bquote(GO~symbol("\336")~Number~of~responders~symbol("\263")~.(temp$cutoff[1])),xpd=T,adj=0,cex=0.8,col=rgb(0, 0.7,0))
          text(min(RR),120,bquote(NOGO~symbol("\336")~Number~of~responders~symbol("\74")~.(temp$cutoff[1])),xpd=T,adj=0,cex=0.8,col=rgb(0.9, 0,0))
          #text(min(RR),115,paste0('Warning: GO and NOGO zones are overlaped at some sample sizes'),xpd=T,adj=0,cex=0.5,col=rgb(0.9,0.6,0))
          #text(min(RR),110,paste0('The zones are classfied by cutoff of ', overlap.option),xpd=T,adj=0,cex=0.5,col=rgb(0.9,0.6,0))

        }
        if(direction=="Less"){
          text(min(RR),130,bquote(GO~symbol("\336")~Number~of~responders~symbol("\243")~.(temp$cutoff[1])),xpd=T,adj=0,cex=0.8,col=rgb(0, 0.7,0))
          text(min(RR),120,bquote(NOGO~symbol("\336")~Number~of~responders~symbol("\76")~.(temp$cutoff[1])),xpd=T,adj=0,cex=0.8,col=rgb(0.9, 0,0))
          #text(min(RR),115,paste0('Warning: GO and NOGO zones are overlaped at some sample sizes'),xpd=T,adj=0,cex=0.5,col=rgb(0.9,0.6,0))
          #text(min(RR),110,paste0('The zones are classfied by cutoff of ', overlap.option),xpd=T,adj=0,cex=0.5,col=rgb(0.9,0.6,0))

        }
      }

      if(overlap.flag==1&unsatisfied.flag==0&overlap.option=='NOGO'){
        if(direction=="Greater"){
          text(min(RR),130,bquote(GO~symbol("\336")~Number~of~responders~symbol("\263")~.(temp$cutoff[2])),xpd=T,adj=0,cex=0.8,col=rgb(0, 0.7,0))
          text(min(RR),120,bquote(NOGO~symbol("\336")~Number~of~responders~symbol("\74")~.(temp$cutoff[2])),xpd=T,adj=0,cex=0.8,col=rgb(0.9, 0,0))
          #text(min(RR),115,paste0('Warning: GO and NOGO zones are overlaped at some sample sizes'),xpd=T,adj=0,cex=0.5,col=rgb(0.9,0.6,0))
          #text(min(RR),110,paste0('The zones are classfied by cutoff of ', overlap.option),xpd=T,adj=0,cex=0.5,col=rgb(0.9,0.6,0))
        }
        if(direction=="Less"){
          text(min(RR),130,bquote(GO~symbol("\336")~Number~of~responders~symbol("\243")~.(temp$cutoff[2])),xpd=T,adj=0,cex=0.8,col=rgb(0, 0.7,0))
          text(min(RR),120,bquote(NOGO~symbol("\336")~Number~of~responders~symbol("\76")~.(temp$cutoff[2])),xpd=T,adj=0,cex=0.8,col=rgb(0.9, 0,0))
          #text(min(RR),115,paste0('Warning: GO and NOGO zones are overlaped at some sample sizes'),xpd=T,adj=0,cex=0.5,col=rgb(0.9,0.6,0))
          #text(min(RR),110,paste0('The zones are classfied by cutoff of ', overlap.option),xpd=T,adj=0,cex=0.5,col=rgb(0.9,0.6,0))
        }
      }
      if(unsatisfied.flag==1){
        text(min(RR),130,paste0('Warning: We could not find classification of zones to satisfy your both GO and NOGO criterions at some sample sizes'),xpd=T,adj=0,cex=0.5,col=rgb(0.9,0,0))
        text(min(RR),120,paste0('Please check and modify your desicion rule or the range of sample size!'),xpd=T,adj=0,cex=0.5,col=rgb(0.9,0,0))

      }
    }
  }
  ################### end of first if condition #######
  if(plot.figure==TRUE){

    delta=RRseq
    p_nogo=nogo_prob_plot*100
    p_grey=inconclusive_prob_plot*100
    p_go=go_prob_plot*100

    if (true_go_cutoff==true_nogo_cutoff){

      par(mfrow=c(1,1),pty="m",bty="l",cex=1.4)

      plot(delta,p_nogo,xlab=expression(paste("True ",RR,sep="")),
           ylab="Prob of GO/ NOGO/ Inconclusive (%)",xlim=range(delta),
           ylim=c(0,100),type="n",axes=F)

      axis(1, at=RRseq, labels=T)
      axis(2, at=seq(0,100,10),labels=T)
      #text(rep(-0.6,6),seq(0,100,20),seq(0,100,20),adj=1,xpd=T)
      box()

      points(delta,p_nogo,type="b",pch=16,col=rgb(0.9,0,0),lwd=3,lty=1)
      #points(delta,p_grey,type="b",pch=4,col=rgb(0.9,0.6,0),lwd=3,lty=6)
      points(delta,p_go,type="b",pch=2,col=rgb(0,0.7,0),lwd=3,lty=2)

      text(min(RR),140,bquote(n==~.(n)),xpd=T,adj=0,cex=0.8,col=rgb(0, 0,0))
      if(overlap.flag==0&unsatisfied.flag==0){
        if(direction=='Greater'){
          text(min(RR),130,bquote(GO~symbol("\336")~Number~of~responders~symbol("\263")~.(temp$cutoff[1])),xpd=T,adj=0,cex=0.8,col=rgb(0, 0.7,0))
          text(min(RR),120,bquote(NOGO~symbol("\336")~Number~of~responders~symbol("\74")~.(temp$cutoff[2])),xpd=T,adj=0,cex=0.8,col=rgb(0.9, 0,0))
          #text(min(RR),110,bquote(Inconlusive~symbol("\336")~anything~between~GO~and~NOGO),xpd=T,adj=0,cex=0.8,col=rgb(0.9,0.6,0))
        }
        if(direction=='Less'){
          text(min(RR),130,bquote(GO~symbol("\336")~Number~of~responders~symbol("\243")~.(temp$cutoff[1])),xpd=T,adj=0,cex=0.8,col=rgb(0, 0.7,0))
          text(min(RR),120,bquote(NOGO~symbol("\336")~Number~of~responders~symbol("\76")~.(temp$cutoff[2])),xpd=T,adj=0,cex=0.8,col=rgb(0.9, 0,0))
          #text(min(RR),110,bquote(Inconclusive~symbol("\336")~anything~between~GO~and~NOGO),xpd=T,adj=0,cex=0.8,col=rgb(0.9,0.6,0))
        }
      }
      if(overlap.flag==1&unsatisfied.flag==0&overlap.option=='GO'){
        if(direction=="Greater"){
          text(min(RR),130,bquote(GO~symbol("\336")~Number~of~responders~symbol("\263")~.(temp$cutoff[1])),xpd=T,adj=0,cex=0.8,col=rgb(0, 0.7,0))
          text(min(RR),120,bquote(NOGO~symbol("\336")~Number~of~responders~symbol("\74")~.(temp$cutoff[1])),xpd=T,adj=0,cex=0.8,col=rgb(0.9, 0,0))
          #text(min(RR),115,paste0('Warning: GO and NOGO zones are overlaped at some sample sizes'),xpd=T,adj=0,cex=0.5,col=rgb(0.9,0.6,0))
          #text(min(RR),110,paste0('The zones are classfied by cutoff of ', overlap.option),xpd=T,adj=0,cex=0.5,col=rgb(0.9,0.6,0))

        }
        if(direction=="Less"){
          text(min(RR),130,bquote(GO~symbol("\336")~Number~of~responders~symbol("\243")~.(temp$cutoff[1])),xpd=T,adj=0,cex=0.8,col=rgb(0, 0.7,0))
          text(min(RR),120,bquote(NOGO~symbol("\336")~Number~of~responders~symbol("\76")~.(temp$cutoff[1])),xpd=T,adj=0,cex=0.8,col=rgb(0.9, 0,0))
          #text(min(RR),115,paste0('Warning: GO and NOGO zones are overlaped at some sample sizes'),xpd=T,adj=0,cex=0.5,col=rgb(0.9,0.6,0))
          #text(min(RR),110,paste0('The zones are classfied by cutoff of ', overlap.option),xpd=T,adj=0,cex=0.5,col=rgb(0.9,0.6,0))

        }
      }

      if(overlap.flag==1&unsatisfied.flag==0&overlap.option=='NOGO'){
        if(direction=="Greater"){
          text(min(RR),130,bquote(GO~symbol("\336")~Number~of~responders~symbol("\263")~.(temp$cutoff[2])),xpd=T,adj=0,cex=0.8,col=rgb(0, 0.7,0))
          text(min(RR),120,bquote(NOGO~symbol("\336")~Number~of~responders~symbol("\74")~.(temp$cutoff[2])),xpd=T,adj=0,cex=0.8,col=rgb(0.9, 0,0))
          #text(min(RR),115,paste0('Warning: GO and NOGO zones are overlaped at some sample sizes'),xpd=T,adj=0,cex=0.5,col=rgb(0.9,0.6,0))
          #text(min(RR),110,paste0('The zones are classfied by cutoff of ', overlap.option),xpd=T,adj=0,cex=0.5,col=rgb(0.9,0.6,0))
        }
        if(direction=="Less"){
          text(min(RR),130,bquote(GO~symbol("\336")~Number~of~responders~symbol("\243")~.(temp$cutoff[2])),xpd=T,adj=0,cex=0.8,col=rgb(0, 0.7,0))
          text(min(RR),120,bquote(NOGO~symbol("\336")~Number~of~responders~symbol("\76")~.(temp$cutoff[2])),xpd=T,adj=0,cex=0.8,col=rgb(0.9, 0,0))
          #text(min(RR),115,paste0('Warning: GO and NOGO zones are overlaped at some sample sizes'),xpd=T,adj=0,cex=0.5,col=rgb(0.9,0.6,0))
          #text(min(RR),110,paste0('The zones are classfied by cutoff of ', overlap.option),xpd=T,adj=0,cex=0.5,col=rgb(0.9,0.6,0))
        }
      }
      if(unsatisfied.flag==1){
        text(min(RR),130,paste0('Warning: We could not find classification of zones to satisfy your both GO and NOGO criterions at some sample sizes'),xpd=T,adj=0,cex=0.5,col=rgb(0.9,0,0))
        text(min(RR),120,paste0('Please check and modify your desicion rule or the range of sample size!'),xpd=T,adj=0,cex=0.5,col=rgb(0.9,0,0))

      }

    }
  }

##########################################################################################


  return(list(go_prob=go_prob,nogo_prob=nogo_prob,
              inconclusive_prob=inconclusive_prob,
              overlap.flag=overlap.flag,overlap.option=overlap.option,
              unsatisfied.flag=unsatisfied.flag,cutoff=temp$cutoff,true_cutoff=c(true_go_cutoff,true_nogo_cutoff))
       )


}

#######################################%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Fix_SS_SAB_Bin_Prob(n=46,a=1/3,b=1,
#                     CT1.go=0.18,
#                     false.go.CT1=TRUE,
#                     FGR.CT1=0.8,
#                     CT1.nogo=0.18,
#                     false.nogo.CT1=TRUE,
#                     FNGR.CT1=0.8,
#                     CT2.go=0.3,
#                     false.go.CT2=FALSE,
#                     FGR.CT2=0.5,
#                     CT2.nogo=0.3,
#                     false.nogo.CT2=FALSE,
#                     FNGR.CT2=0.5,
#                     RR=c(0.1,0.65),overlap.option='GO',
#                     plot.figure=TRUE,
#                     method='Frequentist',direction='Greater',logic.go='and',logic.nogo='and')
#
#
# Vary_SS_SAB_Bin_Prob(nmin=10,nmax=300,a=1/3,b=1,
#                                CT1.go=0.18,
#                                false.go.CT1=TRUE,FGR.CT1=0.8,
#                                CT1.nogo=0.18,
#                                false.nogo.CT1=TRUE,FNGR.CT1=0.8,
#                                CT2.go=0.3,
#                                false.go.CT2=TRUE, FGR.CT2=0.5,
#                                CT2.nogo=0.3,
#                                false.nogo.CT2=TRUE, FNGR.CT2=0.5,
#                                RR=c(0.25),overlap.option='GO',
#                                plot.cutoff=TRUE,plot.prob=TRUE,
#                                method="Bayesian",direction="Greater",logic.go='and',logic.nogo='and')
#

######################################################################
Vary_SS_SAB_Bin_Prob<-function(nmin=10,nmax=300,a=1/3,b=1,
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
                               method="Bayesian",direction="Greater",logic.go='and',logic.nogo='and'){
  nseq=unique(round(c(seq(nmin,(nmax+nmin)/2,length=6)[-6],seq((nmax+nmin)/2,nmax,length=6))))
  go_prob<-matrix(NA,ncol=length(nseq),nrow=length(RR))
  nogo_prob<-matrix(NA,ncol=length(nseq),nrow=length(RR))
  inconclusive_prob<-matrix(NA,ncol=length(nseq),nrow=length(RR))
  go_cutoff<-rep(NA,length(nseq))
  nogo_cutoff<-rep(NA,length(nseq))
  index=1
  n_unsatisfied=NA
  n_overlap=NA
  for(i in nseq){
    temp=SAB_Bin_Cutoff(n=i,a=a,b=b,
                        CT1.go=CT1.go,
                        false.go.CT1 = false.go.CT1, FGR.CT1=FGR.CT1,
                        CT1.nogo =CT1.nogo,
                        false.nogo.CT1 = false.nogo.CT1, FNGR.CT1 = FNGR.CT1,
                        CT2.go=CT2.go,
                        false.go.CT2 = false.go.CT2, FGR.CT2=FGR.CT2,
                        CT2.nogo = CT2.nogo,
                        false.nogo.CT2 = false.nogo.CT2, FNGR.CT2 = FNGR.CT2,
                        method=method,direction=direction,logic.go=logic.go,logic.nogo=logic.nogo)
    if(direction=='Greater'){
      if(all(temp$flag==0)&temp$overlap==0){
        go_prob[,index]=1-pbinom(temp$cutoff[1]-1,i,RR)
        nogo_prob[,index]=pbinom(temp$cutoff[2]-1,i,RR)
        inconclusive_prob[,index]=1-go_prob[,index]-nogo_prob[,index]
        go_cutoff[index]=temp$cutoff[1]
        nogo_cutoff[index]=temp$cutoff[2]
      }
      if(all(temp$flag==0)&temp$overlap==1&overlap.option=="GO"){
        go_prob[,index]=1-pbinom(temp$cutoff[1]-1,i,RR)
        nogo_prob[,index]=pbinom(temp$cutoff[1]-1,i,RR)
        inconclusive_prob[,index]=0
        go_cutoff[index]=temp$cutoff[1]
        nogo_cutoff[index]=temp$cutoff[2]
        n_overlap=c(n_overlap,i)
      }
      if(all(temp$flag==0)&temp$overlap==1&overlap.option=="NOGO"){
        go_prob[,index]=1-pbinom(temp$cutoff[2]-1,i,RR)
        nogo_prob[,index]=pbinom(temp$cutoff[2]-1,i,RR)
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
        go_prob[,index]=pbinom(temp$cutoff[1],i,RR)
        nogo_prob[,index]=1-pbinom(temp$cutoff[2],i,RR)
        inconclusive_prob[,index]=1-go_prob[,index]-nogo_prob[,index]
        go_cutoff[index]=temp$cutoff[1]
        nogo_cutoff[index]=temp$cutoff[2]
      }
      if(all(temp$flag==0)&temp$overlap==1&overlap.option=="GO"){
        go_prob[,index]=pbinom(temp$cutoff[1],i,RR)
        nogo_prob[,index]=1-pbinom(temp$cutoff[1],i,RR)
        inconclusive_prob[,index]=1-go_prob[,index]-nogo_prob[,index]
        go_cutoff[index]=temp$cutoff[1]
        nogo_cutoff[index]=temp$cutoff[2]
        n_overlap=c(n_overlap,i)
      }
      if(all(temp$flag==0)&temp$overlap==1&overlap.option=="NOGO"){
        go_prob[,index]=pbinom(temp$cutoff[2],i,RR)
        nogo_prob[,index]=1-pbinom(temp$cutoff[2],i,RR)
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
      text((nmin+nmax)/2,120,paste0("True RR=",round(RR[j],3)),xpd=T,adj=0.5,cex=0.8)
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
    p.obs_go=go_cutoff/nseq
    p.obs_nogo=nogo_cutoff/nseq
    P_go=p.obs_go*100
    P_nogo=p.obs_nogo*100
    par(mfrow=c(1,1),pty="m",bty="l",cex=1.4)
    ylim_max=min(max(c(P_go,P_nogo,round(CT1.go*100,3),round(CT2.go*100,3),round(CT1.nogo*100,3),round(CT2.nogo*100,3)),na.rm=TRUE)+10,100,na.rm=TRUE)
    ylim_min=max(min(c(P_go,P_nogo,round(CT1.go*100,3),round(CT2.go*100,3),round(CT1.nogo*100,3),round(CT2.nogo*100,3)),na.rm=TRUE)-10,0,na.rM=TRUE)
    plot(NA,NA,xlab='Sample size',ylab="Observed RR(%)",xlim=c(max(0,range(nseq)[1]-10),range(nseq)[2]+10),ylim=c(ylim_min,ylim_max),type="n",axes=F,col=rgb(1,0,0),lty=1,lwd=2)
    axis(1, at=nseq,labels=nseq)
    axis(2, at=seq(0,100,10),labels=T)
    #text(rep(-0.6,6),seq(0,100,20),seq(0,300,20),adj=1,xpd=T)
    box()
    lines(nseq,P_nogo,col=rgb(0.9,0,0),lwd=2)
    lines(nseq,P_go,col=rgb(0,0.7,0),lwd=2)
    legend('bottomright',legend=c("Cut off of GO","Cut off of NOGO"),
           col=c(rgb(0,0.7,0),rgb(0.9,0,0)),
           lwd=c(2,2),
           lty=c(1,1),cex=0.5)
    if(length(n_overlap)!=0){
      text(nmin,ylim_max+10,paste0('Warning: GO and NOGO zones are overlaped at some sample sizes'),xpd=T,adj=0,cex=0.5,col=rgb(0.9,0.6,0))
      #text(nmin,110,paste(n_overlap),xpd=T,adj=0,cex=0.5,col=rgb(0.9,0.6,0))
      #print(n_overlap)
      text(nmin,ylim_max+5,paste0('The zones are classfied by cutoff of ', overlap.option),xpd=T,adj=0,cex=0.5,col=rgb(0.9,0.6,0))
    }
    if(length(n_unsatisfied)!=0){
      text(nmin,ylim_max+10,paste0('Warning: We could not find classification of zones to satisfy your both GO and NOGO criterions at some sample sizes'),xpd=T,adj=0,cex=0.5,col=rgb(0.9,0,0))
      text(nmin,ylim_max+5,paste0('Please check and modify your desicion rule or the range of sample size!'),xpd=T,adj=0,cex=0.5,col=rgb(0.9,0,0))

    }


  }
  return(list(n_unsatisfied,n_overlap))


}
Interim_SAB<-function(num_interim=3,interim_n=c(50,100,150),
                      CT1.go=c(0.25,0.25,0.25),
                      false.go.CT1=c(TRUE,TRUE,TRUE),FGR.CT1=c(0.25,0.25,0.25),
                      CT1.nogo=c(0.25,0.25,0.25),
                      false.nogo.CT1=c(TRUE,TRUE,TRUE),FNGR.CT1=c(0.25,0.25,0.25),
                      CT2.go=c(0.25,0.25,0.25),
                      false.go.CT2=c(TRUE,TRUE,TRUE),FGR.CT2=c(0.25,0.25,0.25),
                      CT2.nogo=c(0.25,0.25,0.25),
                      false.nogo.CT2=c(TRUE,TRUE,TRUE),FNGR.CT2=c(0.25,0.25,0.25),
                      overlap.option=c('GO','GO','GO'),logic.go=c('and','and','and'),
                      logic.nogo=c('and','and','and'),
                      method='Bayesian',direction="Greater",nsim_IA=10000,seed.num=369,
                      task=c('Futility','Superiority','Futility and superiority'),
                      a=1/3,b=1,RR=c(0.25)){
  interim_n=sort(interim_n)
  go_cutoff<-rep(NA,num_interim)
  nogo_cutoff<-rep(NA,num_interim)
  true_go_cutoff<-rep(NA,num_interim)
  true_nogo_cutoff<-rep(NA,num_interim)
  overlap<-rep(NA,num_interim)
  for(i in 1:num_interim){
    temp<-SAB_Bin_Cutoff(n=interim_n[i],a=a,b=b,
                         CT1.go=CT1.go[i],
                         false.go.CT1=false.go.CT1[i],FGR.CT1=FGR.CT1[i],
                         CT1.nogo=CT1.nogo[i],
                         false.nogo.CT1=false.nogo.CT1[i],FNGR.CT1=FNGR.CT1[i],
                         CT2.go=CT2.go[i],
                         false.go.CT2=false.go.CT2[i], FGR.CT2=FGR.CT2[i],
                         CT2.nogo=CT2.nogo[i],
                         false.nogo.CT2=false.nogo.CT2[i], FNGR.CT2=FNGR.CT2[i],
                         method=method,direction=direction[i],logic.go=logic.go[i],
                         logic.nogo=logic.nogo[i]
    )
    go_cutoff[i]<-temp$cutoff[1]
    nogo_cutoff[i]<-temp$cutoff[2]
    overlap[i]<-temp$overlap
    ###SAB####
    if(direction[i]=='Greater'){
      true_go_cutoff[i]<-ifelse(overlap[i]==0,go_cutoff[i],(overlap.option[i]=='GO')*go_cutoff[i]+(overlap.option[i]=='NOGO')*(nogo_cutoff[i]))
      true_nogo_cutoff[i]<-ifelse(overlap[i]==0,nogo_cutoff[i],(overlap.option[i]=='GO')*(go_cutoff[i])+(overlap.option[i]=='NOGO')*(nogo_cutoff[i]))
    }
    if(direction[i]=='Less'){
      true_go_cutoff[i]<-ifelse(overlap[i]==0,go_cutoff[i],(overlap.option[i]=='GO')*go_cutoff[i]+(overlap.option[i]=='NOGO')*(nogo_cutoff[i]))
      true_nogo_cutoff[i]<-ifelse(overlap[i]==0,nogo_cutoff[i],(overlap.option[i]=='GO')*(go_cutoff[i])+(overlap.option[i]=='NOGO')*(nogo_cutoff[i]))
    }
    ####
  }
  temptable=c()
  for(meanindex in 1:length(RR)){

    set.seed(seed.num)
    diff_interim_n<-diff(interim_n)
    generate_n<-c(interim_n[1],diff_interim_n)
    ###SAB####
    sim_data<-sapply(generate_n,rbinom,n=nsim_IA,prob=RR[meanindex])
    ######
    cum_sim_data<-t(apply(sim_data,1,cumsum))
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
        # "<=" is replaced with "<"
        inconclusive_matrix[,j]<-rep(1,nsim_IA)-go_matrix[,j]-nogo_matrix[,j]
      }
      if(direction[j]=='Less'){
        go_matrix[,j]<-cum_sim_data[,j]<=true_go_cutoff[j]
        nogo_matrix[,j]<-cum_sim_data[,j]>true_nogo_cutoff[j] 
        # ">=" is replaced with ">"
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
      table[2,j]=RR[meanindex]
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
    table[2,num_interim+1]=RR[meanindex]
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
    tablerowname<-c('Sample size','True RR','Task','Success','To next interim/final or inconclusive',
                    'Stop','Superority/Futility zone','Warning',
                    'Cut off for NOGO rule',
                    'Cut off for GO rule')
    table<-cbind(tablerowname,rep(meanindex,10),table)
    colnames(table)<-c(" ",'Setting',tablecolname)
    temptable=rbind(temptable,table)

  }
  return(temptable)
}

############@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
# Fix_SS_SAB_Bin_Prob(n=25,a=1/3,b=1,
#                               CT1.go=0.23,
#                               false.go.CT1=TRUE,FGR.CT1=0.3,
#                               CT1.nogo=0.23,
#                               false.nogo.CT1=TRUE,FNGR.CT1=0.2,
#                               CT2.go=0.3,
#                               false.go.CT2=FALSE, FGR.CT2=0.5,
#                               CT2.nogo=0.3,
#                               false.nogo.CT2=FALSE, FNGR.CT2=0.5,
#                               RR=c(0.15,0.6),overlap.option='GO',
#                               plot.figure=TRUE,
#                               method='Frequentist',direction='Greater',logic.go='and',logic.nogo='and')
#

