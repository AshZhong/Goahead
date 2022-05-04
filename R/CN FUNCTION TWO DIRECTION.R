library(shiny)


CN_OC_table<-function(n=c(100,200),
                      prior.mean.1=1/3,prior.sd.1=1,
                      prior.k.1=1,prior.a.1=1,prior.b.1=1,
                      sd.1=1,ssq.1=1,
                      CT1.go.1=0.25,
                      false.go.CT1.1=TRUE,FGR.CT1.1=0.25,
                      CT1.nogo.1=0.25,
                      false.nogo.CT1.1=TRUE,FNGR.CT1.1=0.25,
                      CT2.go.1=0.3,
                      false.go.CT2.1=TRUE, FGR.CT2.1=0.5,
                      CT2.nogo.1=0.3,
                      false.nogo.CT2.1=TRUE, FNGR.CT2.1=0.5,
                      direction.1='Greater',
                      logic.go.1='or',logic.nogo.1='and',
                      overlap.option.1='GO',
                      prior.mean.2=1/3,prior.sd.2=1,
                      prior.k.2=1,prior.a.2=1,prior.b.2=1,
                      sd.2=1,ssq.2=1,
                      CT1.go.2=0.2,
                      false.go.CT1.2=TRUE,FGR.CT1.2=0.25,
                      CT1.nogo.2=0.25,
                      false.nogo.CT1.2=TRUE,FNGR.CT1.2=0.25,
                      CT2.go.2=0.3,
                      false.go.CT2.2=TRUE, FGR.CT2.2=0.5,
                      CT2.nogo.2=0.3,
                      false.nogo.CT2.2=TRUE, FNGR.CT2.2=0.5,
                      direction.2='Greater',
                      logic.go.2='and',logic.nogo.2='and',
                      overlap.option.2='GO',
                      logic.go.points='and',
                      logic.nogo.points='or',
                      fix.var=TRUE,noninfo=TRUE,
                      seed.num=425,nsim=10000,
                      stop.criterion=10^-3,
                      mean=rbind(c(1,2,0.5,0.3),c(1,3,0.3,0.4)),
                      method='Bayesian'){
  temp_table<-c()
  for(j in 1:length(n)){
  temp1=SAN_Normal_Cutoff(n=n[j],
                         prior.mean=prior.mean.1,prior.sd = prior.sd.1,
                         prior.k = prior.k.1,prior.a=prior.a.1,prior.b=prior.b.1,
                         sd=sd.1,ssq=ssq.1,
                         CT1.go=CT1.go.1,
                         false.go.CT1 = false.go.CT1.1, FGR.CT1=FGR.CT1.1,
                         CT1.nogo=CT1.nogo.1,
                         false.nogo.CT1 = false.nogo.CT1.1, FNGR.CT1 = FNGR.CT1.1,
                         CT2.go=CT2.go.1,
                         false.go.CT2 = false.go.CT2.1, FGR.CT2=FGR.CT2.1,
                         CT2.nogo=CT2.nogo.1,
                         false.nogo.CT2 = false.nogo.CT2.1, FNGR.CT2 = FNGR.CT2.1,
                         logic.go=logic.go.1,logic.nogo = logic.nogo.1,
                         direction=direction.1,
                         method=method,
                         fix.var=fix.var,noninfo=noninfo,
                         seed.num=seed.num,
                         stop.criterion=stop.criterion)
  temp2=SAN_Normal_Cutoff(n=n[j],
                          prior.mean=prior.mean.2,prior.sd = prior.sd.2,
                          prior.k = prior.k.2,prior.a=prior.a.2,prior.b=prior.b.2,
                          sd=sd.2,ssq=ssq.2,
                          CT1.go=CT1.go.2,
                          false.go.CT1 = false.go.CT1.2, FGR.CT1=FGR.CT1.2,
                          CT1.nogo=CT1.nogo.2,
                          false.nogo.CT1 = false.nogo.CT1.2, FNGR.CT1 = FNGR.CT1.2,
                          CT2.go=CT2.go.2,
                          false.go.CT2 = false.go.CT2.2, FGR.CT2=FGR.CT2.2,
                          CT2.nogo=CT2.nogo.2,
                          false.nogo.CT2 = false.nogo.CT2.2, FNGR.CT2 = FNGR.CT2.2,
                          logic.go=logic.go.2,logic.nogo = logic.nogo.2,
                          direction=direction.2,
                          method=method,
                          fix.var=fix.var,noninfo=noninfo,
                          seed.num=seed.num,
                          stop.criterion=stop.criterion)


    true_go_cutoff1<-ifelse(temp1$overlap==0,temp1$cutoff[1],(overlap.option.1=='GO')*temp1$cutoff[1]+(overlap.option.1=='NOGO')*(temp1$cutoff[2]))
    true_nogo_cutoff1<-ifelse(temp1$overlap==0,temp1$cutoff[2],(overlap.option.1=='GO')*(temp1$cutoff[1])+(overlap.option.1=='NOGO')*(temp1$cutoff[2]))
    true_go_cutoff2<-ifelse(temp2$overlap==0,temp2$cutoff[1],(overlap.option.2=='GO')*temp2$cutoff[1]+(overlap.option.2=='NOGO')*(temp2$cutoff[2]))
    true_nogo_cutoff2<-ifelse(temp2$overlap==0,temp2$cutoff[2],(overlap.option.2=='GO')*(temp2$cutoff[1])+(overlap.option.2=='NOGO')*(temp2$cutoff[2]))

    mean<-as.matrix(mean)
    num_setting<-dim(mean)[1]

    go_prob<-rep(NA,num_setting)
    nogo_prob<-rep(NA,num_setting)
    inconclusive_prob<-rep(NA,num_setting)
    table_RR<-rep(NA,num_setting)
    table_true_go_cutoff<-rep(NA,num_setting)
    table_true_nogo_cutoff<-rep(NA,num_setting)
    table_go_cutoff<-rep(NA,num_setting)
    table_nogo_cutoff<-rep(NA,num_setting)
    table_warning<-rep(NA,num_setting)
    for(index.setting in 1: num_setting){
      rho=mean[index.setting,3]
      sigma<-matrix(NA,ncol=2,nrow=2)
      sigma[1,1]<-sd.1^2/n[j]
      sigma[2,2]<-sd.2^2/n[j]
      sigma[1,2]<-rho*sd.1*sd.2/n[j]
      sigma[2,1]<-rho*sd.1*sd.2/n[j]
      if(logic.go.points=='and'){
        if(direction.1=='Greater'&direction.2=='Greater'){
          go_prob[index.setting]<-pmvnorm(lower=c(true_go_cutoff1,true_go_cutoff2),upper=c(Inf,Inf),mean=mean[index.setting,1:2],sigma=sigma)
          table_true_go_cutoff[index.setting]<-HTML('X>=',round(true_go_cutoff1,3),logic.go.points,' Y>=',round(true_go_cutoff2,3))
          table_go_cutoff[index.setting]<-HTML('X>=',round(temp1$cutoff[1],3),logic.go.points,' Y>=',round(temp2$cutoff[1],3))
        }
        if(direction.1=='Greater'&direction.2=='Less'){
          go_prob[index.setting]<-pmvnorm(lower=c(true_go_cutoff1,-Inf),upper=c(Inf,true_go_cutoff2),mean=mean[index.setting,1:2],sigma=sigma)
          table_true_go_cutoff[index.setting]<-HTML('X>=',round(true_go_cutoff1,3),logic.go.points,' Y<=',round(true_go_cutoff2,3))
          table_go_cutoff[index.setting]<-HTML('X>=',round(temp1$cutoff[1],3),logic.go.points,' Y<=',round(temp2$cutoff[1],3))
        }
        if(direction.1=='Less'&direction.2=='Greater'){
          go_prob[index.setting]<-pmvnorm(lower=c(-Inf,true_go_cutoff2),upper=c(true_go_cutoff1,Inf),mean=mean[index.setting,1:2],sigma=sigma)
          table_true_go_cutoff[index.setting]<-HTML('X<=',round(true_go_cutoff1,3),logic.go.points,' Y>=',round(true_go_cutoff2,3))
          table_go_cutoff[index.setting]<-HTML('X<=',round(temp1$cutoff[1],3),logic.go.points,' Y>=',round(temp2$cutoff[1],3))
        }
        if(direction.1=='Less'&direction.2=='Less'){
          go_prob[index.setting]<-pmvnorm(lower=c(-Inf,-Inf),upper=c(true_go_cutoff1,true_go_cutoff2),mean=mean[index.setting,1:2],sigma=sigma)
          table_true_go_cutoff[index.setting]<-HTML('X<=',round(true_go_cutoff1,3),logic.go.points,' Y<=',round(true_go_cutoff2,3))
          table_go_cutoff[index.setting]<-HTML('X<=',round(temp1$cutoff[1],3),logic.go.points,' Y<=',round(temp2$cutoff[1],3))

        }
      }
      if(logic.go.points=='or'){
        if(direction.1=='Greater'&direction.2=='Greater'){
          go_prob[index.setting]<-1-pmvnorm(lower=c(-Inf,-Inf),upper=c(true_go_cutoff1,true_go_cutoff2),mean=mean[index.setting,1:2],sigma=sigma)
          table_true_go_cutoff[index.setting]<-HTML('X>=',round(true_go_cutoff1,3),logic.go.points,' Y>=',round(true_go_cutoff2,3))
          table_go_cutoff[index.setting]<-HTML('X>=',round(temp1$cutoff[1],3),logic.go.points,' Y>=',round(temp2$cutoff[1],3))

        }
        if(direction.1=='Greater'&direction.2=='Less'){
          go_prob[index.setting]<-1-pmvnorm(lower=c(-Inf,true_go_cutoff2),upper=c(true_go_cutoff1,Inf),mean=mean[index.setting,1:2],sigma=sigma)
          table_true_go_cutoff[index.setting]<-HTML('X>=',round(true_go_cutoff1,3),logic.go.points,' Y<=',round(true_go_cutoff2,3))
          table_go_cutoff[index.setting]<-HTML('X>=',round(temp1$cutoff[1],3),logic.go.points,' Y<=',round(temp2$cutoff[1],3))
        }
        if(direction.1=='Less'&direction.2=='Greater'){
          go_prob[index.setting]<-1-pmvnorm(lower=c(true_go_cutoff1,-Inf),upper=c(Inf,true_go_cutoff2),mean=mean[index.setting,1:2],sigma=sigma)
          table_true_go_cutoff[index.setting]<-HTML('X<=',round(true_go_cutoff1,3),logic.go.points,' Y>=',round(true_go_cutoff2,3))
          table_go_cutoff[index.setting]<-HTML('X<=',round(temp1$cutoff[1],3),logic.go.points,' Y>=',round(temp2$cutoff[1],3))
        }
        if(direction.1=='Less'&direction.2=='Less'){
          go_prob[index.setting]<-1-pmvnorm(lower=c(true_go_cutoff1,true_go_cutoff2),upper=c(Inf,Inf),mean=mean[index.setting,1:2],sigma=sigma)
          table_true_go_cutoff[index.setting]<-HTML('X<=',round(true_go_cutoff1,3),logic.go.points,' Y<=',round(true_go_cutoff2,3))
          table_go_cutoff[index.setting]<-HTML('X<=',round(temp1$cutoff[1],3),logic.go.points,' Y<=',round(temp2$cutoff[1],3))
        }
      }

      if(logic.nogo.points=='and'){
        if(direction.1=='Greater'&direction.2=='Greater'){
          nogo_prob[index.setting]<-pmvnorm(lower=c(-Inf,-Inf),upper=c(true_nogo_cutoff1,true_nogo_cutoff2),mean=mean[index.setting,1:2],sigma=sigma)
          table_true_nogo_cutoff[index.setting]<-HTML('X<',round(true_nogo_cutoff1,3),logic.nogo.points,' Y<',round(true_nogo_cutoff2,3))
          table_nogo_cutoff[index.setting]<-HTML('X<',round(temp1$cutoff[2],3),logic.nogo.points,' Y<',round(temp2$cutoff[2],3))

        }
        if(direction.1=='Greater'&direction.2=='Less'){
          nogo_prob[index.setting]<-pmvnorm(lower=c(-Inf,true_nogo_cutoff2),upper=c(true_nogo_cutoff1,Inf),mean=mean[index.setting,1:2],sigma=sigma)
          table_true_nogo_cutoff[index.setting]<-HTML('X<',round(true_nogo_cutoff1,3),logic.nogo.points,' Y>',round(true_nogo_cutoff2,3))
          table_nogo_cutoff[index.setting]<-HTML('X<',round(temp1$cutoff[2],3),logic.nogo.points,' Y>',round(temp2$cutoff[2],3))

        }
        if(direction.1=='Less'&direction.2=='Greater'){
          nogo_prob[index.setting]<-pmvnorm(lower=c(true_nogo_cutoff1,-Inf),upper=c(true_nogo_cutoff1,true_nogo_cutoff2),mean=mean[index.setting,1:2],sigma=sigma)
          table_true_nogo_cutoff[index.setting]<-HTML('X>',round(true_nogo_cutoff1,3),logic.nogo.points,' Y<',round(true_nogo_cutoff2,3))
          table_nogo_cutoff[index.setting]<-HTML('X>',round(temp1$cutoff[2],3),logic.nogo.points,' Y<',round(temp2$cutoff[2],3))

        }
        if(direction.1=='Less'&direction.2=='Less'){
          nogo_prob[index.setting]<-pmvnorm(lower=c(true_nogo_cutoff1,true_nogo_cutoff2),upper=c(Inf,Inf),mean=mean[index.setting,1:2],sigma=sigma)
          table_true_nogo_cutoff[index.setting]<-HTML('X>',round(true_nogo_cutoff1,3),logic.nogo.points,' Y>',round(true_nogo_cutoff2,3))
          table_nogo_cutoff[index.setting]<-HTML('X>',round(temp1$cutoff[2],3),logic.nogo.points,' Y>',round(temp2$cutoff[2],3))

        }
      }

      if(logic.nogo.points=='or'){
        if(direction.1=='Greater'&direction.2=='Greater'){
          nogo_prob[index.setting]<-1-pmvnorm(lower=c(true_nogo_cutoff1,true_nogo_cutoff2),upper=c(Inf,Inf),mean=mean[index.setting,1:2],sigma=sigma)
          table_true_nogo_cutoff[index.setting]<-HTML('X<',round(true_nogo_cutoff1,3),logic.nogo.points,' Y<',round(true_nogo_cutoff2,3))
          table_nogo_cutoff[index.setting]<-HTML('X<',round(temp1$cutoff[2],3),logic.nogo.points,' Y<',round(temp2$cutoff[2],3))

        }
        if(direction.1=='Greater'&direction.2=='Less'){
          nogo_prob[index.setting]<-1-pmvnorm(lower=c(true_nogo_cutoff1,-Inf),upper=c(true_nogo_cutoff1,true_nogo_cutoff2),mean=mean[index.setting,1:2],sigma=sigma)
          table_true_nogo_cutoff[index.setting]<-HTML('X<',round(true_nogo_cutoff1,3),logic.nogo.points,' Y>',round(true_nogo_cutoff2,3))
          table_nogo_cutoff[index.setting]<-HTML('X<',round(temp1$cutoff[2],3),logic.nogo.points,' Y>',round(temp2$cutoff[2],3))


        }
        if(direction.1=='Less'&direction.2=='Greater'){
          nogo_prob[index.setting]<-1-pmvnorm(lower=c(-Inf,true_nogo_cutoff2),upper=c(true_nogo_cutoff1,Inf),mean=mean[index.setting,1:2],sigma=sigma)
          table_true_nogo_cutoff[index.setting]<-HTML('X>',round(true_nogo_cutoff1,3),logic.nogo.points,' Y<',round(true_nogo_cutoff2,3))
          table_nogo_cutoff[index.setting]<-HTML('X>',round(temp1$cutoff[2],3),logic.nogo.points,' Y<',round(temp2$cutoff[2],3))


        }
        if(direction.1=='Less'&direction.2=='Less'){
          nogo_prob[index.setting]<-1-pmvnorm(lower=c(-Inf,-Inf),upper=c(true_nogo_cutoff1,true_nogo_cutoff2),mean=mean[index.setting,1:2],sigma=sigma)
          table_true_nogo_cutoff[index.setting]<-HTML('X>',round(true_nogo_cutoff1,3),logic.nogo.points,' Y>',round(true_nogo_cutoff2,3))
          table_nogo_cutoff[index.setting]<-HTML('X>',round(temp1$cutoff[2],3),logic.nogo.points,' Y>',round(temp2$cutoff[2],3))

        }
      }

      inconclusive_prob[index.setting]<-1-go_prob[index.setting]-nogo_prob[index.setting]
      table_RR[index.setting]<-paste0(mean[index.setting,],collapse = ',')
      table_warning[index.setting]<-'None'
      if(temp1$overlap==1&temp2$overlap==1){
        table_warning[index.setting]<-HTML('GO/NOGO zones for Endpoint 1 are overlapped, classifed by ',overlap.option.1,' zone.<br>',
                                           'GO/NOGO zones for Endpoint 2 are overlapped, classifed by ',overlap.option.2,' zone.')
      }
      if(temp1$overlap==1&temp2$overlap==0){
        table_warning[index.setting]<-HTML('GO/NOGO zones for Endpoint 1 are overlapped, classifed by ',overlap.option.1,' zone')
      }
      if(temp1$overlap==0&temp2$overlap==1){
        table_warning[index.setting]<-HTML('GO/NOGO zones for Endpoint 2 are overlapped, classifed by ',overlap.option.2,' zone')
      }

    }
    table_n<-rep(n[j],num_setting)
    table_go_prob<-round(go_prob,3)
    table_nogo_prob<-round(nogo_prob,3)
    table_inconclusive_prob<-round(inconclusive_prob,3)
    table_zone<-paste0(table_true_go_cutoff,'/',table_true_nogo_cutoff)
    temp<-cbind(table_n,table_RR,table_go_prob,table_nogo_prob,table_inconclusive_prob,table_zone,table_warning,table_go_cutoff,table_nogo_cutoff)
    temp_table<-rbind(temp_table,temp)
  }
  colnames(temp_table)<-c('Sample size',HTML('&mu;<sub>1</sub>,&mu;<sub>2</sub>,Corr'),'Pr(GO)','Pr(NOGO)','Pr(Inconclusive)','GO/NOGO zones','Warning','Cut off for GO','Cut off for NOGO')
  return(temp_table)
}


CN_OC_table_interim<-function(num_interim=3,n_interim=c(100,200,300),
                              CT1.go.1=rep(0.5,3),
                              false.go.CT1.1=c(TRUE,TRUE,FALSE),FGR.CT1.1=rep(0.95,3),
                              CT1.nogo.1=rep(0.5,3),
                              false.nogo.CT1.1=c(TRUE,TRUE,FALSE),FNGR.CT1.1=rep(0.25,3),
                              CT2.go.1=rep(0.5,3),
                              false.go.CT2.1=rep(TRUE,3), FGR.CT2.1=rep(0.95,3),
                              CT2.nogo.1=rep(0.4,3),
                              false.nogo.CT2.1=rep(TRUE,3), FNGR.CT2.1=rep(0.5,3),
                              direction.1=rep('Greater',3),
                              logic.go.1=rep('and',3),logic.nogo.1=rep('or',3),
                              overlap.option.1=rep('GO',3),
                              CT1.go.2=rep(0.25,3),
                              false.go.CT1.2=c(TRUE,TRUE,FALSE),FGR.CT1.2=rep(0.95,3),
                              CT1.nogo.2=rep(0.25,3),
                              false.nogo.CT1.2=c(TRUE,TRUE,FALSE),FNGR.CT1.2=rep(0.25,3),
                              CT2.go.2=rep(0.3,3),
                              false.go.CT2.2=rep(TRUE,3), FGR.CT2.2=rep(0.95,3),
                              CT2.nogo.2=rep(0.3,3),
                              false.nogo.CT2.2=rep(TRUE,3), FNGR.CT2.2=rep(0.5,3),
                              direction.2=rep('Greater',3),
                              logic.go.2=rep('and',3),logic.nogo.2=rep('or',3),
                              overlap.option.2=rep('GO',3),
                              logic.go.points=c('and','or','and'),
                              logic.nogo.points=c('or','or','and'),
                              task=c('Futility','Superiority','Futility and superiority'),
                              prior.mean.1=1/3,prior.sd.1=1,
                              prior.k.1=1,prior.a.1=1,prior.b.1=1,
                              sd.1=1,ssq.1=1,
                              a.2=1/3,b.2=1,
                              prior.mean.2=1/3,prior.sd.2=1,
                              prior.k.2=1,prior.a.2=1,prior.b.2=1,
                              sd.2=1,ssq.2=1,
                              fix.var=TRUE,noninfo=TRUE,
                              nsim=10000,
                              stop.criterion=10^-3,
                              mean=rbind(c(0.5,0.4,0.3),c(1,3,-0.1)),
                              method='Bayesian',nsim_IA=10000,seed.num=425){


  go_matrix<-array(NA,c(nsim_IA,num_interim,dim(mean)[1]))
  nogo_matrix<-array(NA,c(nsim_IA,num_interim,dim(mean)[1]))
  inconclusive_matrix<-array(NA,c(nsim_IA,num_interim,dim(mean)[1]))
  IA_go_matrix<-array(NA,c(nsim_IA,num_interim,dim(mean)[1]))
  n_interim<-sort(n_interim)
  num_setting<-dim(mean)[1]
  table_RR<-array(NA,c(num_interim,num_setting))
  table_true_go_cutoff<-array(NA,c(num_interim,num_setting))
  table_true_nogo_cutoff<-array(NA,c(num_interim,num_setting))
  table_go_cutoff<-array(NA,c(num_interim,num_setting))
  table_nogo_cutoff<-array(NA,c(num_interim,num_setting))
  table_warning<-array(NA,c(num_interim,num_setting))
  table_zone<-array(NA,c(num_interim,num_setting))
  table_n<-array(NA,c(num_interim,num_setting))
  for(interim.index in 1:num_interim){
    temp1=SAN_Normal_Cutoff(n=n_interim[interim.index],
                            prior.mean=prior.mean.1,prior.sd = prior.sd.1,
                            prior.k = prior.k.1,prior.a=prior.a.1,prior.b=prior.b.1,
                            sd=sd.1,ssq=ssq.1,
                            CT1.go=CT1.go.1[interim.index],
                            false.go.CT1 = false.go.CT1.1[interim.index], FGR.CT1=FGR.CT1.1[interim.index],
                            CT1.nogo=CT1.nogo.1[interim.index],
                            false.nogo.CT1 = false.nogo.CT1.1[interim.index], FNGR.CT1 = FNGR.CT1.1[interim.index],
                            CT2.go=CT2.go.1[interim.index],
                            false.go.CT2 = false.go.CT2.1[interim.index], FGR.CT2=FGR.CT2.1[interim.index],
                            CT2.nogo=CT2.nogo.1[interim.index],
                            false.nogo.CT2 = false.nogo.CT2.1[interim.index], FNGR.CT2 = FNGR.CT2.1[interim.index],
                            direction=direction.1[interim.index],
                            logic.go = logic.go.1[interim.index],logic.nogo=logic.nogo.1[interim.index],
                            method=method,
                            fix.var=fix.var,noninfo=noninfo,
                            seed.num=seed.num,
                            stop.criterion=stop.criterion)

      true_go_cutoff1<-ifelse(temp1$overlap==0,temp1$cutoff[1],(overlap.option.1[interim.index]=='GO')*temp1$cutoff[1]+(overlap.option.1[interim.index]=='NOGO')*(temp1$cutoff[2]))
      true_nogo_cutoff1<-ifelse(temp1$overlap==0,temp1$cutoff[2],(overlap.option.1[interim.index]=='GO')*(temp1$cutoff[1])+(overlap.option.1[interim.index]=='NOGO')*(temp1$cutoff[2]))



    temp2=SAN_Normal_Cutoff(n=n_interim[interim.index],
                            prior.mean=prior.mean.2,prior.sd = prior.sd.2,
                            prior.k = prior.k.2,prior.a=prior.a.2,prior.b=prior.b.2,
                            sd=sd.2,ssq=ssq.2,
                            CT1.go=CT1.go.2[interim.index],
                            false.go.CT1 = false.go.CT1.2[interim.index], FGR.CT1=FGR.CT1.2[interim.index],
                            CT1.nogo=CT1.nogo.2[interim.index],
                            false.nogo.CT1 = false.nogo.CT1.2[interim.index], FNGR.CT1 = FNGR.CT1.2[interim.index],
                            CT2.go=CT2.go.2[interim.index],
                            false.go.CT2 = false.go.CT2.2[interim.index], FGR.CT2=FGR.CT2.2[interim.index],
                            CT2.nogo=CT2.nogo.2[interim.index],
                            false.nogo.CT2 = false.nogo.CT2.2[interim.index], FNGR.CT2 = FNGR.CT2.2[interim.index],
                            direction=direction.2[interim.index],
                            logic.go = logic.go.2[interim.index],logic.nogo=logic.nogo.2[interim.index],
                            method=method,
                            fix.var=fix.var,noninfo=noninfo,
                            seed.num=seed.num,
                            stop.criterion=stop.criterion)


      true_go_cutoff2<-ifelse(temp2$overlap==0,temp2$cutoff[1],(overlap.option.2[interim.index]=='GO')*temp2$cutoff[1]+(overlap.option.2[interim.index]=='NOGO')*(temp2$cutoff[2]))
      true_nogo_cutoff2<-ifelse(temp2$overlap==0,temp2$cutoff[2],(overlap.option.2[interim.index]=='GO')*(temp2$cutoff[1])+(overlap.option.2[interim.index]=='NOGO')*(temp2$cutoff[2]))




    for(index.setting in 1: num_setting){


      rho=mean[index.setting,3]
      set.seed(seed.num)
      sigma<-matrix(NA,ncol=2,nrow=2)
      sigma[1,1]<-sd.1^2/n_interim[interim.index]
      sigma[2,2]<-sd.2^2/n_interim[interim.index]
      sigma[1,2]<-rho*sd.1*sd.2/n_interim[interim.index]
      sigma[2,1]<-rho*sd.1*sd.2/n_interim[interim.index]

      sim_data<-matrix(NA,nrow=nsim_IA,ncol=2)

      set.seed(seed.num)
      sim_mean=rmvnorm(nsim_IA,mean=mean[index.setting,1:2],sigma=sigma)
      sumx=sim_mean[,1] ####mean of x
      sumy=sim_mean[,2] ####mean of y



      if(logic.go.points[interim.index]=='and'){
        if(direction.1[interim.index]=='Greater'&direction.2[interim.index]=='Greater'){
          go_matrix[,interim.index,index.setting]<-(sumx>=true_go_cutoff1)&(sumy>=true_go_cutoff2)
          table_true_go_cutoff[interim.index,index.setting]<-HTML('X>=',round(true_go_cutoff1,3),logic.go.points[interim.index],' Y>=',round(true_go_cutoff2,3))
          table_go_cutoff[interim.index,index.setting]<-HTML('X>=',round(temp1$cutoff[1],3),logic.go.points[interim.index],' Y>=',round(temp2$cutoff[1],3))
        }
        if(direction.1[interim.index]=='Greater'&direction.2[interim.index]=='Less'){
          go_matrix[,interim.index,index.setting]<-(sumx>=true_go_cutoff1)&(sumy<=true_go_cutoff2)
          table_true_go_cutoff[interim.index,index.setting]<-HTML('X>=',round(true_go_cutoff1,3),logic.go.points[interim.index],' Y<=',round(true_go_cutoff2,3))
          table_go_cutoff[interim.index,index.setting]<-HTML('X>=',round(temp1$cutoff[1],3),logic.go.points[interim.index],' Y<=',round(temp2$cutoff[1],3))
        }
        if(direction.1[interim.index]=='Less'&direction.2[interim.index]=='Greater'){
          go_matrix[,interim.index,index.setting]<-(sumx<=true_go_cutoff1)&(sumy>=true_go_cutoff2)
          table_true_go_cutoff[interim.index,index.setting]<-HTML('X<=',round(true_go_cutoff1,3),logic.go.points[interim.index],' Y>=',round(true_go_cutoff2,3))
          table_go_cutoff[interim.index,index.setting]<-HTML('X<=',round(temp1$cutoff[1],3),logic.go.points[interim.index],' Y>=',round(temp2$cutoff[1],3))
        }
        if(direction.1[interim.index]=='Less'&direction.2[interim.index]=='Less'){
          go_matrix[,interim.index,index.setting]<-(sumx<=true_go_cutoff1)&(sumy<=true_go_cutoff2)
          table_true_go_cutoff[interim.index,index.setting]<-HTML('X<=',round(true_go_cutoff1,3),logic.go.points[interim.index],' Y<=',round(true_go_cutoff2,3))
          table_go_cutoff[interim.index,index.setting]<-HTML('X<=',round(temp1$cutoff[1],3),logic.go.points[interim.index],' Y<=',round(temp2$cutoff[1],3))

        }
      }
      if(logic.go.points[interim.index]=='or'){
        if(direction.1[interim.index]=='Greater'&direction.2[interim.index]=='Greater'){
          go_matrix[,interim.index,index.setting]<-(sumx>=true_go_cutoff1)|(sumy>=true_go_cutoff2)
          table_true_go_cutoff[interim.index,index.setting]<-HTML('X>=',round(true_go_cutoff1,3),logic.go.points[interim.index],' Y>=',round(true_go_cutoff2,3))
          table_go_cutoff[interim.index,index.setting]<-HTML('X>=',round(temp1$cutoff[1],3),logic.go.points[interim.index],' Y>=',round(temp2$cutoff[1],3))
        }
        if(direction.1[interim.index]=='Greater'&direction.2[interim.index]=='Less'){
          go_matrix[,interim.index,index.setting]<-(sumx>=true_go_cutoff1)|(sumy<=true_go_cutoff2)
          table_true_go_cutoff[interim.index,index.setting]<-HTML('X>=',round(true_go_cutoff1,3),logic.go.points[interim.index],' Y<=',round(true_go_cutoff2,3))
          table_go_cutoff[interim.index,index.setting]<-HTML('X>=',round(temp1$cutoff[1],3),logic.go.points[interim.index],' Y<=',round(temp2$cutoff[1],3))
        }
        if(direction.1[interim.index]=='Less'&direction.2[interim.index]=='Greater'){
          go_matrix[,interim.index,index.setting]<-(sumx<=true_go_cutoff1)|(sumy>=true_go_cutoff2)
          table_true_go_cutoff[interim.index,index.setting]<-HTML('X<=',round(true_go_cutoff1,3),logic.go.points[interim.index],' Y>=',round(true_go_cutoff2,3))
          table_go_cutoff[interim.index,index.setting]<-HTML('X<=',round(temp1$cutoff[1],3),logic.go.points[interim.index],' Y>=',round(temp2$cutoff[1],3))
        }
        if(direction.1[interim.index]=='Less'&direction.2[interim.index]=='Less'){
          go_matrix[,interim.index,index.setting]<-(sumx<=true_go_cutoff1)|(sumy<=true_go_cutoff2)
          table_true_go_cutoff[interim.index,index.setting]<-HTML('X<=',round(true_go_cutoff1,3),logic.go.points[interim.index],' Y<=',round(true_go_cutoff2,3))
          table_go_cutoff[interim.index,index.setting]<-HTML('X<=',round(temp1$cutoff[1],3),logic.go.points[interim.index],' Y<=',round(temp2$cutoff[1],3))

        }
      }

      if(logic.nogo.points[interim.index]=='and'){
        if(direction.1[interim.index]=='Greater'&direction.2[interim.index]=='Greater'){
          nogo_matrix[,interim.index,index.setting]<-(sumx<=true_nogo_cutoff1)&(sumy<=true_nogo_cutoff2)
          table_true_nogo_cutoff[interim.index,index.setting]<-HTML('X<',round(true_nogo_cutoff1,3),logic.nogo.points[interim.index],' Y<',round(true_nogo_cutoff2,3))
          table_nogo_cutoff[interim.index,index.setting]<-HTML('X<',round(temp1$cutoff[2],3),logic.nogo.points[interim.index],' Y<',round(temp2$cutoff[2],3))

        }
        if(direction.1[interim.index]=='Greater'&direction.2[interim.index]=='Less'){
          nogo_matrix[,interim.index,index.setting]<-(sumx<=true_nogo_cutoff1)&(sumy>=true_nogo_cutoff2)
          table_true_nogo_cutoff[interim.index,index.setting]<-HTML('X<',round(true_nogo_cutoff1,3),logic.nogo.points[interim.index],' Y>',round(true_nogo_cutoff2,3))
          table_nogo_cutoff[interim.index,index.setting]<-HTML('X<',round(temp1$cutoff[2],3),logic.nogo.points[interim.index],' Y>',round(temp2$cutoff[2],3))

        }
        if(direction.1[interim.index]=='Less'&direction.2[interim.index]=='Greater'){
          nogo_matrix[,interim.index,index.setting]<-(sumx>=true_nogo_cutoff1)&(sumy<=true_nogo_cutoff2)
          table_true_nogo_cutoff[interim.index,index.setting]<-HTML('X>',round(true_nogo_cutoff1,3),logic.nogo.points[interim.index],' Y<',round(true_nogo_cutoff2,3))
          table_nogo_cutoff[interim.index,index.setting]<-HTML('X>',round(temp1$cutoff[2],3),logic.nogo.points[interim.index],' Y<',round(temp2$cutoff[2],3))

        }
        if(direction.1[interim.index]=='Less'&direction.2[interim.index]=='Less'){
          nogo_matrix[,interim.index,index.setting]<-(sumx>=true_nogo_cutoff1)&(sumy>=true_nogo_cutoff2)
          table_true_nogo_cutoff[interim.index,index.setting]<-HTML('X>',round(true_nogo_cutoff1,3),logic.nogo.points[interim.index],' Y>',round(true_nogo_cutoff2,3))
          table_nogo_cutoff[interim.index,index.setting]<-HTML('X>',round(temp1$cutoff[2],3),logic.nogo.points[interim.index],' Y>',round(temp2$cutoff[2],3))

        }
      }

      if(logic.nogo.points[interim.index]=='or'){
        if(direction.1[interim.index]=='Greater'&direction.2[interim.index]=='Greater'){
          nogo_matrix[,interim.index,index.setting]<-(sumx<=true_nogo_cutoff1)|(sumy<=true_nogo_cutoff2)
          table_true_nogo_cutoff[interim.index,index.setting]<-HTML('X<',round(true_nogo_cutoff1,3),logic.nogo.points[interim.index],' Y<',round(true_nogo_cutoff2,3))
          table_nogo_cutoff[interim.index,index.setting]<-HTML('X<',round(temp1$cutoff[2],3),logic.nogo.points[interim.index],' Y<',round(temp2$cutoff[2],3))

        }
        if(direction.1[interim.index]=='Greater'&direction.2[interim.index]=='Less'){
          nogo_matrix[,interim.index,index.setting]<-(sumx<=true_nogo_cutoff1)|(sumy>=true_nogo_cutoff2)
          table_true_nogo_cutoff[interim.index,index.setting]<-HTML('X<',round(true_nogo_cutoff1,3),logic.nogo.points[interim.index],' Y>',round(true_nogo_cutoff2,3))
          table_nogo_cutoff[interim.index,index.setting]<-HTML('X<',round(temp1$cutoff[2],3),logic.nogo.points[interim.index],' Y>',round(temp2$cutoff[2],3))

        }
        if(direction.1[interim.index]=='Less'&direction.2[interim.index]=='Greater'){
          nogo_matrix[,interim.index,index.setting]<-(sumx>=true_nogo_cutoff1)|(sumy<=true_nogo_cutoff2)
          table_true_nogo_cutoff[interim.index,index.setting]<-HTML('X>',round(true_nogo_cutoff1,3),logic.nogo.points[interim.index],' Y<',round(true_nogo_cutoff2,3))
          table_nogo_cutoff[interim.index,index.setting]<-HTML('X>',round(temp1$cutoff[2],3),logic.nogo.points[interim.index],' Y<',round(temp2$cutoff[2],3))

        }
        if(direction.1[interim.index]=='Less'&direction.2[interim.index]=='Less'){
          nogo_matrix[,interim.index,index.setting]<-(sumx>=true_nogo_cutoff1)|(sumy>=true_nogo_cutoff2)
          table_true_nogo_cutoff[interim.index,index.setting]<-HTML('X>',round(true_nogo_cutoff1,3),logic.nogo.points[interim.index],' Y>',round(true_nogo_cutoff2,3))
          table_nogo_cutoff[interim.index,index.setting]<-HTML('X>',round(temp1$cutoff[2],3),logic.nogo.points[interim.index],' Y>',round(temp2$cutoff[2],3))

        }
      }



      inconclusive_matrix[,interim.index,index.setting]<-1-go_matrix[,interim.index,index.setting]-nogo_matrix[,interim.index,index.setting]

      if(task[interim.index]=='Futility'){
        IA_go_matrix[,interim.index,index.setting]=inconclusive_matrix[,interim.index,index.setting]+go_matrix[,interim.index,index.setting]
        table_zone[interim.index,index.setting]<-paste0('-/',table_true_nogo_cutoff[interim.index,index.setting])
      }
      if(task[interim.index]=='Superiority'){
        IA_go_matrix[,interim.index,index.setting]=inconclusive_matrix[,interim.index,index.setting]+nogo_matrix[,interim.index,index.setting]
        table_zone[interim.index,index.setting]<-paste0(table_true_go_cutoff[interim.index,index.setting],'/-')
      }
      if(task[interim.index]=='Futility and superiority'){
        IA_go_matrix[,interim.index,index.setting]=inconclusive_matrix[,interim.index,index.setting]
        table_zone[interim.index,index.setting]<-paste0(table_true_go_cutoff[interim.index,index.setting],'/',table_true_nogo_cutoff[interim.index,index.setting])
      }


      table_RR[interim.index,index.setting]<-paste0(mean[index.setting,],collapse = ',')
      table_warning[interim.index,index.setting]<-'None'
      if(temp1$overlap==1&temp2$overlap==1){
        table_warning[interim.index,index.setting]<-HTML('GO/NOGO zones for Endpoint 1 are overlapped, classifed by ',overlap.option.1[interim.index],' zone.<br>',
                                                         'GO/NOGO zones for Endpoint 2 are overlapped, classifed by ',overlap.option.2[interim.index],' zone.')
      }
      if(temp1$overlap==1&temp2$overlap==0){
        table_warning[interim.index,index.setting]<-HTML('GO/NOGO zones for Endpoint 1 are overlapped, classifed by ',overlap.option.1[interim.index],' zone')
      }
      if(temp1$overlap==0&temp2$overlap==1){
        table_warning[interim.index,index.setting]<-HTML('GO/NOGO zones for Endpoint 2 are overlapped, classifed by ',overlap.option.2[interim.index],' zone')
      }
      table_n[interim.index,index.setting]<-n_interim[interim.index]

    }
  }


  temptable=c()
  for(meanindex in 1:num_setting){
    table<-matrix(NA,ncol=num_interim+1,nrow=10)
    cum_IA_go_matrix<-t(apply(IA_go_matrix[,,meanindex],1,cumprod))
    for(j in 1:(num_interim)){
      table[1,j]=table_n[j,meanindex]
      table[2,j]=paste0(mean[meanindex,],collapse = ',')
      table[3,j]=task[j]
      if(j==1){
        if(task[j]=='Superiority'|task[j]=='Futility and superiority'){
          table[4,j]=round(sum(go_matrix[,j,meanindex]==1)/nsim_IA,3)}else{table[4,j]=0}
        if(task[j]=='Futility'|task[j]=='Futility and superiority'){
          table[6,j]=round(sum(nogo_matrix[,j,meanindex]==1)/nsim_IA,3)
        }else{table[6,j]=0}
      }else{
        if(task[j]=='Superiority'|task[j]=='Futility and superiority'){
          table[4,j]=round(sum(go_matrix[,j,meanindex]==1&cum_IA_go_matrix[,j-1]==1)/nsim_IA,3)}else{table[4,j]=0}
        if(task[j]=='Futility'|task[j]=='Futility and superiority'){
          table[6,j]=round(sum(nogo_matrix[,j,meanindex]==1&cum_IA_go_matrix[,j-1]==1)/nsim_IA,3)
        }else{table[6,j]=0}
      }
      table[5,j]=round(sum(cum_IA_go_matrix[,j]==1)/nsim_IA,3)
      table[7,j]=table_zone[j,meanindex]
      table[10,j]=table_go_cutoff[j,meanindex]
      table[9,j]=table_nogo_cutoff[j,meanindex]
      table[8,j]=table_warning[j,meanindex]

    }
    expectss<-round(sum(as.numeric(table[1,1:num_interim])*(c(as.numeric(table[4,1:num_interim-1])+as.numeric(table[6,1:num_interim-1]),as.numeric(table[5,num_interim-1])))),3)
    table[1,num_interim+1]=HTML(paste0(expectss,' (expected)'))
    table[2,num_interim+1]=paste0(mean[meanindex,],collapse = ',')
    table[3,num_interim+1]=''
    table[4,num_interim+1]=round(sum(as.numeric(table[4,1:num_interim])),3)
    table[5,num_interim+1]=round(as.numeric(table[5,num_interim]),3)
    table[6,num_interim+1]=round(sum(as.numeric(table[6,1:num_interim])),3)
    table[7,num_interim+1]=''
    table[8,num_interim+1]=''
    table[9,num_interim+1]=''
    table[10,num_interim+1]=''
    tablecolname<-c(paste0('Interim analysis ',1:(num_interim-1)),'Final analysis',"Summary")
    tablerowname<-c('Sample size',HTML('&mu;<sub>1</sub>,&mu;<sub>2</sub>,Corr'),'Task','Success','To next interim/final or inconclusive',
                    'Stop','Superority/Futility zone','Warning',
                    'Cut off for NOGO rule',
                    'Cut off for GO rule')
    table<-cbind(tablerowname,rep(meanindex,10),table)
    colnames(table)<-c(" ",'Setting',tablecolname)
    temptable=rbind(temptable,table)
  }
  return(temptable)
}

#CN_OC_table()
#CN_OC_table_interim()
