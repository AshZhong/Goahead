library(shiny)

# pmf of multinomial distribution for 2x2 contingency table
CB_prob<-function(x,prob){
  a00=prob[1]
  a01=prob[2]
  a10=prob[3]
  a11=prob[4]
  n=sum(x)
  n1=n-x[1]
  n2=n1-x[2]
  return(choose(n,x[1])*choose(n1,x[2])*choose(n2,x[3])*(a00^x[1])*(a01^x[2])*(a10^x[3])*(a11^x[4]))
}


CB_OC_table<-function(n=c(100,200),
                      a.1=1/3,b.1=1,
                      CT1.go.1=0.25,
                      false.go.CT1.1=TRUE,FGR.CT1.1=0.25,
                      CT1.nogo.1=0.25,
                      false.nogo.CT1.1=TRUE,FNGR.CT1.1=0.25,
                      CT2.go.1=0.3,
                      false.go.CT2.1=TRUE, FGR.CT2.1=0.5,
                      CT2.nogo.1=0.3,
                      false.nogo.CT2.1=TRUE, FNGR.CT2.1=0.5,
                      direction.1='Greater',
                      logic.go.1='and',logic.nogo.1='and',
                      overlap.option.1='GO',
                      a.2=1/3,b.2=1,
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
                      RR=rbind(c(0.2,0.3,0.4,0.1),c(0.3,0.3,0,0.4)),
                      method='Bayesian'){

  ############Arguments#########
  ###Endpoint 1###
  #n: sample size
  #CT1.go.1: Criterion value 1 in GO rule
  #false.go.CT1.1:'TRUE' or 'FALSE', whether include CT1.go.1 in GO rule
  #FGR.CT1.1: False GO risk (1-confidence of GO for CT1.go.1)
  #CT2.go.1: Criterion value 2 in G rule
  #false.go.CT2.1:'TRUE' or 'FALSE', whether include CT2.go.1 in GO rule
  #FGR.CT2.1: False GO risk (1-confidence of GO for CT2.go.1)
  ###logic.go.1:'and' or 'or' relationship between two GO criterions for Endpoint 1

  #CT1.nogo.1: Criterion value 1 in NOGO rule
  #false.nogo.CT1.1:'TRUE' or 'FALSE', whether include CT1.nogo.1 in NOGO rule
  #FNGR.CT1.1: False NOGO risk for CT1.nogo.1
  #CT2.nogo.1: Criterion value 2 in NOGO rule
  #false.nogo.CT2.1:'TRUE' or 'FALSE', whether include CT2.nogo.1 in NOGO rule
  #FNGR.CT2.1: False NOGO risk FOR CT2.go.1
  ###logic.nogo.1:'and' or 'or' relationship between two NOGO criterions for Endpoint 1

  #overlap.option.1:'GO' or 'NOGO', the dominated rule for Endpoint 1
  #direction.1:'Greater' or 'Less', the direciton of comparision for Endpoint 1

  ###Endpoint 2###
  #CT1.go.2: Criterion value 1 in GO rule
  #false.go.CT1.2:'TRUE' or 'FALSE', whether include CT1.go.1 in GO rule
  #FGR.CT1.2: False GO risk (1-confidence of GO for CT1.go.1)
  #CT2.go.2: Criterion value 2 in G rule
  #false.go.CT2.2:'TRUE' or 'FALSE', whether include CT2.go.1 in GO rule
  #FGR.CT2.2: False GO risk (1-confidence of GO for CT2.go.1)
  ###logic.go.2:'and' or 'or' relationship between two GO criterions for endpoint 2

  #CT1.nogo.2: Criterion value 1 in NOGO rule
  #false.nogo.CT1.2:'TRUE' or 'FALSE', whether include CT1.nogo.1 in NOGO rule
  #FNGR.CT1.2: False NOGO risk for CT1.nogo.1
  #CT2.nogo.2: Criterion value 2 in NOGO rule
  #false.nogo.CT2.2:'TRUE' or 'FALSE', whether include CT2.nogo.1 in NOGO rule
  #FNGR.CT2.2: False NOGO risk FOR CT2.go.1
  ###logic.nogo.2:'and' or 'or' relationship between two NOGO criterions for Endpoint 2
  #overlap.option.2:'GO' or 'NOGO', the dominated rule for Endpoint 2
  #direction.2:'Greater' or 'Less', the direciton of comparision for Endpoint 2

  ###
  #a.1,b.1: beta prior for RR in endpoint 1
  #a.2,b.2: beta prior for RR in endpoint 2
  #logic.go.points:'and' or 'or', the relationship between GO criterions for two Endpoints,
  #logic.nogo.points:'and' or 'or', the relationship between GO criterions for two Endpoints
  #RR: matrix, the underlying true setting. each row is a setting, includes 4 elements: (PR(X=0,Y=0),PR(X=0,Y=1),PR(X=1,Y=0),PR(X=1,Y=1))

  #method:"Bayesian" or "Frequentist", the inference method

  ###Return: OC table
  ###########End of Arguments#########################
  temp_table<-c()
  for(jj in 1:length(n)){
    temp1=SAB_Bin_Cutoff(n=n[jj],a=a.1,b=b.1,
                         CT1.go=CT1.go.1,
                         false.go.CT1 = false.go.CT1.1, FGR.CT1=FGR.CT1.1,
                         CT1.nogo=CT1.nogo.1,
                         false.nogo.CT1 = false.nogo.CT1.1, FNGR.CT1 = FNGR.CT1.1,
                         CT2.go=CT2.go.1,
                         false.go.CT2 = false.go.CT2.1, FGR.CT2=FGR.CT2.1,
                         CT2.nogo=CT2.nogo.1,
                         false.nogo.CT2 = false.nogo.CT2.1, FNGR.CT2 = FNGR.CT2.1,
                         direction=direction.1,logic.go = logic.go.1,logic.nogo=logic.nogo.1,
                         method=method)
    if(direction.1=='Greater'){
      true_go_cutoff1<-ifelse(temp1$overlap==0,temp1$cutoff[1],(overlap.option.1=='GO')*temp1$cutoff[1]+(overlap.option.1=='NOGO')*(temp1$cutoff[2]))
      true_nogo_cutoff1<-ifelse(temp1$overlap==0,temp1$cutoff[2],(overlap.option.1=='GO')*(temp1$cutoff[1])+(overlap.option.1=='NOGO')*(temp1$cutoff[2]))
    }
    if(direction.1=='Less'){
      true_go_cutoff1<-ifelse(temp1$overlap==0,temp1$cutoff[1],(overlap.option.1=='GO')*temp1$cutoff[1]+(overlap.option.1=='NOGO')*(temp1$cutoff[2]))
      true_nogo_cutoff1<-ifelse(temp1$overlap==0,temp1$cutoff[2],(overlap.option.1=='GO')*(temp1$cutoff[1])+(overlap.option.1=='NOGO')*(temp1$cutoff[2]))
    }
  
  
    temp2=SAB_Bin_Cutoff(n=n[jj],a=a.2,b=b.2,
                         CT1.go=CT1.go.2,
                         false.go.CT1 = false.go.CT1.2, FGR.CT1=FGR.CT1.2,
                         CT1.nogo=CT1.nogo.2,
                         false.nogo.CT1 = false.nogo.CT1.2, FNGR.CT1 = FNGR.CT1.2,
                         CT2.go=CT2.go.2,
                         false.go.CT2 = false.go.CT2.2, FGR.CT2=FGR.CT2.2,
                         CT2.nogo=CT2.nogo.2,
                         false.nogo.CT2 = false.nogo.CT2.2, FNGR.CT2 = FNGR.CT2.2,
                         direction=direction.2,logic.go = logic.go.2,logic.nogo=logic.nogo.2,
                         method=method)
  
    if(direction.2=='Greater'){
      true_go_cutoff2<-ifelse(temp2$overlap==0,temp2$cutoff[1],(overlap.option.2=='GO')*temp2$cutoff[1]+(overlap.option.2=='NOGO')*(temp2$cutoff[2]))
      true_nogo_cutoff2<-ifelse(temp2$overlap==0,temp2$cutoff[2],(overlap.option.2=='GO')*(temp2$cutoff[1])+(overlap.option.2=='NOGO')*(temp2$cutoff[2]))
    }
    if(direction.2=='Less'){
      true_go_cutoff2<-ifelse(temp2$overlap==0,temp2$cutoff[1],(overlap.option.2=='GO')*temp2$cutoff[1]+(overlap.option.2=='NOGO')*(temp2$cutoff[2]))
      true_nogo_cutoff2<-ifelse(temp2$overlap==0,temp2$cutoff[2],(overlap.option.2=='GO')*(temp2$cutoff[1])+(overlap.option.2=='NOGO')*(temp2$cutoff[2]))
    }
  
    num_combination=choose(n[jj]+3,3) # total number of possible outcomes
    sumx=rep(NA,num_combination)
    sumy=rep(NA,num_combination)
    RR<-as.matrix(RR)
    num_setting<-dim(RR)[1] # number of interested settings
  
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
      a.00=RR[index.setting,1]
      a.01=RR[index.setting,2]
      a.10=RR[index.setting,3]
      a.11=RR[index.setting,4]
      prob_four=c(a.00,a.10,a.01,a.11)
      if((any(prob_four<0))|(any(prob_four>1))){
        stop(paste0('The setting ',index.setting, ' is not reasonable, please check it'))
      }
      k=1
      prob_table<-rep(NA,num_combination)
      for(i in 0:n[jj]){
        for(j in 0:(n[jj]-i)){
          for(q in 0:(n[jj]-i-j)){
            sumx[k]<-n[jj]-i-j
            sumy[k]<-n[jj]-i-q
            prob_table[k]=CB_prob(c(i,j,q,n[jj]-i-j-q),RR[index.setting,])
            k=k+1
          }
        }
      }
  
      if(logic.go.points=='and'){
        if(direction.1=='Greater'&direction.2=='Greater'){
          go_prob[index.setting]<-sum(prob_table[(sumx>=true_go_cutoff1)&(sumy>=true_go_cutoff2)])
          table_true_go_cutoff[index.setting]<-HTML('X>=',true_go_cutoff1,logic.go.points,' Y>=',true_go_cutoff2)
          table_go_cutoff[index.setting]<-HTML('X>=',temp1$cutoff[1],logic.go.points,' Y>=',temp2$cutoff[1])
        }
        if(direction.1=='Greater'&direction.2=='Less'){
          go_prob[index.setting]<-sum(prob_table[(sumx>=true_go_cutoff1)&(sumy<=true_go_cutoff2)])
          table_true_go_cutoff[index.setting]<-HTML('X>=',true_go_cutoff1,logic.go.points,' Y<=',true_go_cutoff2)
          table_go_cutoff[index.setting]<-HTML('X>=',temp1$cutoff[1],logic.go.points,' Y<=',temp2$cutoff[1])
        }
        if(direction.1=='Less'&direction.2=='Greater'){
          go_prob[index.setting]<-sum(prob_table[(sumx<=true_go_cutoff1)&(sumy>=true_go_cutoff2)])
          table_true_go_cutoff[index.setting]<-HTML('X<=',true_go_cutoff1,logic.go.points,' Y>=',true_go_cutoff2)
          table_go_cutoff[index.setting]<-HTML('X<=',temp1$cutoff[1],logic.go.points,' Y>=',temp2$cutoff[1])
        }
        if(direction.1=='Less'&direction.2=='Less'){
          go_prob[index.setting]<-sum(prob_table[(sumx<=true_go_cutoff1)&(sumy<=true_go_cutoff2)])
          table_true_go_cutoff[index.setting]<-HTML('X<=',true_go_cutoff1,logic.go.points,' Y<=',true_go_cutoff2)
          table_go_cutoff[index.setting]<-HTML('X<=',temp1$cutoff[1],logic.go.points,' Y<=',temp2$cutoff[1])
  
        }
      }
      if(logic.go.points=='or'){
        if(direction.1=='Greater'&direction.2=='Greater'){
          go_prob[index.setting]<-sum(prob_table[(sumx>=true_go_cutoff1)|(sumy>=true_go_cutoff2)])
          table_true_go_cutoff[index.setting]<-HTML('X>=',true_go_cutoff1,logic.go.points,' Y>=',true_go_cutoff2)
          table_go_cutoff[index.setting]<-HTML('X>=',temp1$cutoff[1],logic.go.points,' Y>=',temp2$cutoff[1])
  
        }
        if(direction.1=='Greater'&direction.2=='Less'){
          go_prob[index.setting]<-sum(prob_table[(sumx>=true_go_cutoff1)|(sumy<=true_go_cutoff2)])
          table_true_go_cutoff[index.setting]<-HTML('X>=',true_go_cutoff1,logic.go.points,' Y<=',true_go_cutoff2)
          table_go_cutoff[index.setting]<-HTML('X>=',temp1$cutoff[1],logic.go.points,' Y<=',temp2$cutoff[1])
        }
        if(direction.1=='Less'&direction.2=='Greater'){
          go_prob[index.setting]<-sum(prob_table[(sumx<=true_go_cutoff1)|(sumy>=true_go_cutoff2)])
          table_true_go_cutoff[index.setting]<-HTML('X<=',true_go_cutoff1,logic.go.points,' Y>=',true_go_cutoff2)
          table_go_cutoff[index.setting]<-HTML('X<=',temp1$cutoff[1],logic.go.points,' Y>=',temp2$cutoff[1])
        }
        if(direction.1=='Less'&direction.2=='Less'){
          go_prob[index.setting]<-sum(prob_table[(sumx<=true_go_cutoff1)|(sumy<=true_go_cutoff2)])
          table_true_go_cutoff[index.setting]<-HTML('X<=',true_go_cutoff1,logic.go.points,' Y<=',true_go_cutoff2)
          table_go_cutoff[index.setting]<-HTML('X<=',temp1$cutoff[1],logic.go.points,' Y<=',temp2$cutoff[1])
        }
      }
  
      if(logic.nogo.points=='and'){
        if(direction.1=='Greater'&direction.2=='Greater'){
          nogo_prob[index.setting]<-sum(prob_table[(sumx<true_nogo_cutoff1)&(sumy<true_nogo_cutoff2)])
          table_true_nogo_cutoff[index.setting]<-HTML('X<',true_nogo_cutoff1,logic.nogo.points,' Y<',true_nogo_cutoff2)
          table_nogo_cutoff[index.setting]<-HTML('X<',temp1$cutoff[2],logic.nogo.points,' Y<',temp2$cutoff[2])
  
        }
        if(direction.1=='Greater'&direction.2=='Less'){
          nogo_prob[index.setting]<-sum(prob_table[(sumx<true_nogo_cutoff1)&(sumy>true_nogo_cutoff2)])
          table_true_nogo_cutoff[index.setting]<-HTML('X<',true_nogo_cutoff1,logic.nogo.points,' Y>',true_nogo_cutoff2)
          table_nogo_cutoff[index.setting]<-HTML('X<',temp1$cutoff[2],logic.nogo.points,' Y>',temp2$cutoff[2])
  
        }
        if(direction.1=='Less'&direction.2=='Greater'){
          nogo_prob[index.setting]<-sum(prob_table[(sumx>true_nogo_cutoff1)&(sumy<true_nogo_cutoff2)])
          table_true_nogo_cutoff[index.setting]<-HTML('X>',true_nogo_cutoff1,logic.nogo.points,' Y<',true_nogo_cutoff2)
          table_nogo_cutoff[index.setting]<-HTML('X>',temp1$cutoff[2],logic.nogo.points,' Y<',temp2$cutoff[2])
  
        }
        if(direction.1=='Less'&direction.2=='Less'){
          nogo_prob[index.setting]<-sum(prob_table[(sumx>true_nogo_cutoff1)&(sumy>true_nogo_cutoff2)])
          table_true_nogo_cutoff[index.setting]<-HTML('X<',true_nogo_cutoff1,logic.nogo.points,' Y>',true_nogo_cutoff2)
          table_nogo_cutoff[index.setting]<-HTML('X<',temp1$cutoff[2],logic.nogo.points,' Y>',temp2$cutoff[2])
  
        }
      }
  
      if(logic.nogo.points=='or'){
        if(direction.1=='Greater'&direction.2=='Greater'){
          nogo_prob[index.setting]<-sum(prob_table[(sumx<true_nogo_cutoff1)|(sumy<true_nogo_cutoff2)])
          table_true_nogo_cutoff[index.setting]<-HTML('X<',true_nogo_cutoff1,logic.nogo.points,' Y<',true_nogo_cutoff2)
          table_nogo_cutoff[index.setting]<-HTML('X<',temp1$cutoff[2],logic.nogo.points,' Y<',temp2$cutoff[2])
  
        }
        if(direction.1=='Greater'&direction.2=='Less'){
          nogo_prob[index.setting]<-sum(prob_table[(sumx<true_nogo_cutoff1)|(sumy>true_nogo_cutoff2)])
          table_true_nogo_cutoff[index.setting]<-HTML('X<',true_nogo_cutoff1,logic.nogo.points,' Y>',true_nogo_cutoff2)
          table_nogo_cutoff[index.setting]<-HTML('X<',temp1$cutoff[2],logic.nogo.points,' Y>',temp2$cutoff[2])
  
  
        }
        if(direction.1=='Less'&direction.2=='Greater'){
          nogo_prob[index.setting]<-sum(prob_table[(sumx>true_nogo_cutoff1)|(sumy<true_nogo_cutoff2)])
          table_true_nogo_cutoff[index.setting]<-HTML('X>',true_nogo_cutoff1,logic.nogo.points,' Y<',true_nogo_cutoff2)
          table_nogo_cutoff[index.setting]<-HTML('X>',temp1$cutoff[2],logic.nogo.points,' Y<',temp2$cutoff[2])
  
  
        }
        if(direction.1=='Less'&direction.2=='Less'){
          nogo_prob[index.setting]<-sum(prob_table[(sumx>true_nogo_cutoff1)|(sumy>true_nogo_cutoff2)])
          table_true_nogo_cutoff[index.setting]<-HTML('X>',true_nogo_cutoff1,logic.nogo.points,' Y>',true_nogo_cutoff2)
          table_nogo_cutoff[index.setting]<-HTML('X>',temp1$cutoff[2],logic.nogo.points,' Y>',temp2$cutoff[2])
  
        }
      }
  
      inconclusive_prob[index.setting]<-1-go_prob[index.setting]-nogo_prob[index.setting]
  
      table_RR[index.setting]<-paste0(round(RR[index.setting,],3),collapse = ',')
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
    table_n<-rep(n[jj],num_setting)
    table_go_prob<-round(go_prob,3)
    table_nogo_prob<-round(nogo_prob,3)
    table_inconclusive_prob<-round(inconclusive_prob,3)
    table_zone<-paste0(table_true_go_cutoff,'/',table_true_nogo_cutoff)
    temp<-cbind(table_n,table_RR,table_go_prob,table_nogo_prob,table_inconclusive_prob,table_zone,table_warning,table_go_cutoff,table_nogo_cutoff)
    temp_table<-rbind(temp_table,temp)
  }
  colnames(temp_table)<-c('Sample size','P<sub>00</sub>,P<sub>01</sub>,P<sub>10</sub>,P<sub>11</sub>','Pr(GO)','Pr(NOGO)','Pr(Inconclusive)','GO/NOGO zones','Warning','Cut off for GO','Cut off for NOGO')
  return(temp_table)
}


CB_OC_table_interim<-function(num_interim=3,n_interim=c(100,200,300),
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
                              a.1=1/3,b.1=1,
                              a.2=1/3,b.2=1,
                              RR=rbind(c(0.2,0.3,0.4,0.1),c(0,0,0.5,0.5)),
                              method='Bayesian',nsim_IA=10000,seed.num=425){

  ############Arguments#########
  #num_interim:number of interim analysis and final analysis
  #n_interim: sample size in each interim analysis and final analysis
  ###Endpoint 1###
  #CT1.go.1: Criterion value 1 in GO rule
  #false.go.CT1.1:'TRUE' or 'FALSE', whether include CT1.go.1 in GO rule
  #FGR.CT1.1: False GO risk (1-confidence of GO for CT1.go.1)
  #CT2.go.1: Criterion value 2 in G rule
  #false.go.CT2.1:'TRUE' or 'FALSE', whether include CT2.go.1 in GO rule
  #FGR.CT2.1: False GO risk (1-confidence of GO for CT2.go.1)
  ###logic.go.1:'and' or 'or' relationship between two GO criterions for Endpoint 1

  #CT1.nogo.1: Criterion value 1 in NOGO rule
  #false.nogo.CT1.1:'TRUE' or 'FALSE', whether include CT1.nogo.1 in NOGO rule
  #FNGR.CT1.1: False NOGO risk for CT1.nogo.1
  #CT2.nogo.1: Criterion value 2 in NOGO rule
  #false.nogo.CT2.1:'TRUE' or 'FALSE', whether include CT2.nogo.1 in NOGO rule
  #FNGR.CT2.1: False NOGO risk FOR CT2.go.1
  ###logic.nogo.1:'and' or 'or' relationship between two NOGO criterions for Endpoint 1

  #overlap.option.1:'GO' or 'NOGO', the dominated rule for Endpoint 1
  #direction.1:'Greater' or 'Less', the direciton of comparision for Endpoint 1

  ###Endpoint 2###
  #CT1.go.2: Criterion value 1 in GO rule
  #false.go.CT1.2:'TRUE' or 'FALSE', whether include CT1.go.1 in GO rule
  #FGR.CT1.2: False GO risk (1-confidence of GO for CT1.go.1)
  #CT2.go.2: Criterion value 2 in G rule
  #false.go.CT2.2:'TRUE' or 'FALSE', whether include CT2.go.1 in GO rule
  #FGR.CT2.2: False GO risk (1-confidence of GO for CT2.go.1)
  ###logic.go.2:'and' or 'or' relationship between two GO criterions for endpoint 2

  #CT1.nogo.2: Criterion value 1 in NOGO rule
  #false.nogo.CT1.2:'TRUE' or 'FALSE', whether include CT1.nogo.1 in NOGO rule
  #FNGR.CT1.2: False NOGO risk for CT1.nogo.1
  #CT2.nogo.2: Criterion value 2 in NOGO rule
  #false.nogo.CT2.2:'TRUE' or 'FALSE', whether include CT2.nogo.1 in NOGO rule
  #FNGR.CT2.2: False NOGO risk FOR CT2.go.1
  ###logic.nogo.2:'and' or 'or' relationship between two NOGO criterions for Endpoint 2
  #overlap.option.2:'GO' or 'NOGO', the dominated rule for Endpoint 2
  #direction.2:'Greater' or 'Less', the direciton of comparision for Endpoint 2

  ###
  #logic.go.points:'and' or 'or', the relationship between GO criterions for two Endpoints,
  #logic.nogo.points:'and' or 'or', the relationship between GO criterions for two Endpoints
  #RR: matrix, the underlying true setting. each row is a setting, includes 4 elements: (PR(X=0,Y=0),PR(X=0,Y=1),PR(X=1,Y=0),PR(X=1,Y=1))
    #method:"Bayesian" or "Frequentist", the inference method

  ###Return: OC table
  ###########End of Arguments#########################
  go_matrix<-array(NA,c(nsim_IA,num_interim,dim(RR)[1]))
  nogo_matrix<-array(NA,c(nsim_IA,num_interim,dim(RR)[1]))
  inconclusive_matrix<-array(NA,c(nsim_IA,num_interim,dim(RR)[1]))
  IA_go_matrix<-array(NA,c(nsim_IA,num_interim,dim(RR)[1]))
  n_interim<-sort(n_interim)
  num_setting<-dim(RR)[1]
  table_RR<-array(NA,c(num_interim,num_setting))
    table_true_go_cutoff<-array(NA,c(num_interim,num_setting))
    table_true_nogo_cutoff<-array(NA,c(num_interim,num_setting))
    table_go_cutoff<-array(NA,c(num_interim,num_setting))
    table_nogo_cutoff<-array(NA,c(num_interim,num_setting))
    table_warning<-array(NA,c(num_interim,num_setting))
    table_zone<-array(NA,c(num_interim,num_setting))
    table_n<-array(NA,c(num_interim,num_setting))
  for(interim.index in 1:num_interim){
    temp1=SAB_Bin_Cutoff(n=n_interim[interim.index],a=a.1,b=b.1,
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
                         method=method)
    if(direction.1[interim.index]=='Greater'){
      true_go_cutoff1<-ifelse(temp1$overlap==0,temp1$cutoff[1],(overlap.option.1[interim.index]=='GO')*temp1$cutoff[1]+(overlap.option.1[interim.index]=='NOGO')*(temp1$cutoff[2]))
      true_nogo_cutoff1<-ifelse(temp1$overlap==0,temp1$cutoff[2],(overlap.option.1[interim.index]=='GO')*(temp1$cutoff[1])+(overlap.option.1[interim.index]=='NOGO')*(temp1$cutoff[2]))
    }
    if(direction.1[interim.index]=='Less'){
      true_go_cutoff1<-ifelse(temp1$overlap==0,temp1$cutoff[1],(overlap.option.1[interim.index]=='GO')*temp1$cutoff[1]+(overlap.option.1[interim.index]=='NOGO')*(temp1$cutoff[2]))
      true_nogo_cutoff1<-ifelse(temp1$overlap==0,temp1$cutoff[2],(overlap.option.1[interim.index]=='GO')*(temp1$cutoff[1])+(overlap.option.1[interim.index]=='NOGO')*(temp1$cutoff[2]))
    }


    temp2=SAB_Bin_Cutoff(n=n_interim[interim.index],a=a.2,b=b.2,
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
                         method=method)

    if(direction.2[interim.index]=='Greater'){
      true_go_cutoff2<-ifelse(temp2$overlap==0,temp2$cutoff[1],(overlap.option.2[interim.index]=='GO')*temp2$cutoff[1]+(overlap.option.2[interim.index]=='NOGO')*(temp2$cutoff[2]))
      true_nogo_cutoff2<-ifelse(temp2$overlap==0,temp2$cutoff[2],(overlap.option.2[interim.index]=='GO')*(temp2$cutoff[1])+(overlap.option.2[interim.index]=='NOGO')*(temp2$cutoff[2]))
    }
    if(direction.2[interim.index]=='Less'){
      true_go_cutoff2<-ifelse(temp2$overlap==0,temp2$cutoff[1],(overlap.option.2[interim.index]=='GO')*temp2$cutoff[1]+(overlap.option.2[interim.index]=='NOGO')*(temp2$cutoff[2]))
      true_nogo_cutoff2<-ifelse(temp2$overlap==0,temp2$cutoff[2],(overlap.option.2[interim.index]=='GO')*(temp2$cutoff[1])+(overlap.option.2[interim.index]=='NOGO')*(temp2$cutoff[2]))
    }



    for(index.setting in 1: num_setting){

      a.00=RR[index.setting,1]
      a.01=RR[index.setting,2]
      a.10=RR[index.setting,3]
      a.11=RR[index.setting,4]
      prob_four=c(a.00,a.01,a.10,a.11)
      if((any(prob_four<0))|(any(prob_four>1))){
        stop(paste0('The setting ',index.setting, ' is not reasonable, please check it'))
      }
      set.seed(seed.num)
      sim_data=rmultinom(nsim_IA,size=n_interim[interim.index],prob=prob_four)
      sumx=sim_data[3,]+sim_data[4,]
      sumy=sim_data[2,]+sim_data[4,]


      if(logic.go.points[interim.index]=='and'){
        if(direction.1[interim.index]=='Greater'&direction.2[interim.index]=='Greater'){
          go_matrix[,interim.index,index.setting]<-(sumx>=true_go_cutoff1)&(sumy>=true_go_cutoff2)
          table_true_go_cutoff[interim.index,index.setting]<-HTML('X>=',true_go_cutoff1,logic.go.points[interim.index],' Y>=',true_go_cutoff2)
          table_go_cutoff[interim.index,index.setting]<-HTML('X>=',temp1$cutoff[1],logic.go.points[interim.index],' Y>=',temp2$cutoff[1])
        }
        if(direction.1[interim.index]=='Greater'&direction.2[interim.index]=='Less'){
          go_matrix[,interim.index,index.setting]<-(sumx>=true_go_cutoff1)&(sumy<=true_go_cutoff2)
          table_true_go_cutoff[interim.index,index.setting]<-HTML('X>=',true_go_cutoff1,logic.go.points[interim.index],' Y<=',true_go_cutoff2)
          table_go_cutoff[interim.index,index.setting]<-HTML('X>=',temp1$cutoff[1],logic.go.points[interim.index],' Y<=',temp2$cutoff[1])
        }
        if(direction.1[interim.index]=='Less'&direction.2[interim.index]=='Greater'){
          go_matrix[,interim.index,index.setting]<-(sumx<=true_go_cutoff1)&(sumy>=true_go_cutoff2)
          table_true_go_cutoff[interim.index,index.setting]<-HTML('X<=',true_go_cutoff1,logic.go.points[interim.index],' Y>=',true_go_cutoff2)
          table_go_cutoff[interim.index,index.setting]<-HTML('X<=',temp1$cutoff[1],logic.go.points[interim.index],' Y>=',temp2$cutoff[1])
        }
        if(direction.1[interim.index]=='Less'&direction.2[interim.index]=='Less'){
          go_matrix[,interim.index,index.setting]<-(sumx<=true_go_cutoff1)&(sumy<=true_go_cutoff2)
          table_true_go_cutoff[interim.index,index.setting]<-HTML('X<=',true_go_cutoff1,logic.go.points[interim.index],' Y<=',true_go_cutoff2)
          table_go_cutoff[interim.index,index.setting]<-HTML('X<=',temp1$cutoff[1],logic.go.points[interim.index],' Y<=',temp2$cutoff[1])

        }
      }
      if(logic.go.points[interim.index]=='or'){
        if(direction.1[interim.index]=='Greater'&direction.2[interim.index]=='Greater'){
          go_matrix[,interim.index,index.setting]<-(sumx>=true_go_cutoff1)|(sumy>=true_go_cutoff2)
          table_true_go_cutoff[interim.index,index.setting]<-HTML('X>=',true_go_cutoff1,logic.go.points[interim.index],' Y>=',true_go_cutoff2)
          table_go_cutoff[interim.index,index.setting]<-HTML('X>=',temp1$cutoff[1],logic.go.points[interim.index],' Y>=',temp2$cutoff[1])
        }
        if(direction.1[interim.index]=='Greater'&direction.2[interim.index]=='Less'){
          go_matrix[,interim.index,index.setting]<-(sumx>=true_go_cutoff1)|(sumy<=true_go_cutoff2)
          table_true_go_cutoff[interim.index,index.setting]<-HTML('X>=',true_go_cutoff1,logic.go.points[interim.index],' Y<=',true_go_cutoff2)
          table_go_cutoff[interim.index,index.setting]<-HTML('X>=',temp1$cutoff[1],logic.go.points[interim.index],' Y<=',temp2$cutoff[1])
        }
        if(direction.1[interim.index]=='Less'&direction.2[interim.index]=='Greater'){
          go_matrix[,interim.index,index.setting]<-(sumx<=true_go_cutoff1)|(sumy>=true_go_cutoff2)
          table_true_go_cutoff[interim.index,index.setting]<-HTML('X<=',true_go_cutoff1,logic.go.points[interim.index],' Y>=',true_go_cutoff2)
          table_go_cutoff[interim.index,index.setting]<-HTML('X<=',temp1$cutoff[1],logic.go.points[interim.index],' Y>=',temp2$cutoff[1])
        }
        if(direction.1[interim.index]=='Less'&direction.2[interim.index]=='Less'){
          go_matrix[,interim.index,index.setting]<-(sumx<=true_go_cutoff1)|(sumy<=true_go_cutoff2)
          table_true_go_cutoff[interim.index,index.setting]<-HTML('X<=',true_go_cutoff1,logic.go.points[interim.index],' Y<=',true_go_cutoff2)
          table_go_cutoff[interim.index,index.setting]<-HTML('X<=',temp1$cutoff[1],logic.go.points[interim.index],' Y<=',temp2$cutoff[1])

        }
      }

      if(logic.nogo.points[interim.index]=='and'){
        if(direction.1[interim.index]=='Greater'&direction.2[interim.index]=='Greater'){
          nogo_matrix[,interim.index,index.setting]<-(sumx<true_nogo_cutoff1)&(sumy<true_nogo_cutoff2)
          table_true_nogo_cutoff[interim.index,index.setting]<-HTML('X<',true_nogo_cutoff1,logic.nogo.points[interim.index],' Y<',true_nogo_cutoff2)
          table_nogo_cutoff[interim.index,index.setting]<-HTML('X<',temp1$cutoff[2],logic.nogo.points[interim.index],' Y<',temp2$cutoff[2])

        }
        if(direction.1[interim.index]=='Greater'&direction.2[interim.index]=='Less'){
          nogo_matrix[,interim.index,index.setting]<-(sumx<true_nogo_cutoff1)&(sumy>true_nogo_cutoff2)
          table_true_nogo_cutoff[interim.index,index.setting]<-HTML('X<',true_nogo_cutoff1,logic.nogo.points[interim.index],' Y>',true_nogo_cutoff2)
          table_nogo_cutoff[interim.index,index.setting]<-HTML('X<',temp1$cutoff[2],logic.nogo.points[interim.index],' Y>',temp2$cutoff[2])

        }
        if(direction.1[interim.index]=='Less'&direction.2[interim.index]=='Greater'){
          nogo_matrix[,interim.index,index.setting]<-(sumx>true_nogo_cutoff1)&(sumy<true_nogo_cutoff2)
          table_true_nogo_cutoff[interim.index,index.setting]<-HTML('X>',true_nogo_cutoff1,logic.nogo.points[interim.index],' Y<',true_nogo_cutoff2)
          table_nogo_cutoff[interim.index,index.setting]<-HTML('X>',temp1$cutoff[2],logic.nogo.points[interim.index],' Y<',temp2$cutoff[2])

        }
        if(direction.1[interim.index]=='Less'&direction.2[interim.index]=='Less'){
          nogo_matrix[,interim.index,index.setting]<-(sumx>true_nogo_cutoff1)&(sumy>true_nogo_cutoff2)
          table_true_nogo_cutoff[interim.index,index.setting]<-HTML('X>',true_nogo_cutoff1,logic.nogo.points[interim.index],' Y>',true_nogo_cutoff2)
          table_nogo_cutoff[interim.index,index.setting]<-HTML('X>',temp1$cutoff[2],logic.nogo.points[interim.index],' Y>',temp2$cutoff[2])

        }
      }

      if(logic.nogo.points[interim.index]=='or'){
        if(direction.1[interim.index]=='Greater'&direction.2[interim.index]=='Greater'){
          nogo_matrix[,interim.index,index.setting]<-(sumx<true_nogo_cutoff1)|(sumy<true_nogo_cutoff2)
          table_true_nogo_cutoff[interim.index,index.setting]<-HTML('X<',true_nogo_cutoff1,logic.nogo.points[interim.index],' Y<',true_nogo_cutoff2)
          table_nogo_cutoff[interim.index,index.setting]<-HTML('X<',temp1$cutoff[2],logic.nogo.points[interim.index],' Y<',temp2$cutoff[2])

        }
        if(direction.1[interim.index]=='Greater'&direction.2[interim.index]=='Less'){
          nogo_matrix[,interim.index,index.setting]<-(sumx<true_nogo_cutoff1)|(sumy>true_nogo_cutoff2)
          table_true_nogo_cutoff[interim.index,index.setting]<-HTML('X<',true_nogo_cutoff1,logic.nogo.points[interim.index],' Y>',true_nogo_cutoff2)
          table_nogo_cutoff[interim.index,index.setting]<-HTML('X<',temp1$cutoff[2],logic.nogo.points[interim.index],' Y>',temp2$cutoff[2])

        }
        if(direction.1[interim.index]=='Less'&direction.2[interim.index]=='Greater'){
          nogo_matrix[,interim.index,index.setting]<-(sumx>true_nogo_cutoff1)|(sumy<true_nogo_cutoff2)
          table_true_nogo_cutoff[interim.index,index.setting]<-HTML('X>',true_nogo_cutoff1,logic.nogo.points[interim.index],' Y<',true_nogo_cutoff2)
          table_nogo_cutoff[interim.index,index.setting]<-HTML('X>',temp1$cutoff[2],logic.nogo.points[interim.index],' Y<',temp2$cutoff[2])

        }
        if(direction.1[interim.index]=='Less'&direction.2[interim.index]=='Less'){
          nogo_matrix[,interim.index,index.setting]<-(sumx>true_nogo_cutoff1)|(sumy>true_nogo_cutoff2)
          table_true_nogo_cutoff[interim.index,index.setting]<-HTML('X>',true_nogo_cutoff1,logic.nogo.points[interim.index],' Y>',true_nogo_cutoff2)
          table_nogo_cutoff[interim.index,index.setting]<-HTML('X>',temp1$cutoff[2],logic.nogo.points[interim.index],' Y>',temp2$cutoff[2])

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


      table_RR[interim.index,index.setting]<-paste0(RR[index.setting,],collapse = ',')
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
     table[2,j]=paste0(round(RR[meanindex,],3),collapse = ',')
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
   table[2,num_interim+1]=paste0(round(RR[meanindex,],3),collapse = ',')
   table[3,num_interim+1]=''
   table[4,num_interim+1]=round(sum(as.numeric(table[4,1:num_interim])),3)
   table[5,num_interim+1]=round(as.numeric(table[5,num_interim]),3)
   table[6,num_interim+1]=round(sum(as.numeric(table[6,1:num_interim])),3)
   table[7,num_interim+1]=''
   table[8,num_interim+1]=''
   table[9,num_interim+1]=''
   table[10,num_interim+1]=''
   tablecolname<-c(paste0('Interim analysis ',1:(num_interim-1)),'Final analysis',"Summary")
   tablerowname<-c('Sample size','P<sub>00</sub>,P<sub>01</sub>,P<sub>10</sub>,P<sub>11</sub>','Task','Success','To next interim/final or inconclusive',
                   'Stop','Superority/Futility zone','Warning',
                   'Cut off for NOGO rule',
                   'Cut off for GO rule')
   table<-cbind(tablerowname,rep(meanindex,10),table)
   colnames(table)<-c(" ",'Setting',tablecolname)
   temptable=rbind(temptable,table)
 }
  return(temptable)
}


#CB_OC_table()
#CB_OC_table_interim()

