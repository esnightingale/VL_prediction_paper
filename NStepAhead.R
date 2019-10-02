stepaheadN<-function(model,
                     start,
                     type="first",
                     n=3,
                     prob=0.9){
  
  obs<-model$stsObj@observed
  nblock<-dim(obs)[2]
  t.max<-dim(obs)[1]
  t.min<-model$control$subset[1]
  preds<-matrix(0,nrow=1, ncol=nblock)
  disp<-vector(length=1)  #matrix(0,nrow=1, ncol=nblock)
  quants<-matrix(0,nrow=1, ncol=nblock)

  mod<-update(model,subset=t.min:start)     #refit model up to starting month
  
  t<-0
  
  while ((start+t+n)<=t.max){     # while the month to be predicted is within observed period
    
    #print(t)
    
    t.fit<-start+t    # time point up to which model is refitted
    t.pred<-start+t+n     # time point to be predicted
    window<-c(t.fit,t.pred)
    print(window)
    
    if (type=="rolling"){mod<-update(model,subset=t.min:t.fit)}     #refit model up to current month
    
    predmom<-predictive_moments(mod,t_condition=t.fit,lgt=n)    # forecast n months ahead of current fit
    mu<-t(as.matrix(predmom[["mu_matrix"]][n,]))    #predicted means, n x n.block. Select prediction n months ahead
   #var<-t(as.matrix(predmom[["var_matrix"]][n,])) #future variances
    #psi<-as.matrix(((var/mu)-1)/mu)
    row.names(mu)<-t.pred
    #row.names(psi)<-t.pred
    psi<-exp(-mod$coefficients[grepl("disp",names(mod$coefficients))])     # pull dispersion out of model fit
    
    preds<-rbind(preds,mu)
    disp<-c(disp,psi)
    
    quants<-rbind(quants,qnbinom(p=prob, size=1/disp, mu=mu))
    
    t <- t+1}
  
  preds<-preds[-1,]   #remove redundant rows
  disp<-disp[-1]
  quants<-quants[-1,]
  
  return(list(preds,disp,quants))
}
