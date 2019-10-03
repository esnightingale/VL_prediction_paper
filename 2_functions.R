# Author:         Emily S Nightingale
# Institutions:   London Schoool of Hygiene and Tropical Medicine, London, UK
# Date Published: October 2019
################################################################################
# User-defined functions required for model assessment, plotting and 
# forecasting
################################################################################
################################################################################


################################################################################
# For random model draw
################################################################################

random.model<-function(sts.object, n, tp=c(48,71)){
  
  seas1<-addSeason2formula(~1,S=1, period=sts.object@freq)
  seas1t<-addSeason2formula(~1+t,S=1, period=sts.object@freq)
  seas1cov<-addSeason2formula(~1+logpopdens,S=1, period=sts.object@freq)
  seas1tcov<-addSeason2formula(~1+t+logpopdens,S=1, period=sts.object@freq)
  END.form.list<-list(~1,~1+t,~1+t+logpopdens,~1+logpopdens,seas1,seas1t,seas1cov,seas1tcov)
  END.offset.list<-list(NULL,population(sts.object))
  AR.form.list<-list(~1,~1+t,~1+t+logpopdens,~1+logpopdens,seas1,seas1t,seas1cov,seas1tcov)
  AR.lag.list<-c(1:12)
  NE.form.list<-list(~1,~1+t,~1+t+logpopdens,~1+logpopdens,seas1,seas1t,seas1cov,seas1tcov)
  NE.lag.list<-c(1:7)  #using powerlaw syntax with maxlag=1 same as neighbourhood(sts.object)==1
  fam.list <- list("Poisson","NegBin1",as.factor(df_wide$State),as.factor(df_wide$District))
  
  n.opt<-sapply(list(END.form.list,END.offset.list,AR.form.list,AR.lag.list,NE.form.list,NE.lag.list,fam.list),length)
  
  models <- list()
  osa.first <- list()
  osa.roll <- list()
  
  while (length(models) < n){
    
    rand.opt <- sapply(n.opt,FUN=sample,size=1)
    
    END.form <- END.form.list[[rand.opt[1]]]
    END.offset <- END.offset.list[[rand.opt[2]]]
    AR.form <- AR.form.list[[rand.opt[3]]]
    AR.lag <- AR.lag.list[[rand.opt[4]]]
    NE.form <- NE.form.list[[rand.opt[5]]]
    NE.lag <- NE.lag.list[[rand.opt[6]]]
    fam <- fam.list[[rand.opt[7]]]
    
    
    # print(paste0("END = ",END.form,", offset = ",END.offset,", AR = ",AR.form,", lag = ",AR.lag,", NE = ",NE.form,", sp.lag = ",NE.lag,", family = ",fam))
    
    control <- list(end = list(f = END.form, offset=END.offset), 
                    ar = list(f = AR.form),
                    ne = list(f = NE.form, weights=W_powerlaw(maxlag = NE.lag)), #something wrong with this?
                    subset=subset,
                    max_lag=AR.lag,
                    family = fam)
    
    # Fit models with error catching since some combinations lead to convergence issues
    model <- tryCatch(profile_par_lag(sts.object, control=control), error = function(e) e, warning = function(w) w) #, finally = function(){ models <- list.append(models, model) })
    osa1<-tryCatch(oneStepAhead_hhh4lag(model,tp=tp,type="first",which.start="current",keep.estimates=T), error = function(e) e, warning = function(w) w)
    osa2<-tryCatch(oneStepAhead_hhh4lag(model,tp=tp,type="rolling",which.start="current",keep.estimates=T), error = function(e) e, warning = function(w) w)
    # if(inherits(model, "error")) next 
    # if(inherits(model, "warning")) next   
    if(!(inherits(model, "error")|inherits(model, "warning")|inherits(osa1, "error")|inherits(osa1, "warning")|inherits(osa2, "error")|inherits(osa2, "warning"))){
      print(summary(model))
      models <- list.append(models, model) 
      osa.first<-list.append(osa.first,osa1)
      osa.roll<-list.append(osa.roll,osa2)
    }
    
  }
  return(list(models,osa.first,osa.roll))
}


################################################################################
# For 3_selection.R
################################################################################

# Comparison of spatial lags (determine maximum neighbour order at which block 
# time series influence each other; fit according to a powerlaw)
comp.nblag<-function(model,tp,type,SCORES=c("logs", "rps", "dss", "ses")){
  
  AIC_nblag<-vector(length=7)
  AIC_nblag[1]<-AIC(model)
  scores_nblag<-matrix(nrow=7,ncol=4)
  calibp_nblag<-vector(length=7)
  calibz_nblag<-vector(length=7)
  
  ctl<-model$control
  ctl$ne$weights<-(neighbourhood(stsobj)==1)
  model<-profile_par_lag(model$stsObj,ctl)
  
  osa.list<-list(oneStepAhead_hhh4lag(model, tp=tp, type = type, which.start = "current"))
  scores_nblag[1,]<-colMeans(scores(osa.all[[1]],which=SCORES))
  calibp_nblag[1]<-calibrationTest(osa.all[[1]], which="rps")[["p.value"]]
  calibz_nblag[1]<-calibrationTest(osa.all[[1]], which="rps")[["statistic"]]
  
  for (i in 2:7){
    
    ctl2<-model$control
    ctl2$ne$weights=W_powerlaw(maxlag=i)
    mod<-profile_par_lag(model$stsObj,ctl2)
    AIC_nblag[i]<-AIC(mod)
    osa<-oneStepAhead_hhh4lag(mod, tp=tp, type = type, which.start = "current")
    scores_nblag[i,]<-colMeans(scores(osa,which=SCORES))
    calib<-calibrationTest(osa, which="rps", individual=T)
    calibp_nblag[i]<-calib[["p.value"]]
    calibz_nblag[i]<-calib[["statistic"]]
    osa.list<-list.append(osa.list,osa)
  }
  
  results<-list(AIC_nblag,scores_nblag,calibp_nblag,calibz_nblag)  
  return(list(osa.all,results))
}

# Comparison of temporal lags (determine maximum lag in time at which to fit 
# distributed lags)
comp.arlag<-function(model, tp, type, SCORES=c("logs", "rps", "dss", "ses")){

  subset<-13:(tp[2]+1)
  
  AIC_arlag<-vector(length=12)
  scores_arlag<-matrix(nrow=12,ncol=4)
  calibp_arlag<-vector(length=12)
  calibz_arlag<-vector(length=12)

  mod<-update(model,subset=subset)
  osa<-oneStepAhead_hhh4lag(mod, tp=tp, type = type, which.start = "current")
  osa.list<-list(osa)
  
  AIC_arlag[1]<-AIC(mod)
  scores_arlag[1,]<-colMeans(scores(osa,which=SCORES))
  calibp_arlag[1]<-calibrationTest(osa, which="rps")[["p.value"]]
  calibz_arlag[1]<-calibrationTest(osa, which="rps")[["statistic"]]
  
  for (i in 2:12){
    ctl<-mod$control
    ctl$max_lag = i
    
    mod<-profile_par_lag(stsobj,control=ctl)
    osa<-oneStepAhead_hhh4lag(mod, tp=tp, type = type, which.start = "current")
    osa.list<-list.append(osa.list,osa)
    AIC_arlag[i]<-AIC(mod)
    scores_arlag[i,]<-colMeans(scores(osa,which=SCORES))
    calibp_arlag[i]<-calibrationTest(osa, which="rps")[["p.value"]]
    calibz_arlag[i]<-calibrationTest(osa, which="rps")[["statistic"]]
  }
  
  
  result<-list(AIC_arlag,scores_arlag,calibp_arlag,calibz_arlag)
  return(result)
}

#------------------------------------------------------------------------------#
# Plotting AIC and scoring rules for inspection of optimal lags 
plot.ar <- function(result,file){
  SCORES = c("logs", "rps", "dss", "ses")
  png(filename = paste0(file,"_1.png"), height=330, width=1000)
  par(mfrow=c(1,3))
  plot(c(1:12),result[[1]],type="l", col="red", ylab=" ", main="AIC", xlab="max temporal lag")
  plot(c(1:12),result[[4]],type="l", col="forestgreen", ylab=" ", main="Calibration test - Z", xlab="max temporal lag")
  plot(c(1:12),result[[3]],type="l", col="orange", ylab=" ", main="Calibration test - p-value", xlab="max temporal lag")
  dev.off()
  png(filename = paste0(file,"_2.png"), height=500, width=600)
  par(mfrow=c(2,2))
  for (i in 1:4){
    plot(c(1:12),result[[2]][,i], type="l", main=SCORES[i], xlab="max temporal lag", ylab="")
  }
  dev.off()
}

plot.nb <- function(result,file){
  SCORES = c("logs", "rps", "dss", "ses")
  png(filename = paste0(file,"_1.png"), height=330, width=1000)
  par(mfrow=c(1,3))
  plot(c(1:7),result[[2]][[1]],type="l", col="red", ylab=" ", main="AIC", xlab="max spatial lag")
  plot(c(1:7),result[[2]][[4]],type="l", col="forestgreen", ylab=" ", main="Calibration test - Z", xlab="max spatial lag")
  plot(c(1:7),result[[2]][[3]],type="l", col="orange", ylab=" ", main="Calibration test - p-value", xlab="max spatial lag")
  dev.off()
  png(filename = paste0(file,"_2.png"), height=500, width=600)
  par(mfrow=c(2,2))
  for (i in 1:4){
    plot(c(1:7),result[[2]][[2]][,i], type="l", main=SCORES[i], xlab="max spatial lag", ylab="")
  }
  dev.off()
}

#------------------------------------------------------------------------------#
# Fit all models in a list and calculate fit/prediction metrics using a 
# one-step-ahead prediction approach ("type" defines either fixed or rolling 
# updated parameter estimates).
modelassess <- function(model.list, time.period,type,SCORES=c("logs", "rps", "dss", "ses")){
  
  model.calibr <- vector(length=2)
  
  AIC <- sapply(model.list, AIC)
  model.preds <- lapply(model.list, oneStepAhead_hhh4lag, tp = time.period, 
                        type = type, which.start = "current", 
                        keep.estimates=T)
  all.scores <- lapply(model.preds, scores, which = SCORES, individual = T)
  all.calibr <- lapply(model.preds, calibrationTest, which = "rps", individual = T)
  model.scores <- t(sapply(all.scores, colMeans, dims=2))
  for (m in 1:length(model.list)){
    model.calibr <- rbind(model.calibr,
                          c(all.calibr[[m]]$p.value, all.calibr[[m]]$statistic))}
  
  model.calibr <- model.calibr[-1,]
  return(list(AIC=AIC,pred=model.preds,scores=model.scores, calib=model.calibr))
}

# Calculate overall average scores for each model, across all blocks and predicted
# time points
# score_tidy<-function(x){
#   nmod<-length(x[[3]])
#   avg.scores<-matrix(nrow=nmod,ncol=4)
#   for(i in 1:nmod){
#     avg.scores[i,]<-apply(x[[3]][[i]],MARGIN = 3, FUN="mean", na.rm=T)
#   }
#   colnames(avg.scores)<-c("logs","RPS","DSS","SES")
#   #print(avg.scores)
#   return(avg.scores)}

# Permuatation test of scores between models
permut.test<-function(null.model,new.model, tp, type){
  
  SCORES <- c("logs", "rps", "dss", "ses")
  
  models2compare <- c('null.model', 'new.model')
  models.pred <- lapply(mget(models2compare), oneStepAhead_hhh4lag, tp = tp, type = type, which.start = "current", keep.estimates=T)
  models.scores <- lapply(models.pred, scores, which = SCORES, individual = T)
permut.test <- sapply(SCORES, function (score) permutationTest(
  models.scores$null.model[, , score],
  models.scores$new.model[, , score],
  nPermutation = 10000))
print(permut.test)
# permut.test.rps <- as.numeric(permut.test[,2])
# return(permut.test.rps)
}


################################################################################
# For 4_evaluation.R
# Forecasting n-months-ahead with rolling window and evaluating predictive power
################################################################################

# Not set up to work for models with spatially varying dispersion parameter
stepaheadN<-function(model,
                     start,
                     type="first",
                     n=3){
  
  obs<-model$stsObj@observed
  nblock<-dim(obs)[2]
  t.max<-dim(obs)[1]
  t.min<-model$control$subset[1]
  preds<-matrix(0,nrow=1, ncol=nblock)
  disp<-vector(length=1)  #matrix(0,nrow=1, ncol=nblock)
  prob<-matrix(0,nrow=1, ncol=nblock)
  cprob<-matrix(0,nrow=1, ncol=nblock)
  
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
    var<-t(as.matrix(predmom[["var_matrix"]][n,])) #future variances
    #psi<-as.matrix(((var/mu)-1)/mu)
    row.names(mu)<-t.pred
    #row.names(psi)<-t.pred
    psi<-exp(-mod$coefficients[grepl("disp",names(mod$coefficients))])     # pull dispersion out of model fit
    
    preds<-rbind(preds,mu)
    disp<-c(disp,psi)
    
    if (mod$control$family=="Poisson"){     # Calculate probabilities of observed values according to the predicted distribution
      p<-dpois(obs[t.pred,],lambda=mu)
      c<-ppois(obs[t.pred,],lambda=mu)
      # }else if (length(psi)>1){
      #  psi2<-as.numeric(mod$control$family)
      #  psi2[psi2=="JHARKHAND"]<-psi[1]
      #   psi2[psi2=="BIHAR"]<-psi[2]
      #   p<-dnbinom(obs[t.pred,],size=1/psi2,mu=mu)
      #   c<-pnbinom(obs[t.pred,],size=1/psi2,mu=mu)
    }else{
      p<-dnbinom(obs[t.pred,],size=1/psi,mu=mu)
      c<-pnbinom(obs[t.pred,],size=1/psi,mu=mu)
    }
    
    
    prob<-rbind(prob,p)
    cprob<-rbind(cprob,c)
    
    t<-t+1
  }
  
  preds<-preds[-1,]   #remove redundant rows
  disp<-disp[-1]
  prob<-prob[-1,]
  cprob<-cprob[-1,]
  
  return(list(preds,disp,prob,cprob))
}

#------------------------------------------------------------------------------#
# Calculate n-step-ahead predicted quantile intervals
# Output from stepaheadN: list(preds,disp,prob,cprob)
# Need pred mean and disp to specify the distribution and pull out quantile values

predquants<-function(pred.output,probs=list(0.1,0.25,0.45,0.55,0.75,0.9)){
  
  if(class(pred.output) == "list"){
    mu<-pred.output[[1]]
    disp<-pred.output[[2]]
  }else if(class(pred.output) == "oneStepAhead"){
    mu<-pred.output$pred
    disp<-exp(-pred.output$psi)
  }
  
  qs<-lapply(probs,FUN=qnbinom,size=1/disp,mu=mu)
  
  return(qs)
}

# predquantsOSA<-function(osa.obj,probs=list(0.1,0.25,0.45,0.55,0.75,0.9)){
#   
#   mu<-osa.obj$pred
#   disp<-exp(-osa.obj$psi)
#   
#   if(model$control$family=="Poisson"){qs<-lapply(probs,FUN=qpois,lambda=mu)
#   }else{qs<-lapply(probs,FUN=qnbinom,size=1/disp,mu=mu)}
#   
#   return(qs)
# }


#------------------------------------------------------------------------------#
# Calculate empirical coverage probabilities over three quantile intervals output by 
# predquants
covprob<-function(quants,obs){
  
  score_array<-array(dim=c(dim(obs),3))
  score_qwd_array<-array(dim=c(dim(obs),3))
  score<-vector(length=3)
  score_qwd<-vector(length=3)
  
  score_array[,,1]<-(quants[[1]]<=obs & obs<=quants[[6]])
  score_array[,,2]<-(quants[[2]]<=obs & obs<=quants[[5]])
  score_array[,,3]<-(quants[[3]]<=obs & obs<=quants[[4]])
  
  score_qwd_array[,,1]<-quants[[6]]-quants[[1]]
  score_qwd_array[,,2]<-quants[[5]]-quants[[2]]
  score_qwd_array[,,3]<-quants[[4]]-quants[[3]]
  
  for (i in 1:3){score_array[,,i]<-abs(apply(score_array[,,i],2,FUN=as.numeric)-1)
  score[i]<-mean(score_array[,,i])
  score_qwd[i]<-mean(score_qwd_array[,,i])
  }
  
  return(list(score,score_qwd,score_array,score_qwd_array))
  
} 


################################################################################
# For 6_figures.R 
# Plot continuous measure on region map 

mapplot<-function(shapefile,data,value,legend_title,pal = "A"){
  
  theme_map <- function (base_size = 12, base_family = "") {
    theme_gray(base_size = base_size, base_family = base_family) %+replace% 
      theme(
        axis.line=element_blank(),
        axis.text.x=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks=element_blank(),
        axis.ticks.length=unit(0.3, "lines"),
        axis.ticks.margin=unit(0.5, "lines"),
        axis.title.x=element_blank(),
        axis.title.y=element_blank(),
        legend.background=element_rect(fill="white", colour=NA),
        legend.key=element_rect(colour="white"),
        legend.key.size=unit(1.2, "lines"),
        legend.position="right",
        legend.text=element_text(size=rel(0.8)),
        legend.title=element_text(size=rel(0.8), face="bold", hjust=0),
        panel.background=element_blank(),
        panel.border=element_blank(),
        panel.grid.major=element_blank(),
        panel.grid.minor=element_blank(),
        panel.spacing=unit(0, "lines"),
        plot.background=element_blank(),
        plot.margin=unit(c(1, 1, 0.5, 0.5), "lines"),
        plot.title=element_text(size=rel(1.2)),
        strip.background=element_rect(fill="grey90", colour="grey50"),
        strip.text.x=element_text(size=rel(0.8)),
        strip.text.y=element_text(size=rel(0.8), angle=-90) 
      )   
  }
  
  shp.df <- fortify(shapefile) # transform shapefile into dataframe
  shp.df <- merge(data,shp.df,by="id") # merge the spatial data and the value to be plotted
  shp.df$value <- eval(substitute(value),shp.df)
  
  plot.map <- ggplot(data = shp.df, aes(long, lat, group = group)) 
  plot.map <- plot.map + geom_polygon(aes(fill = value))
  plot.map <- plot.map + geom_path(colour = 'black')
  plot.map <- plot.map + scale_fill_viridis_c(option = pal, name = legend_title) 
  plot.map <- plot.map + coord_equal() + theme_map() #+ geom_text(data=shp.df, aes(x=long, y=lat, label=Block, group=Block), size=0.5)
  return(plot.map)

}

# Set up data frame and plot predicted quantiles
make.df <- function(preds,quants,blk,shift){
  df.pred <- data.frame(Month=c(1:obsmth),Obs=cases[1:obsmth,blk], 
                        Pred=c(rep(NA,shift),preds[,blk]),
                        q10=c(rep(NA,shift),quants[[1]][1,,blk]),
                        q90=c(rep(NA,shift),quants[[2]][1,,blk]),
                        q25=c(rep(NA,shift),quants[[1]][2,,blk]),
                        q75=c(rep(NA,shift),quants[[2]][2,,blk]),
                        q45=c(rep(NA,shift),quants[[1]][3,,blk]),
                        q55=c(rep(NA,shift),quants[[2]][3,,blk]))
  df.pred$correct1 <- "Y"
  df.pred$correct1[df.pred$Obs>df.pred$q90 | df.pred$Obs<df.pred$q10] <- "N"
  df.pred$correct2 <- "Y"
  df.pred$correct2[df.pred$Obs>df.pred$q75 | df.pred$Obs<df.pred$q25] <- "N"
  df.pred$correct3 <- "Y"
  df.pred$correct3[df.pred$Obs>df.pred$q55 | df.pred$Obs<df.pred$q45] <- "N"
  
  return(df.pred)
}

quantplot <- function(df,var,title,legend=F){
  df2 <- df[1:9]
  df2$correct <- df[,names(df)==var]
  df3<-melt(df2[,c(1:4,6,8,10)],id.vars = c("Month","Obs","Pred","correct"))
  names(df3)[5:6] <- c("quant1", "min")
  df3[7:8] <- melt(df2[,c(1:3,5,7,9)],id.vars = c("Month","Obs","Pred"))[,4:5]
  names(df3)[7:8] <- c("quant2", "max")
  
  p <- ggplot(df3, aes(Month,Obs), ylab = "No. reported cases") +
    geom_ribbon(aes(ymin=min, ymax=max, fill=quant1), show.legend = T, 
                na.rm = T, alpha=0.5) +
    geom_line(aes(Month, Pred), na.rm = T, col = "white",lwd = 1.2) +
    geom_point(aes(shape = correct)) +
    scale_shape_manual(values = c(4, 19)) +
    scale_fill_manual(values=c("gold","orange","red"), 
                      name = "Predicted quantiles",
                      labels=c("10 - 90%", "25 - 75%", "45 - 55%")) +
    theme(legend.title = element_text(size = 10), 
          legend.text = element_text(size = 10),
          legend.justification=c(0,1), 
          legend.position=c(0.1,1),
          text = element_text(size = 14),
          axis.text = element_text(size = 14)) +
    ylab("No. reported cases") +
    guides(shape = FALSE) +
    ggtitle(title)  
  
  if(legend == T){return(p)}
  else{return(p + theme(legend.position="none"))}
  
}