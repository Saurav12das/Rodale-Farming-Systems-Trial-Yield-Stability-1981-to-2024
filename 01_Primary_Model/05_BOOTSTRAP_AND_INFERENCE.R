# Parametric bootstrap for mean yield, annual variability and system contrasts.
# Simulated datasets retain the fitted covariance and time-dependence structures.
# Simultaneous intervals account for three planned contrasts within each crop.
v9_md5 <- function() unname(tools::md5sum(c(v9_file(YIELD_FILE),v9_file(WEATHER_FILE))))
v9_init_checkpoint <- function(path,required,md5) {
  if(!file.exists(path)) return(data.frame())
  z<-tryCatch(read.csv(path,stringsAsFactors=FALSE),error=function(e)data.frame())
  if(!nrow(z)) return(data.frame())
  if(!all(required %in% names(z))) stop('Old checkpoint incompatible with current code: ',path)
  if(any(z$Input_MD5_1!=md5[1]|z$Input_MD5_2!=md5[2]))
    stop('Input files changed since saved bootstrap checkpoint: ',path)
  z
}
v9_boot_main <- function(pd,fit,nboot,tag,seedoffset=0L) {
  dest<-v9_output(paste0('S5_BOOT_draws_',tag,'.csv'),'bootstrap')
  statfile<-v9_output(paste0('S5_BOOT_status_',tag,'.csv'),'bootstrap')
  md5<-v9_md5();kind<-fit$kind;serial<-fit$serial
  required<-c('Replicate','Input_MD5_1','Input_MD5_2','System')
  draws<-v9_init_checkpoint(dest,required,md5)
  status<-v9_init_checkpoint(statfile,c('Replicate','Input_MD5_1','Input_MD5_2','Success'),md5)
  if(nrow(draws)&&any(draws$Model_tag!=tag)) stop('Checkpoint tag changed')
  done<-if(nrow(status)) unique(status$Replicate) else integer()
  if(any(done>nboot)) message('Checkpoint contains more replicates than requested; summary uses first ',nboot)
  for(b in seq_len(nboot)) {
    if(b%in%done) next
    set.seed(as.integer(SEED+seedoffset+b))
    sim<-pd; sim$y<-v9_simulate_y(pd,fit)
    f<-v9_fit(sim,kind,serial,'REML',warm=fit$par,fast=TRUE,multi=FALSE)
    if(!isTRUE(f$ok)||f$convergence!=0)
      f<-v9_fit(sim,kind,serial,'REML',warm=fit$par,fast=FALSE,multi=TRUE)
    okay<-isTRUE(f$ok)&&f$convergence==0
    if(okay){
      q<-v9_standardized(sim,f)
      q$Replicate<-b;q$Model_tag<-tag
      q$Input_MD5_1<-md5[1];q$Input_MD5_2<-md5[2]
      draws<-rbind(draws,q)
    }
    status<-rbind(status,data.frame(Replicate=b,Crop=pd$crop,Model_tag=tag,
      Input_MD5_1=md5[1],Input_MD5_2=md5[2],
      Success=okay,Boundary=if(isTRUE(f$ok)) v9_boundary(f) else 'FAILED',
      Optimizer=if(isTRUE(f$ok)) f$convergence else NA_integer_))
    if(b%%CHECKPOINT_EVERY==0L||b==nboot){
      v9_save(draws,dest);v9_save(status,statfile)
      cat('  BOOT ',tag,': ',b,'/',nboot,' (successful ',sum(status$Success),')\n',sep='')
    }
  }
  draws<-draws[draws$Replicate %in% seq_len(nboot),,drop=FALSE]
  status<-status[status$Replicate %in% seq_len(nboot),,drop=FALSE]
  nsuccess<-sum(status$Success)
  if(nsuccess<.95*nboot) warning(tag,' bootstrap success rate below 95%.')
  list(draws=draws,status=status,nsuccess=nsuccess,nrequested=nboot)
}
v9_boot_summary <- function(point,draws,nboot,tag) {
  out<-list()
  metrics<-c('Mean_Mg_ha','Latent_residual_SD_Mg_ha','Detrended_total_SD_Mg_ha',
    'Detrended_total_CV_percent','Trend_including_total_SD_Mg_ha',
    'Trend_including_total_CV_percent')
  for(s in unique(point$System)) for(m in metrics) {
    p<-point[point$System==s,,drop=FALSE]
    d<-draws[draws$System==s,,drop=FALSE]
    z<-d[[m]];z<-z[is.finite(z)]
    ci<-if(length(z)>30) quantile(z,c(.025,.975),names=FALSE) else c(NA_real_,NA_real_)
    out[[length(out)+1L]]<-data.frame(Crop=p$Crop,System=s,Metric=m,Estimate=p[[m]],
      CI_low=ci[1],CI_high=ci[2],N_success=length(z),N_requested=nboot,
      Model_tag=tag,CI_type='pointwise conditional parametric percentile')
  }
  v9_bind(out)
}
v9_pairwise_draws <- function(draws){
  out<-list()
  for(b in sort(unique(draws$Replicate))) {
    z<-draws[draws$Replicate==b,,drop=FALSE]
    if(nrow(z)!=3L) next
    pp<-v9_pairwise(z)
    pp$Replicate<-b
    out[[length(out)+1L]]<-pp
  }
  v9_bind(out)
}
v9_pairwise_intervals <- function(point,draws,nboot,tag) {
  pts<-v9_pairwise(point)
  bd<-v9_pairwise_draws(draws)
  if(!nrow(bd)) return(list(summary=data.frame(),draws=bd))
  rows<-list();metrics<-unique(pts$Metric)
  for(m in metrics){
    p<-pts[pts$Metric==m,,drop=FALSE]
    b<-bd[bd$Metric==m,,drop=FALSE]
    pairs<-p$Comparison
    vals<-sapply(pairs,function(cp){
      zz<-b[b$Comparison==cp,,drop=FALSE]
      zz$Estimate[match(sort(unique(b$Replicate)),zz$Replicate)]
    })
    if(is.null(dim(vals))) vals<-matrix(vals,ncol=length(pairs))
    colnames(vals)<-pairs
    estim<-p$Estimate[match(pairs,p$Comparison)]
    isratio<-all(p$Type=='ratio')
    zmat<-if(isratio) log(vals) else vals
    pest<-if(isratio) log(estim) else estim
    enull<-if(isratio) rep(0,length(estim)) else p$Null
    sdvec<-apply(zmat,2,sd,na.rm=TRUE)
    # Maximum-statistic correction for the three pairwise contrasts.
    if(nrow(zmat)>=30&&all(is.finite(sdvec)&sdvec>1e-12)) {
      zcent<-sweep(zmat,2,pest,'-')
      zcent<-sweep(zcent,2,sdvec,'/')
      max_t<-apply(abs(zcent),1,max)
      max_t<-max_t[is.finite(max_t)]
      crit<-if(length(max_t)>=30) as.numeric(quantile(max_t,FAMILYWISE_CONFIDENCE)) else NA_real_
    } else {crit<-NA_real_;max_t<-numeric()}
    for(j in seq_along(pairs)) {
      x<-zmat[,j];x<-x[is.finite(x)]
      ci<-if(length(x)>=30) as.numeric(quantile(x,c(.025,.975))) else c(NA_real_,NA_real_)
      cisim<-if(is.finite(crit)&&is.finite(sdvec[j])) pest[j]+c(-1,1)*crit*sdvec[j]
        else c(NA_real_,NA_real_)
      pmax_t<-if(length(max_t)>10&&is.finite(sdvec[j])&&sdvec[j]>0)
        (1+sum(max_t>=abs((pest[j]-enull[j])/sdvec[j])))/(1+length(max_t)) else NA_real_
      if(isratio){ci<-exp(ci);cisim<-exp(cisim)}
      rows[[length(rows)+1L]]<-data.frame(Crop=p$Crop[j],Comparison=pairs[j],
        Metric=m,Type=p$Type[j],Null=p$Null[j],Estimate=estim[j],
        Pointwise_CI_low=ci[1],Pointwise_CI_high=ci[2],
        Simultaneous_maxT_CI_low=cisim[1],Simultaneous_maxT_CI_high=cisim[2],
        Bootstrap_maxT_adjusted_p_approx=pmax_t,
        Family='three system contrasts within crop and metric',
        N_boot=length(x),N_requested=nboot,Model_tag=tag)
    }
  }
  list(summary=v9_bind(rows),draws=bd)
}
v9_run_selected_bootstrap <- function(models) {
  if(!RUN_BOOTSTRAP) return(invisible(NULL))
  nb<-if(RUN_MODE=='FULL') MAIN_BOOT_FULL else MAIN_BOOT_SMOKE
  alls<-list();pairs<-list();valid<-list(); sens<-list();
  for(j in seq_along(CROPS)){
    cr<-CROPS[j]; pd<-models$pds[[cr]]; fit<-models$fits[[cr]]
    tag<-paste(cr,pd$trend,fit$kind,fit$serial,sep='_')
    cat('\nBOOTSTRAP main: ',tag,' | ',nb,' draws\n',sep='')
    bs<-v9_boot_main(pd,fit,nb,tag,seedoffset=100000L*j)
    pt<-v9_standardized(pd,fit)
    alls[[cr]]<-v9_boot_summary(pt,bs$draws,nb,tag)
    qq<-v9_pairwise_intervals(pt,bs$draws,nb,tag)
    pairs[[cr]]<-qq$summary
    v9_save(qq$draws,v9_output(paste0('S5_BOOT_pairwise_draws_',cr,'.csv'),'bootstrap'))
    valid[[cr]]<-data.frame(Crop=cr,Model=tag,N_requested=nb,N_success=bs$nsuccess,
      N_failed=nb-bs$nsuccess,Percent_success=100*bs$nsuccess/nb,
      Interpretation=if(RUN_MODE=='SMOKE') 'PRELIMINARY SMOKE TEST' else 'FULL')
  }
  v9_save(v9_bind(alls),v9_output('MAIN_04_bootstrap_mean_SD_CV_95CI.csv'))
  v9_save(v9_bind(pairs),v9_output('MAIN_05_pairwise_differences_ratios_simultaneous_95CI.csv'))
  v9_save(v9_bind(valid),v9_output('S5_04_primary_bootstrap_validation.csv'))
  # Additional corn UN + AR(1) bootstrap estimates and intervals.
  if(RUN_CORN_UN_BOOTSTRAP){
    cr<-'corn';pd<-models$pds[[cr]]
    alt<-v9_fit(pd,'UN','AR1','REML')
    if(isTRUE(alt$ok)&&alt$convergence==0){
      nalt<-if(RUN_MODE=='FULL') SENSITIVITY_BOOT_FULL else SENSITIVITY_BOOT_SMOKE
      tag<-'corn_UN_AR1_key_sensitivity'
      bs<-v9_boot_main(pd,alt,nalt,tag,seedoffset=900000L)
      p<-v9_standardized(pd,alt)
      si<-v9_boot_summary(p,bs$draws,nalt,tag)
      pp<-v9_pairwise_intervals(p,bs$draws,nalt,tag)
      v9_save(si,v9_output('S5_05_corn_UN_AR1_stability_bootstrap_95CI.csv'))
      v9_save(pp$summary,v9_output('S5_06_corn_UN_AR1_pairwise_bootstrap_95CI.csv'))
      v9_save(data.frame(Crop='corn',Covariance='UN',Serial='AR1',N_requested=nalt,
        N_success=bs$nsuccess),v9_output('S5_07_corn_UN_AR1_bootstrap_validation.csv'))
    } else warning('Corn UN AR1 sensitivity fit failed; no alternative bootstrap intervals.')
  }
}
v9_run_omnibus <- function(models) {
  if(!RUN_BOOTSTRAP) return(invisible(NULL))
  B<-if(RUN_MODE=='FULL') OMNIBUS_BOOT_FULL else OMNIBUS_BOOT_SMOKE
  out<-list(); md5<-v9_md5()
  for(j in seq_along(CROPS)){
    cr<-CROPS[j]; pd<-models$pds[[cr]]; sp<-PRIMARY[[cr]]
    full<-models$fits[[cr]]
    null<-v9_fit(pd,'B_EQUAL',sp$serial,'REML')
    if(!isTRUE(null$ok)||null$convergence!=0){
      warning('Omnibus null model failed for ',cr)
      next
    }
    lr_obs<-max(0,2*(full$logLik-null$logLik))
    path<-v9_output(paste0('S5_OMNIBUS_null_simulations_',cr,'.csv'),'bootstrap')
    cur<-v9_init_checkpoint(path,c('Replicate','Input_MD5_1','Input_MD5_2'),md5)
    done<-if(nrow(cur)) unique(cur$Replicate) else integer()
    for(b in seq_len(B)){
      if(b%in%done) next
      set.seed(as.integer(SEED+2000000L+j*100000L+b))
      sim<-pd;sim$y<-v9_simulate_y(pd,null)
      f0<-v9_fit(sim,'B_EQUAL',sp$serial,'REML',warm=null$par,fast=TRUE,multi=FALSE)
      f1<-v9_fit(sim,'B',sp$serial,'REML',warm=full$par,fast=TRUE,multi=FALSE)
      if(!isTRUE(f0$ok)||f0$convergence!=0)
        f0<-v9_fit(sim,'B_EQUAL',sp$serial,'REML',warm=null$par,fast=FALSE,multi=TRUE)
      if(!isTRUE(f1$ok)||f1$convergence!=0)
        f1<-v9_fit(sim,'B',sp$serial,'REML',warm=full$par,fast=FALSE,multi=TRUE)
      valid<-isTRUE(f0$ok)&&isTRUE(f1$ok)&&f0$convergence==0&&f1$convergence==0
      cur<-rbind(cur,data.frame(Replicate=b,Input_MD5_1=md5[1],Input_MD5_2=md5[2],
        Success=valid,LR_stat=if(valid) max(0,2*(f1$logLik-f0$logLik)) else NA_real_))
      if(b%%CHECKPOINT_EVERY==0L||b==B) v9_save(cur,path)
    }
    ok<-cur[cur$Replicate%in%seq_len(B)&cur$Success&is.finite(cur$LR_stat),,drop=FALSE]
    p<-(1+sum(ok$LR_stat>=lr_obs))/(1+nrow(ok))
    out[[cr]]<-data.frame(Crop=cr,Trend=pd$trend,Covariance='B',Serial=sp$serial,
      Null='equal system-specific latent variances; common shared-year variance',
      LR_observed=lr_obs,Approx_chisq_p=pchisq(lr_obs,2,lower.tail=FALSE),
      Parametric_bootstrap_p=p,N_requested=B,N_success=nrow(ok),
      Reliable=if(nrow(ok)>=.95*B) 'YES' else 'CHECK_FAILED_REFITS',
      Note='Test equal latent temporal variances, not equality of mean yields')
  }
  v9_save(v9_bind(out),v9_output('MAIN_06_omnibus_equal_temporal_variance_bootstrap.csv'))
}
