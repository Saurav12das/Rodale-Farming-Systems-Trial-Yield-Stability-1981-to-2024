# Matched-year full- versus reduced-tillage analyses, 2008–2024.
# Each tillage treatment retains its own annual sampling variance.
v9_pd_tillage <- function(m,cr,ss){
  z<-m[m$Crop==cr&m$System==ss,,drop=FALSE]
  z<-z[order(z$Year),,drop=FALSE]
  a<-data.frame(Year=z$Year,Arm='FT',Yield_Mg_ha=z$Yield_Mg_ha_FT,vi=z$vi_FT)
  b<-data.frame(Year=z$Year,Arm='RT',Yield_Mg_ha=z$Yield_Mg_ha_RT,vi=z$vi_RT)
  d<-rbind(a,b)
  d<-d[order(d$Year,match(d$Arm,c('FT','RT'))),,drop=FALSE]
  d$Arm<-factor(d$Arm,levels=c('FT','RT'))
  X<-model.matrix(~0+Arm,data=d)
  list(crop=cr,system=ss,trend='POST2008_MATCHED',weather='NONE',
    d=d,X=X,formula=~0+Arm,sys_names=c('FT','RT'),ng=2L,
    sys=as.integer(d$Arm),y=d$Yield_Mg_ha,vi=d$vi,n=nrow(d),p=ncol(X),
    gaps=abs(outer(d$Year,d$Year,'-')),
    logdet_XtX=2*sum(log(diag(chol(crossprod(X))))))
}
v9_tillage_point <- function(pd,f,year_pairs) {
  mu<-f$beta
  # Estimate the FT and RT means from the matched-year observations.
  i<-match('ArmFT',names(mu));j<-match('ArmRT',names(mu))
  if(anyNA(c(i,j))) stop('Unexpected tillage coefficient names')
  diff<-mu[j]-mu[i]
  vec<-rep(0,length(mu));vec[j]<-1;vec[i]<- -1
  sed<-sqrt(as.numeric(t(vec)%*%f$beta_cov%*%vec))
  sds<-sqrt(diag(f$S))
  cvs<-100*sds/c(mu[i],mu[j])
  data.frame(Crop=pd$crop,System=pd$system,N_paired_years=year_pairs,
    Earliest_matched=min(pd$d$Year),Latest_matched=max(pd$d$Year),
    FT_mean_Mg_ha=mu[i],RT_mean_Mg_ha=mu[j],
    RT_minus_FT_Mg_ha=diff,RT_minus_FT_SE=sed,
    RT_minus_FT_Wald_CI_low=diff-1.96*sed,
    RT_minus_FT_Wald_CI_high=diff+1.96*sed,
    RT_minus_FT_Wald_p=2*pnorm(-abs(diff/sed)),
    FT_temporal_SD_Mg_ha=sds[1],RT_temporal_SD_Mg_ha=sds[2],
    RT_over_FT_SD_ratio=sds[2]/sds[1],
    FT_temporal_CV_percent=cvs[1],RT_temporal_CV_percent=cvs[2],
    RT_over_FT_CV_ratio=cvs[2]/cvs[1],
    Covariance='B',Serial='NONE',
    Boundary=v9_boundary(f))
}
v9_tillage_boot <- function(pd,fit,point,nboot,tag,seedoffset){
  md5<-v9_md5()
  file<-v9_output(paste0('S6_BOOT_tillage_',tag,'.csv'),'bootstrap')
  statusfile<-v9_output(paste0('S6_BOOT_tillage_status_',tag,'.csv'),'bootstrap')
  draws<-v9_init_checkpoint(file,c('Replicate','Input_MD5_1','Input_MD5_2'),md5)
  statuses<-v9_init_checkpoint(statusfile,c('Replicate','Input_MD5_1','Input_MD5_2'),md5)
  done<-if(nrow(statuses)) unique(statuses$Replicate) else integer()
  for(b in seq_len(nboot)){
    if(b%in%done) next
    set.seed(as.integer(SEED+seedoffset+b))
    sim<-pd;sim$y<-v9_simulate_y(pd,fit)
    f<-v9_fit(sim,'B','NONE','REML',warm=fit$par,fast=TRUE,multi=FALSE)
    if(!isTRUE(f$ok)||f$convergence!=0)
      f<-v9_fit(sim,'B','NONE','REML',warm=fit$par,fast=FALSE,multi=TRUE)
    ok<-isTRUE(f$ok)&&f$convergence==0
    if(ok){
      q<-v9_tillage_point(sim,f,point$N_paired_years)
      q$Replicate<-b;q$Input_MD5_1<-md5[1];q$Input_MD5_2<-md5[2]
      draws<-rbind(draws,q)
    }
    statuses<-rbind(statuses,data.frame(Replicate=b,Input_MD5_1=md5[1],
      Input_MD5_2=md5[2],Success=ok))
    if(b%%CHECKPOINT_EVERY==0L||b==nboot){
      v9_save(draws,file);v9_save(statuses,statusfile)
    }
  }
  draws<-draws[draws$Replicate%in%seq_len(nboot),,drop=FALSE]
  statuses<-statuses[statuses$Replicate%in%seq_len(nboot),,drop=FALSE]
  vars<-c('RT_minus_FT_Mg_ha','RT_over_FT_SD_ratio','RT_over_FT_CV_ratio',
          'FT_temporal_SD_Mg_ha','RT_temporal_SD_Mg_ha',
          'FT_temporal_CV_percent','RT_temporal_CV_percent')
  rows<-lapply(vars,function(v){
    x<-draws[[v]];x<-x[is.finite(x)]
    ci<-if(length(x)>=30) as.numeric(quantile(x,c(.025,.975))) else c(NA_real_,NA_real_)
    data.frame(Crop=pd$crop,System=pd$system,Metric=v,
      Estimate=point[[v]],CI_low=ci[1],CI_high=ci[2],
      N_requested=nboot,N_success=length(x),
      Status=if(RUN_MODE=='SMOKE') 'PRELIMINARY' else 'FULL')
  })
  v9_bind(rows)
}
v9_run_tillage <- function(data){
  if(!RUN_POST2008_TILLAGE) return(invisible(NULL))
  m<-data$matched
  coverage<-v9_bind(lapply(CROPS,function(cr) v9_bind(lapply(SYSTEMS,function(ss){
    z<-m[m$Crop==cr & m$System==ss,,drop=FALSE]
    data.frame(Crop=cr,System=ss,N_matched_calendar_years=nrow(z),
      Years=paste(sort(z$Year),collapse=';'))
  }))))
  v9_save(coverage,v9_output('S6_01_tillage_paired_year_coverage.csv'))
  points<-list();cirows<-list();skip<-list();alts<-list();nr<-0
  nboot<-if(RUN_MODE=='FULL') TILLAGE_BOOT_FULL else TILLAGE_BOOT_SMOKE
  for(cr in CROPS) for(ss in SYSTEMS){
    z<-m[m$Crop==cr&m$System==ss,,drop=FALSE]
    n<-nrow(z)
    if(n>=2){
      diff<-z$Yield_Mg_ha_RT-z$Yield_Mg_ha_FT
      qq<-data.frame(Crop=cr,System=ss,N_years=n,
        Paired_mean_RT_minus_FT=mean(diff),
        Paired_year_SD=sd(diff),Paired_SE=sd(diff)/sqrt(n),
        Paired_t_CI_low=mean(diff)-qt(.975,n-1)*sd(diff)/sqrt(n),
        Paired_t_CI_high=mean(diff)+qt(.975,n-1)*sd(diff)/sqrt(n))
      skip[[length(skip)+1L]]<-qq
    }
    if(n<4){
      cat('TILLAGE ',cr,' ',ss,': only ',n,' matched years. No variance model.\n',sep='')
      next
    }
    pd<-v9_pd_tillage(m,cr,ss)
    f<-v9_fit(pd,'B','NONE','REML')
    if(!isTRUE(f$ok)||f$convergence!=0){
      warning('Tillage B model not converged: ',cr,' ',ss)
      next
    }
    pt<-v9_tillage_point(pd,f,n)
    # Unstructured two-treatment covariance for sensitivity analysis.
    fun<-v9_fit(pd,'UN','NONE','REML')
    if(isTRUE(fun$ok)) {
      unpt<-v9_tillage_point(pd,fun,n)
      unpt$Covariance<-'UN'
      unpt$REML_AIC<-fun$AIC
      unpt$Primary_B_REML_AIC<-f$AIC
      alts[[length(alts)+1L]]<-unpt
    }
    points[[length(points)+1L]]<-pt
    nr<-nr+1L
    cat('TILLAGE ',cr,' ',ss,' | matched ',n,' years | RT-FT ',
        round(pt$RT_minus_FT_Mg_ha,3),' Mg/ha | ',v9_boundary(f),'\n',sep='')
    if(RUN_BOOTSTRAP){
      tag<-paste(cr,ss,sep='_')
      cirows[[length(cirows)+1L]]<-v9_tillage_boot(pd,f,pt,nboot,tag,
                                                    seedoffset=4000000L+100000L*nr)
    }
  }
  v9_save(v9_bind(skip),v9_output('S6_02_simple_matched_year_differences_ALL.csv'))
  v9_save(v9_bind(alts),v9_output('S6_03_UN_covariance_tillage_sensitivity.csv'))
  pts<-v9_bind(points)
  if(nrow(pts)) {
    pts$Holm_p_within_crop<-NA_real_
    for(cr in CROPS) {
      i<-which(pts$Crop==cr)
      if(length(i)) pts$Holm_p_within_crop[i]<-p.adjust(pts$RT_minus_FT_Wald_p[i],method='holm')
    }
  }
  v9_save(pts,v9_output('MAIN_07_post2008_FT_RT_mean_and_stability.csv'))
  v9_save(v9_bind(cirows),v9_output('MAIN_08_post2008_FT_RT_bootstrap_CIs.csv'))
}
