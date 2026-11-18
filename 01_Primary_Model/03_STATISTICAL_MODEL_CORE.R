# Gaussian likelihood functions for crop-specific annual yield models.
# Sampling variances, cross-system covariance and calendar-year dependence are explicit.
V9_TRENDS <- c('M1_STEP','M2_LINEAR','M3_STEP_LINEAR','M4_PIECEWISE')
v9_formula <- function(trend,weather='PRIMARY',arms=FALSE) {
  if(arms) return(~0+Arm)
  if(!trend%in%V9_TRENDS) stop('Unknown trend: ',trend)
  tx <- switch(trend,
    M1_STEP='0+System+System:Post2008',
    M2_LINEAR='0+System+System:year_decade',
    M3_STEP_LINEAR='0+System+System:Post2008+System:year_decade',
    M4_PIECEWISE='0+System+System:year_decade+System:hinge_decade')
  wx <- switch(weather,
    PRIMARY='+rain_z+heat_z',
    INTERACTIONS='+System:rain_z+System:heat_z',
    QUADRATIC='+rain_z+heat_z+I(rain_z^2)+I(heat_z^2)',
    NONE='',stop('Bad weather model: ',weather))
  as.formula(paste('~',tx,wx))
}
v9_ref_weather <- function(data,cr,mode=REFERENCE_WEATHER) {
  d<-data$annual[data$annual$Crop==cr,,drop=FALSE]
  w<-data$weather[data$weather$Crop==cr & data$weather$WeatherWindowComplete,,drop=FALSE]
  if(mode=='OBSERVED_CROP_YEARS') w<-w[w$Year %in% unique(d$Year),,drop=FALSE]
  else if(mode!='ALL_ELIGIBLE_WEATHER_YEARS') stop('Unknown reference mode: ',mode)
  if(anyNA(w$Rain_primary)|anyNA(w$Heat_primary)|!nrow(w)) stop('Incomplete reference weather')
  if(anyDuplicated(w$Year)) stop('Duplicate weather years')
  w[order(w$Year),,drop=FALSE]
}
v9_predictor_data <- function(df,cr,center_r,center_h,rain='Rain_primary',heat='Heat_primary') {
  df$Year<-as.numeric(df$Year)
  df$Post2008<-as.numeric(df$Year>=ERA_BREAK)
  df$year_decade<-(df$Year-2002.5)/10
  df$hinge_decade<-pmax(0,df$Year-ERA_BREAK)/10
  df$rain_z<-(df[[rain]]-center_r)/100
  df$heat_z<-(df[[heat]]-center_h)/(if(cr=='wheat') 1 else 10)
  if('System'%in%names(df)) df$System<-factor(df$System,levels=SYSTEMS)
  df
}
v9_prepare <- function(data,cr,trend,weather='PRIMARY',reference_mode=REFERENCE_WEATHER,
                       window='primary',subset_source='COMPOSITE') {
  if(subset_source=='COMPOSITE') d<-data$annual[data$annual$Crop==cr,,drop=FALSE]
  else if(subset_source=='FT_ONLY') {
    d<-data$arms[data$arms$Crop==cr & data$arms$Tillage=='till',,drop=FALSE]
    d<-merge(d,data$weather,by=c('Crop','Year'),all.x=TRUE,sort=FALSE)
    d$Post2008<-as.integer(d$Year>=ERA_BREAK)
    if(anyNA(d$WeatherWindowComplete)|any(!d$WeatherWindowComplete))
      stop('Incomplete FT-only weather')
  } else stop('Unknown source: ',subset_source)
  d<-d[order(d$Year,match(d$System,SYSTEMS)),,drop=FALSE]
  rownames(d)<-NULL
  ref<-v9_ref_weather(data,cr,reference_mode)
  if(window=='primary') {rain<-'Rain_primary';heat<-'Heat_primary'}
  else if(window=='alternative') {rain<-'Rain_alternative';heat<-'Heat_alternative'}
  else stop('Unknown weather window')
  # Center weather indicators using the shared crop-level reference years.
  rc<-mean(ref[[rain]]); hc<-mean(ref[[heat]])
  d<-v9_predictor_data(d,cr,rc,hc,rain,heat)
  ref<-v9_predictor_data(ref,cr,rc,hc,rain,heat)
  form<-v9_formula(trend,weather)
  X<-stats::model.matrix(form,data=d)
  if(qr(X)$rank<ncol(X)) stop(cr,' ',trend,' ',weather,': rank-deficient fixed-effect design')
  if(nrow(d)<=ncol(X)+3L) stop(cr,': too few annual means for fitted predictors')
  if(any(!is.finite(d$Yield_Mg_ha)|!is.finite(d$vi))||any(d$vi<=0)) stop('Bad annual yield/vi')
  n<-nrow(d); p<-ncol(X)
  list(crop=cr,trend=trend,weather=weather,window=window,source=subset_source,
    d=d,ref=ref,year_grid=YEAR_RANGE,formula=form,
    X=X,y=d$Yield_Mg_ha,vi=d$vi,sys=match(as.character(d$System),SYSTEMS),
    sys_names=SYSTEMS,ng=length(SYSTEMS),n=n,p=p,
    gaps=abs(outer(as.numeric(d$Year),as.numeric(d$Year),'-')),
    logdet_XtX=2*sum(log(diag(chol(crossprod(X))))),
    rain_center=rc,heat_center=hc,rain_col=rain,heat_col=heat,
    reference_mode=reference_mode)
}
v9_sigma <- function(theta,kind,k) {
  if(kind=='B') {
    a<-exp(theta[seq_len(k+1)])
    S<-matrix(a[1],k,k)
    diag(S)<-diag(S)+a[seq_len(k)+1L]
    return(S)
  }
  if(kind=='B_EQUAL') {
    a<-exp(theta[1:2]);S<-matrix(a[1],k,k);diag(S)<-diag(S)+a[2]
    return(S)
  }
  if(kind=='DIAG') return(diag(exp(theta[seq_len(k)]),k))
  if(kind=='HCS') {
    ss<-exp(theta[seq_len(k)])
    lo<- -1/(k-1)+1e-5; hi<-.99999
    rho<-lo+(hi-lo)*plogis(theta[k+1L])
    C<-matrix(rho,k,k);diag(C)<-1
    return((ss%o%ss)*C)
  }
  if(kind=='UN') {
    L<-matrix(0,k,k)
    loc<-which(lower.tri(L,diag=TRUE),arr.ind=TRUE)
    if(nrow(loc)!=length(theta[seq_len(nrow(loc))])) stop('Bad UN parameter count')
    for(i in seq_len(nrow(loc))) {
      r<-loc[i,1];c<-loc[i,2]
      L[r,c]<-if(r==c) exp(theta[i]) else theta[i]
    }
    return(tcrossprod(L))
  }
  stop('Unknown covariance type ',kind)
}
v9_ncov <- function(kind,k) switch(kind,B=k+1L,B_EQUAL=2L,DIAG=k,HCS=k+1L,UN=k*(k+1L)/2L,
                                    stop('Unknown covariance type'))
v9_V <- function(pd,theta,kind,serial) {
  nc<-v9_ncov(kind,pd$ng)
  S<-v9_sigma(theta[seq_len(nc)],kind,pd$ng)
  if(serial=='AR1') {
    rho<-tanh(theta[nc+1L])
    K<-rho^pd$gaps
  } else if(serial=='NONE') K<-1*(pd$gaps==0) else stop('Bad serial')
  V<-S[pd$sys,pd$sys]*K
  diag(V)<-diag(V)+pd$vi
  V
}
v9_likelihood <- function(pd,theta,kind,serial,method='REML',details=FALSE) {
  tryCatch({
    V<-v9_V(pd,theta,kind,serial)
    R<-chol(V)
    VX<-backsolve(R,forwardsolve(t(R),pd$X))
    Vy<-backsolve(R,forwardsolve(t(R),pd$y))
    XX<-crossprod(pd$X,VX)
    R2<-chol(XX)
    beta<-as.vector(backsolve(R2,forwardsolve(t(R2),crossprod(pd$X,Vy))))
    names(beta)<-colnames(pd$X)
    ee<-as.vector(pd$y-pd$X%*%beta)
    Ve<-backsolve(R,forwardsolve(t(R),ee))
    quadratic<-sum(ee*Ve)
    ldV<-2*sum(log(diag(R)))
    if(method=='REML') neg2ll<-(pd$n-pd$p)*log(2*pi)+ldV+
      2*sum(log(diag(R2)))+quadratic-pd$logdet_XtX
    else if(method=='ML') neg2ll<-pd$n*log(2*pi)+ldV+quadratic
    else stop('Invalid likelihood method')
    if(!is.finite(neg2ll)) stop('nonfinite log likelihood')
    if(!details) return(neg2ll/2)
    list(nll=neg2ll/2,logLik=-neg2ll/2,beta=beta,XX=XX,S=v9_sigma(
      theta[seq_len(v9_ncov(kind,pd$ng))],kind,pd$ng),V=V,residuals=ee)
  },error=function(e){
    if(details) stop(conditionMessage(e))
    1e12
  })
}
v9_initial <- function(pd,kind) {
  b<-as.vector(qr.solve(pd$X,pd$y))
  ee<-as.vector(pd$y-pd$X%*%b)
  vg<-vapply(seq_len(pd$ng),function(j){
    ix<-which(pd$sys==j)
    if(length(ix)<2L) return(.08)
    max(.05,var(ee[ix])-mean(pd$vi[ix]))
  },numeric(1))
  sh<-max(.02,.5*min(vg))
  un<-pmax(.02,vg-sh)
  if(kind=='B') return(log(c(sh,un)))
  if(kind=='B_EQUAL') return(log(c(sh,mean(un))))
  if(kind=='DIAG') return(log(vg))
  if(kind=='HCS') {
    lo<- -1/(pd$ng-1)+1e-5;hi<-.99999
    frac<-(.5-lo)/(hi-lo)
    return(c(log(sqrt(vg)),qlogis(frac)))
  }
  if(kind=='UN') {
    L<-t(chol(matrix(sh,pd$ng,pd$ng)+diag(un,pd$ng)))
    loc<-which(lower.tri(L,diag=TRUE),arr.ind=TRUE)
    return(vapply(seq_len(nrow(loc)),function(i) {
      r<-loc[i,1];c<-loc[i,2]
      if(r==c) log(L[r,c]) else L[r,c]
    },numeric(1)))
  }
  stop('Unknown kind')
}
v9_bounds <- function(kind,k,serial) {
  if(kind=='HCS'){
    lo<-c(rep(-14,k),-8);hi<-c(rep(12,k),8)
  } else if(kind=='UN'){
    idx<-which(lower.tri(matrix(0,k,k),diag=TRUE),arr.ind=TRUE)
    isdiag<-idx[,1]==idx[,2]
    lo<-ifelse(isdiag,-14,-100);hi<-ifelse(isdiag,12,100)
  } else {lo<-rep(-20,v9_ncov(kind,k));hi<-rep(12,length(lo))}
  if(serial=='AR1') {lo<-c(lo,atanh(-.98));hi<-c(hi,atanh(.98))}
  list(lo=lo,hi=hi)
}
v9_fit <- function(pd,kind='B',serial='NONE',method='REML',
                   warm=NULL,fast=FALSE,multi=TRUE) {
  nc<-v9_ncov(kind,pd$ng)
  if(!is.null(warm)&&length(warm)>=nc) starts0<-head(as.numeric(warm),nc)
  else starts0<-v9_initial(pd,kind)
  lim<-v9_bounds(kind,pd$ng,serial)
  if(serial=='AR1') {
    rs<-if(fast) c(.25,0) else if(multi) c(0,.2,.5,-.2,.75) else c(.25)
    if(!is.null(warm)&&length(warm)==nc+1L)
      rs<-unique(c(tanh(tail(warm,1)),rs))
    starts<-lapply(rs,function(r) c(starts0,atanh(r)))
  } else {
    starts<-list(starts0)
    if(!fast&&multi&&kind=='B'){
      a<-starts0;a[1]<-log(.05);starts[[2]]<-a
    }
  }
  opts<-list()
  for(st in starts){
    st<-pmin(lim$hi-1e-7,pmax(lim$lo+1e-7,st))
    q<-tryCatch(suppressWarnings(optim(st,
       fn=function(theta) v9_likelihood(pd,theta,kind,serial,method),
       method='L-BFGS-B',lower=lim$lo,upper=lim$hi,
       control=list(maxit=if(fast) 300 else 650,factr=3e6))),error=function(e)NULL)
    if(!is.null(q)&&is.finite(q$value)&&q$value<1e11)
      opts[[length(opts)+1L]]<-q
  }
  if(!length(opts)) return(list(ok=FALSE,kind=kind,serial=serial,
                                method=method,reason='No finite fit'))
  good<-which(vapply(opts,function(q)q$convergence==0,logical(1)))
  ch<-if(length(good)) good else seq_along(opts)
  best<-opts[[ch[which.min(vapply(opts[ch],function(q)q$value,numeric(1)))]]]
  deta<-tryCatch(v9_likelihood(pd,best$par,kind,serial,method,TRUE),
                 error=function(e)NULL)
  if(is.null(deta)) return(list(ok=FALSE,kind=kind,serial=serial,
                                method=method,reason='Final V not positive definite'))
  Vb<-chol2inv(chol(deta$XX))
  dimnames(Vb)<-list(names(deta$beta),names(deta$beta))
  list(ok=TRUE,kind=kind,serial=serial,method=method,
    par=best$par,logLik=deta$logLik,AIC=-2*deta$logLik+2*(length(best$par)+pd$p),
    beta=deta$beta,beta_cov=Vb,beta_se=sqrt(pmax(0,diag(Vb))),
    S=deta$S,V=deta$V,residuals=deta$residuals,
    rho=if(serial=='AR1') tanh(tail(best$par,1)) else 0,
    convergence=best$convergence,optim_message=if(is.null(best$message)) '' else best$message,
    n_starts=length(opts),n_parameter=length(best$par)+pd$p)
}
v9_boundary <- function(f) {
  if(!isTRUE(f$ok)) return(paste('FAILED',f$reason))
  x<-character()
  if(f$convergence!=0) x<-c(x,paste('optimizer',f$convergence))
  if(f$kind=='B') {
    a<-exp(head(f$par,nrow(f$S)+1L))
    nm<-c('shared',paste0(seq_len(nrow(f$S)),'_specific'))
    if(min(a)<1e-5) x<-c(x,paste('near_zero',nm[which.min(a)]))
  }
  if(f$kind=='UN'&&any(abs(cov2cor(f$S)[lower.tri(f$S)])>.97)) x<-c(x,'UN_corr_gt_0.97')
  if(f$kind=='HCS'&&any(abs(cov2cor(f$S)[lower.tri(f$S)])>.97)) x<-c(x,'HCS_corr_gt_0.97')
  if(f$serial=='AR1'&&abs(f$rho)>.95) x<-c(x,'AR1_rho_gt_0.95')
  if(!length(x)) 'None' else paste(x,collapse='; ')
}
v9_modelrow <- function(pd,f,stage='') {
  data.frame(Stage=stage,Crop=pd$crop,Trend=pd$trend,Weather=pd$weather,
    Covariance=f$kind,Serial=f$serial,Method=f$method,N_annual=pd$n,
    N_distinct_years=length(unique(pd$d$Year)),
    N_fixed=pd$p,AIC=if(isTRUE(f$ok)) f$AIC else NA_real_,
    logLik=if(isTRUE(f$ok)) f$logLik else NA_real_,
    AR1_rho=if(isTRUE(f$ok)&&f$serial=='AR1') f$rho else NA_real_,
    Convergence=if(isTRUE(f$ok)) f$convergence else NA_integer_,
    Boundary=if(isTRUE(f$ok)) v9_boundary(f) else f$reason,
    stringsAsFactors=FALSE)
}
v9_fit_good <- function(f,label='fit') {
  if(!isTRUE(f$ok)||f$convergence!=0) stop(label,' failed to converge; inspect model diagnostics')
  invisible(f)
}
v9_design <- function(pd,newdata) {
  z<-model.matrix(pd$formula,data=newdata)
  if(!setequal(colnames(z),colnames(pd$X))) stop('Prediction design mismatch')
  z[,colnames(pd$X),drop=FALSE]
}
v9_standardized <- function(pd,f) {
  if(!isTRUE(f$ok)) return(data.frame())
  ans<-list()
  for(s in pd$sys_names) {
    ys<-data.frame(Year=pd$year_grid,System=s,
      Rain_primary=pd$rain_center,Heat_primary=pd$heat_center,
      Rain_alternative=pd$rain_center,Heat_alternative=pd$heat_center)
    ys<-v9_predictor_data(ys,pd$crop,pd$rain_center,pd$heat_center,pd$rain_col,pd$heat_col)
    trend_curve<-as.vector(v9_design(pd,ys)%*%f$beta)
    r<-pd$ref
    r$System<-s
    r$Year<-2002
    r<-v9_predictor_data(r,pd$crop,pd$rain_center,pd$heat_center,pd$rain_col,pd$heat_col)
    Xw<-v9_design(pd,r)
    r$rain_z<-0;r$heat_z<-0
    X0<-v9_design(pd,r)
    weather_effect<-as.vector((Xw-X0)%*%f$beta)
    vw<-mean((weather_effect-mean(weather_effect))^2)
    vt<-mean((trend_curve-mean(trend_curve))^2)
    mean_y<-mean(trend_curve)+mean(weather_effect)
    j<-match(s,pd$sys_names)
    vlat<-f$S[j,j]
    sd_det<-sqrt(max(0,vlat+vw))
    sd_full<-sqrt(max(0,vlat+vw+vt))
    ans[[length(ans)+1L]]<-data.frame(Crop=pd$crop,System=s,
      Trend=pd$trend,Covariance=f$kind,Serial=f$serial,Weather=pd$weather,
      N_observed_years=sum(pd$d$System==s),
      Standardization_years='1981-2024',Weather_reference=pd$reference_mode,
      Mean_Mg_ha=mean_y,
      Latent_residual_SD_Mg_ha=sqrt(max(0,vlat)),
      Weather_variance_Mg2_ha2=vw,Trend_variance_Mg2_ha2=vt,
      Detrended_total_SD_Mg_ha=sd_det,
      Detrended_total_CV_percent=100*sd_det/mean_y,
      Trend_including_total_SD_Mg_ha=sd_full,
      Trend_including_total_CV_percent=100*sd_full/mean_y,
      AR1_rho=if(f$serial=='AR1') f$rho else NA_real_)
  }
  v9_bind(ans)
}
v9_prediction_trajectory <- function(pd,f) {
  rr<-list()
  for(s in pd$sys_names){
    t<-data.frame(Year=YEAR_RANGE,System=s,
      Rain_primary=pd$rain_center,Heat_primary=pd$heat_center,
      Rain_alternative=pd$rain_center,Heat_alternative=pd$heat_center)
    t<-v9_predictor_data(t,pd$crop,pd$rain_center,pd$heat_center,pd$rain_col,pd$heat_col)
    MM<-v9_design(pd,t)
    pred<-as.vector(MM%*%f$beta)
    se<-sqrt(pmax(0,rowSums((MM%*%f$beta_cov)*MM)))
    rr[[s]]<-data.frame(Crop=pd$crop,System=s,Year=t$Year,
      Expected_yield_at_reference_weather=pred,Approx_Wald_SE=se,
      CI_low=pred-1.96*se,CI_high=pred+1.96*se)
  }
  v9_bind(rr)
}
v9_coefrows <- function(pd,f) {
  if(!isTRUE(f$ok)) return(data.frame())
  data.frame(Crop=pd$crop,Trend=pd$trend,Covariance=f$kind,Serial=f$serial,
    Term=names(f$beta),Estimate=as.numeric(f$beta),Wald_SE=as.numeric(f$beta_se),
    Wald_95_Lower=as.numeric(f$beta-1.96*f$beta_se),
    Wald_95_Upper=as.numeric(f$beta+1.96*f$beta_se))
}
v9_trend_coefficients <- function(pd,f){
  out<-list()
  for(s in pd$sys_names){
    makeX<-function(y){
      t<-data.frame(Year=y,System=s,
        Rain_primary=pd$rain_center,Heat_primary=pd$heat_center,
        Rain_alternative=pd$rain_center,Heat_alternative=pd$heat_center)
      t<-v9_predictor_data(t,pd$crop,pd$rain_center,pd$heat_center,pd$rain_col,pd$heat_col)
      as.vector(v9_design(pd,t)[1,])
    }
    slope_pre<-10*(makeX(2000)-makeX(1999))
    slope_post<-10*(makeX(2019)-makeX(2018))
    jump<-makeX(2008)-makeX(2007)-slope_pre/10
    for(qn in c('Pre2008_slope_per_decade','Post2008_slope_per_decade','Level_jump_at_2008')) {
      vec<-switch(qn,Pre2008_slope_per_decade=slope_pre,
        Post2008_slope_per_decade=slope_post,Level_jump_at_2008=jump)
      est<-sum(vec*f$beta)
      se<-sqrt(max(0,as.numeric(t(vec)%*%f$beta_cov%*%vec)))
      out[[length(out)+1L]]<-data.frame(Crop=pd$crop,System=s,Quantity=qn,
        Estimate=est,Wald_SE=se,Wald_low=est-1.96*se,Wald_high=est+1.96*se,
        Fixed_by_model=all(abs(vec)<1e-10))
    }
  }
  v9_bind(out)
}
v9_lags <- function(pd,f,maxlag=8L){
  if(!isTRUE(f$ok)) return(data.frame())
  ans<-list()
  sr<-f$residuals/sqrt(diag(f$V))
  for(s in pd$sys_names){
    ii<-which(pd$d$System==s)
    yy<-pd$d$Year[ii];zz<-sr[ii]
    for(lg in seq_len(maxlag)) {
      a<-which((yy+lg)%in%yy);b<-match(yy[a]+lg,yy)
      n<-length(a)
      corr<-if(n>=4L && sd(zz[a])>0 && sd(zz[b])>0) cor(zz[a],zz[b]) else NA_real_
      ans[[length(ans)+1L]]<-data.frame(Crop=pd$crop,System=s,
        Trend=pd$trend,Covariance=f$kind,Serial=f$serial,
        Calendar_lag_years=lg,N_pairs=n,Residual_correlation=corr)
    }
  }
  v9_bind(ans)
}
v9_implied_correlations <- function(pd,f){
  S<-f$S;C<-cov2cor(S)
  rows<-list()
  for(i in seq_len(pd$ng)) for(j in i:pd$ng){
    rows[[length(rows)+1L]]<-data.frame(Crop=pd$crop,
      System1=pd$sys_names[i],System2=pd$sys_names[j],
      Same_year_latent_covariance=S[i,j],Same_year_latent_correlation=C[i,j],
      Model=f$kind,Serial=f$serial,AR1_rho=f$rho)
  }
  v9_bind(rows)
}
v9_simulate_y <- function(pd,f) {
  U<-chol(f$V)
  as.vector(pd$X%*%f$beta+t(U)%*%rnorm(pd$n))
}
v9_pairwise <- function(st) {
  if(!nrow(st)) return(data.frame())
  pairs<-list(c('LEG','CNV'),c('MNR','CNV'),c('MNR','LEG'))
  out<-list()
  metrics<-c('Mean_Mg_ha','Detrended_total_SD_Mg_ha','Detrended_total_CV_percent',
             'Latent_residual_SD_Mg_ha','Trend_including_total_SD_Mg_ha',
             'Trend_including_total_CV_percent')
  for(pp in pairs){
    a<-st[st$System==pp[1],,drop=FALSE]
    b<-st[st$System==pp[2],,drop=FALSE]
    for(metric in metrics){
      if(metric=='Mean_Mg_ha'){
        val<-a[[metric]]-b[[metric]];type<-'difference';null<-0
      } else {
        val<-a[[metric]]/b[[metric]];type<-'ratio';null<-1
      }
      out[[length(out)+1L]]<-data.frame(Crop=a$Crop,Comparison=paste(pp,collapse='_vs_'),
        Metric=metric,Type=type,Null=null,Estimate=val)
    }
  }
  v9_bind(out)
}
