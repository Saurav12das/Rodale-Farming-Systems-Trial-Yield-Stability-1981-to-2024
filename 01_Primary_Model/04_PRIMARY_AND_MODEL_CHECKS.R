# Crop models, candidate model comparisons, and diagnostic analyses.
v9_descriptives <- function(data) {
  d<-data$annual
  rows<-lapply(split(d,interaction(d$Crop,d$System,drop=TRUE)),function(z)
    data.frame(Crop=z$Crop[1],System=z$System[1],N_years=nrow(z),
      Mean_Mg_ha=mean(z$Yield_Mg_ha),SD_across_observed_year_means=sd(z$Yield_Mg_ha),
      CV_across_observed_year_means=100*sd(z$Yield_Mg_ha)/mean(z$Yield_Mg_ha),
      Min_yield=min(z$Yield_Mg_ha),Max_yield=max(z$Yield_Mg_ha)))
  v9_save(v9_bind(rows),v9_output('S1_13_raw_observed_crop_system_descriptives.csv'))
  wrows<-list();check<-list()
  for(cr in CROPS){
    r<-v9_ref_weather(data,cr)
    z<-d[d$Crop==cr,,drop=FALSE]
    for(s in SYSTEMS) {
      q<-z[z$System==s,,drop=FALSE]
      wrows[[length(wrows)+1L]]<-data.frame(Crop=cr,System=s,
        N_observed_years=nrow(q),Rain_mean=mean(q$Rain_primary),
        Rain_sd=sd(q$Rain_primary),Heat_mean=mean(q$Heat_primary),
        Heat_sd=sd(q$Heat_primary),
        Reference_rain_outside_observed_fraction=mean(r$Rain_primary<min(q$Rain_primary)|
                                                 r$Rain_primary>max(q$Rain_primary)),
        Reference_heat_outside_observed_fraction=mean(r$Heat_primary<min(q$Heat_primary)|
                                                 r$Heat_primary>max(q$Heat_primary)))
    }
  }
  v9_save(v9_bind(wrows),v9_output('S3_01_weather_exposure_and_reference_support.csv'))
}
v9_primary_and_trends <- function(data) {
  fits<-list(); pds<-list(); modelrows<-list(); trendrows<-list();
  pointrows<-list(); coefs<-list(); trajectories<-list(); corrs<-list();
  for(cr in CROPS){
    sp<-PRIMARY[[cr]]
    pd<-v9_prepare(data,cr,sp$trend,'PRIMARY')
    f<-v9_fit(pd,sp$covariance,sp$serial,'REML',multi=TRUE)
    v9_fit_good(f,paste(cr,'primary'))
    fits[[cr]]<-f;pds[[cr]]<-pd
    cat('PRIMARY ',cr,' | ',sp$trend,' B ',sp$serial,' | REML AIC ',round(f$AIC,3),
        ' | AR1 rho ',round(f$rho,3),' | ',v9_boundary(f),'\n',sep='')
    modelrows[[cr]]<-v9_modelrow(pd,f,'PRIMARY')
    pointrows[[cr]]<-v9_standardized(pd,f)
    coefs[[cr]]<-v9_coefrows(pd,f)
    trajectories[[cr]]<-v9_prediction_trajectory(pd,f)
    corrs[[cr]]<-v9_implied_correlations(pd,f)
    trendrows[[cr]]<-v9_trend_coefficients(pd,f)
    # Correlation between the two seasonal weather predictors.
    corr<-cor(pd$d$rain_z,pd$d$heat_z)
    v9_save(data.frame(Crop=cr,Weather_predictors_R=corr,
      VIF_each_predictor=if(abs(corr)<1) 1/(1-corr^2) else Inf),
      v9_output(paste0('S3_06_weather_collinearity_',cr,'.csv')))
    qq<-data.frame(Crop=cr,System=as.character(pd$d$System),Year=pd$d$Year,
      Yield_Mg_ha=pd$y,Fitted_Mg_ha=as.vector(pd$X%*%f$beta),
      Raw_residual=f$residuals,Residual_scale=sqrt(diag(f$V)),
      Standardized_residual=f$residuals/sqrt(diag(f$V)))
    v9_save(qq,v9_output(paste0('S4_residual_diagnostics_',cr,'.csv')))
    # Whitened residuals are used for diagnostic plots.
    rr<-forwardsolve(t(chol(f$V)),f$residuals)
    qq$Approx_whitened_residual<-as.vector(rr)
    v9_save(qq,v9_output(paste0('S4_whitened_residuals_',cr,'.csv')))
    qq$Era<-ifelse(qq$Year>=ERA_BREAK,'2008-2024','1981-2007')
    bd<-v9_bind(lapply(split(qq,interaction(qq$System,qq$Era,drop=TRUE)),function(z)
      data.frame(Crop=cr,System=z$System[1],Era=z$Era[1],N_annual=nrow(z),
       SD_standardized_residual=if(nrow(z)>=2) sd(z$Standardized_residual) else NA_real_,
       Mean_standardized_residual=mean(z$Standardized_residual))))
    v9_save(bd,v9_output(paste0('S4_residual_variance_by_era_',cr,'.csv')))
    v9_save(v9_lags(pd,f),v9_output(paste0('S4_temporal_residual_lags_',cr,'.csv')))
  }
  v9_save(v9_bind(modelrows),v9_output('MAIN_00_primary_model_audit.csv'))
  v9_save(v9_bind(pointrows),v9_output('MAIN_01_primary_point_estimates_ALL_definitions.csv'))
  v9_save(v9_bind(coefs),v9_output('MAIN_02_primary_fixed_effects_weather_coefficients.csv'))
  v9_save(v9_bind(trajectories),v9_output('MAIN_03_selected_expected_yield_trajectories.csv'))
  v9_save(v9_bind(trendrows),v9_output('S2_01_selected_pre_post_2008_slopes_and_steps.csv'))
  v9_save(v9_bind(corrs),v9_output('S4_01_primary_latent_cross_system_correlations.csv'))
  # Record the fit of the configured crop models.
  verify<-list()
  benchmark<-c(corn=363.9809938,soybean=162.6607930,wheat=98.5967319)
  for(cr in CROPS){
    pd<-v9_prepare(data,cr,'M1_STEP','PRIMARY')
    f<-v9_fit(pd,'B','NONE','REML')
    verify[[cr]]<-data.frame(Crop=cr,Expected_previous_AIC=benchmark[[cr]],
      V9_recomputed_AIC=f$AIC,Difference=f$AIC-benchmark[[cr]])
  }
  v9_save(v9_bind(verify),v9_output('S4_02_2008_step_B_noAR_regression_test.csv'))
  if(any(abs(v9_bind(verify)$Difference)>.2))
    warning('V9 does not closely reproduce previous B noAR AIC: inspect source and formula audit.')
  list(fits=fits,pds=pds)
}
v9_compare_models <- function(data,models) {
  rows<-list(); key<-0L
  # Compare yield-trend specifications by ML with the covariance fixed.
  for(cr in CROPS){
    serial<-PRIMARY[[cr]]$serial
    for(tr in V9_TRENDS) for(cv in c('B','UN')) {
      pd<-v9_prepare(data,cr,tr,'PRIMARY')
      f<-v9_fit(pd,cv,serial,'ML')
      key<-key+1L;rows[[key]]<-v9_modelrow(pd,f,'TREND_ML_FIXED_SERIAL')
    }
  }
  trendtab<-v9_bind(rows)
  trendtab$Delta_AIC_within_crop_cov<-NA_real_
  for(cr in CROPS) for(cv in c('B','UN')) {
    ii<-which(trendtab$Crop==cr&trendtab$Covariance==cv)
    if(length(ii)) trendtab$Delta_AIC_within_crop_cov[ii]<-
      trendtab$AIC[ii]-min(trendtab$AIC[ii],na.rm=TRUE)
  }
  v9_save(trendtab,v9_output('S2_02_ML_trend_comparisons_same_serial_weather.csv'))
  # Compare temporal dependence and covariance specifications for each crop.
  alt<-list()
  for(cr in CROPS){
    pd<-models$pds[[cr]]
    for(cv in c('B','UN','HCS','DIAG')) for(sr in c('NONE','AR1')) {
      f<-v9_fit(pd,cv,sr,'REML')
      alt[[length(alt)+1L]]<-v9_modelrow(pd,f,'COVARIANCE_SERIAL_REML')
      if(isTRUE(f$ok)) {
        v9_save(v9_standardized(pd,f),v9_output(
          sprintf('S4_sensitivity_stability_%s_%s_%s.csv',cr,cv,sr)))
        v9_save(v9_implied_correlations(pd,f),v9_output(
          sprintf('S4_sensitivity_latent_correlations_%s_%s_%s.csv',cr,cv,sr)))
      }
    }
  }
  alt<-v9_bind(alt)
  alt$Delta_AIC_within_crop<-NA_real_
  for(cr in CROPS){
    ii<-which(alt$Crop==cr)
    alt$Delta_AIC_within_crop[ii]<-alt$AIC[ii]-min(alt$AIC[ii],na.rm=TRUE)
  }
  v9_save(alt,v9_output('S4_03_REML_covariance_autocorrelation_comparisons.csv'))
  # Weather model comparisons use maximum likelihood.
  wrows<-list(); wrank<-list()
  for(cr in CROPS){
    sp<-PRIMARY[[cr]]
    f0<-v9_fit(v9_prepare(data,cr,sp$trend,'NONE'),sp$covariance,sp$serial,'ML')
    f1<-v9_fit(v9_prepare(data,cr,sp$trend,'PRIMARY'),sp$covariance,sp$serial,'ML')
    f2<-v9_fit(v9_prepare(data,cr,sp$trend,'INTERACTIONS'),sp$covariance,sp$serial,'ML')
    f3<-v9_fit(v9_prepare(data,cr,sp$trend,'QUADRATIC'),sp$covariance,sp$serial,'ML')
    for(nm in c('NONE','PRIMARY','INTERACTIONS','QUADRATIC')){
      f<-switch(nm,NONE=f0,PRIMARY=f1,INTERACTIONS=f2,QUADRATIC=f3)
      pd<-v9_prepare(data,cr,sp$trend,nm)
      wrows[[length(wrows)+1L]]<-v9_modelrow(pd,f,'WEATHER_ML')
    }
    if(isTRUE(f0$ok)&&isTRUE(f1$ok)&&isTRUE(f2$ok)&&isTRUE(f3$ok)){
      lr01<-max(0,2*(f1$logLik-f0$logLik))
      lr12<-max(0,2*(f2$logLik-f1$logLik))
      lr13<-max(0,2*(f3$logLik-f1$logLik))
      wrank[[cr]]<-data.frame(Crop=cr,Comparison=c('pooled vs no weather',
                              'system interactions vs pooled','quadratic weather vs pooled'),
        LR_stat=c(lr01,lr12,lr13),Df=c(2L,4L,2L),
        Asymptotic_p=c(pchisq(lr01,2,lower.tail=FALSE),pchisq(lr12,4,lower.tail=FALSE),
                       pchisq(lr13,2,lower.tail=FALSE)),
        Notes='Exploratory; asymptotic chi-square, not a causal test')
    }
  }
  v9_save(v9_bind(wrows),v9_output('S3_02_ML_weather_and_interactions_comparison.csv'))
  v9_save(v9_bind(wrank),v9_output('S3_03_weather_model_exploratory_LR_tests.csv'))
  # Assess alternative seasonal weather windows.
  ww<-list()
  for(cr in CROPS){
    sp<-PRIMARY[[cr]]
    for(window in c('primary','alternative')) {
      pd<-v9_prepare(data,cr,sp$trend,'PRIMARY',window=window)
      f<-v9_fit(pd,sp$covariance,sp$serial,'REML')
      ww[[length(ww)+1L]]<-data.frame(Window=window,v9_modelrow(pd,f,'WINDOW_REML'))
      if(isTRUE(f$ok)) v9_save(v9_standardized(pd,f),v9_output(
        paste0('S3_weather_window_',window,'_',cr,'_stability.csv')))
    }
  }
  v9_save(v9_bind(ww),v9_output('S3_04_alternative_seasonal_window_fits.csv'))
  # Estimate the effect of alternative reference-weather years.
  rr<-list()
  for(cr in CROPS){
    sp<-PRIMARY[[cr]]
    pd<-v9_prepare(data,cr,sp$trend,'PRIMARY',reference_mode=ALTERNATIVE_REFERENCE)
    # Refit with weather covariates centered on the alternative reference.
    f<-v9_fit(pd,sp$covariance,sp$serial,'REML')
    if(isTRUE(f$ok)) rr[[cr]]<-v9_standardized(pd,f)
  }
  v9_save(v9_bind(rr),v9_output('S3_05_all_weather_year_reference_stability.csv'))
  # Longitudinal sensitivity using the full-tillage series.
  ft<-list()
  for(cr in CROPS){
    sp<-PRIMARY[[cr]]
    pd<-v9_prepare(data,cr,sp$trend,'PRIMARY',subset_source='FT_ONLY')
    f<-v9_fit(pd,sp$covariance,sp$serial,'REML')
    ft[[length(ft)+1L]]<-v9_modelrow(pd,f,'FT_ONLY_SENSITIVITY')
    if(isTRUE(f$ok)) v9_save(v9_standardized(pd,f),v9_output(
      paste0('S5_FT_only_',cr,'_stability.csv')))
  }
  v9_save(v9_bind(ft),v9_output('S5_01_full_tillage_only_model_audit.csv'))
  # Examine sensitivity to the estimated annual sampling variances.
  virows<-list()
  for(cr in CROPS){
    pd<-models$pds[[cr]];sp<-PRIMARY[[cr]]
    for(k in c(.5,1,2)) {
      pd2<-pd;pd2$vi<-pd$vi*k
      f<-v9_fit(pd2,sp$covariance,sp$serial,'REML')
      virows[[length(virows)+1L]]<-data.frame(Sampling_variance_multiplier=k,
        v9_modelrow(pd2,f,'SAMPLING_VARIANCE_SENSITIVITY'))
    }
  }
  v9_save(v9_bind(virows),v9_output('S5_02_sampling_variance_multiplier_sensitivity.csv'))
  invisible(TRUE)
}
v9_influence <- function(data,models) {
  if(!RUN_LEAVE_ONE_YEAR_OUT) return(invisible(NULL))
  out<-list(); full<-v9_bind(lapply(CROPS,function(cr) v9_standardized(models$pds[[cr]],models$fits[[cr]])))
  for(cr in CROPS){
    pd<-models$pds[[cr]]; sp<-PRIMARY[[cr]]
    years<-sort(unique(pd$d$Year))
    for(yr in years){
      sub<-pd
      ix<-which(pd$d$Year!=yr)
      sub$d<-pd$d[ix,,drop=FALSE]
      sub$X<-pd$X[ix,,drop=FALSE]
      sub$y<-pd$y[ix];sub$vi<-pd$vi[ix];sub$sys<-pd$sys[ix]
      sub$gaps<-pd$gaps[ix,ix,drop=FALSE]
      sub$n<-length(ix)
      sub$logdet_XtX<-2*sum(log(diag(chol(crossprod(sub$X)))))
      ff<-v9_fit(sub,sp$covariance,sp$serial,'REML',warm=models$fits[[cr]]$par,fast=TRUE,multi=FALSE)
      if(!isTRUE(ff$ok)||ff$convergence!=0)
         ff<-v9_fit(sub,sp$covariance,sp$serial,'REML',warm=models$fits[[cr]]$par,fast=FALSE,multi=TRUE)
      if(isTRUE(ff$ok)) {
        q<-v9_standardized(sub,ff)
        base<-full[full$Crop==cr,,drop=FALSE]
        q$Removed_year<-yr
        q$Delta_mean<-q$Mean_Mg_ha-base$Mean_Mg_ha[match(q$System,base$System)]
        q$Delta_detrended_SD<-q$Detrended_total_SD_Mg_ha-
          base$Detrended_total_SD_Mg_ha[match(q$System,base$System)]
        q$Converged<-ff$convergence==0
        out[[length(out)+1L]]<-q
      } else out[[length(out)+1L]]<-data.frame(Crop=cr,System=SYSTEMS,
                  Removed_year=yr,Converged=FALSE)
    }
  }
  v9_save(v9_bind(out),v9_output('S5_03_leave_one_calendar_year_out_main_effects.csv'))
}
