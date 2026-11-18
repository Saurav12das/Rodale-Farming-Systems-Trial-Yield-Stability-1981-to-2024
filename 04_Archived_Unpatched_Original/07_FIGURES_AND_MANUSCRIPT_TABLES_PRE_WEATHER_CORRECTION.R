# Tabular summaries, model diagnostics and figures derived from analysis outputs.
v9_read_result <- function(name,folder='tables') {
  path<-v9_output(name,folder)
  if(!file.exists(path)) return(data.frame())
  tryCatch(read.csv(path,stringsAsFactors=FALSE),error=function(e)data.frame())
}
v9_main_table <- function(){
  p<-v9_read_result('MAIN_01_primary_point_estimates_ALL_definitions.csv')
  cis<-v9_read_result('MAIN_04_bootstrap_mean_SD_CV_95CI.csv')
  if(!nrow(p)) return(data.frame())
  cols<-c('Crop','System','N_observed_years','Trend','Covariance','Serial',
    'Mean_Mg_ha','Detrended_total_SD_Mg_ha','Detrended_total_CV_percent',
    'Latent_residual_SD_Mg_ha','Trend_including_total_SD_Mg_ha')
  p<-p[,cols,drop=FALSE]
  for(metric in c('Mean_Mg_ha','Detrended_total_SD_Mg_ha','Detrended_total_CV_percent')){
    p[[paste0(metric,'_CI_low')]]<-NA_real_
    p[[paste0(metric,'_CI_high')]]<-NA_real_
    if(nrow(cis)) {
      z<-cis[cis$Metric==metric,,drop=FALSE]
      k<-match(paste(p$Crop,p$System),paste(z$Crop,z$System))
      p[[paste0(metric,'_CI_low')]]<-z$CI_low[k]
      p[[paste0(metric,'_CI_high')]]<-z$CI_high[k]
    }
  }
  p$Bootstrap_status<-if(RUN_MODE=='SMOKE') 'PRELIMINARY: NOT FOR MANUSCRIPT' else 'FULL'
  v9_save(p,v9_output('MAIN_TABLE_2_productivity_detrended_stability.csv'))
  p
}
v9_theme <- function() ggplot2::theme_bw(base_size=11)+
  ggplot2::theme(panel.grid.minor=ggplot2::element_blank(),
    plot.title=ggplot2::element_text(face='bold',size=12),
    legend.position='bottom',strip.text=ggplot2::element_text(face='bold'))
v9_colors <- c(CNV='#475569',LEG='#2563EB',MNR='#15803D')
v9_png <- function(g,filename,w=9,h=5.5){
  ggplot2::ggsave(v9_output(filename,'figures'),g,width=w,height=h,dpi=260,
                  bg='white',limitsize=FALSE)
}
v9_make_figures <- function(data){
  library(ggplot2)
  # Crop coverage is obtained from the yield data.
  cov<-v9_read_result('S1_12_system_crop_year_coverage.csv')
  if(nrow(cov)){
    cov$Primary_series<- 'Pre-2008 FT; post-2008 50:50 FT/RT where paired'
    v9_save(cov,v9_output('MAIN_TABLE_1_experimental_observation_coverage.csv'))
  }
  desc<-data$annual
  traj<-v9_read_result('MAIN_03_selected_expected_yield_trajectories.csv')
  tab<-v9_main_table()
  if(nrow(traj)){
    p<-ggplot()+
      geom_point(data=desc,aes(Year,Yield_Mg_ha,color=System),alpha=.6,size=1.2)+
      geom_errorbar(data=desc,aes(Year,ymin=Yield_Mg_ha-1.96*SE_Mg_ha,
                         ymax=Yield_Mg_ha+1.96*SE_Mg_ha,color=System),
                    alpha=.24,width=.12,linewidth=.3)+
      geom_line(data=traj,aes(Year,Expected_yield_at_reference_weather,
                        color=System),linewidth=1)+
      facet_wrap(~Crop,scales='free_y',ncol=1)+
      scale_color_manual(values=v9_colors)+
      labs(title='Observed annual yields and modeled expected productivity',
        subtitle='Lines: selected trends at reference weather; points: treatment means',
        x='Harvest year',y='Yield (Mg/ha)',color='Cropping system')+v9_theme()
    v9_png(p,'MAIN_FIGURE_1_annual_yield_trajectories.png',9,9)
  }
  if(nrow(tab)){
    tab$System<-factor(tab$System,levels=SYSTEMS)
    p<-ggplot(tab,aes(System,Mean_Mg_ha,color=System))+geom_point(size=3)+
      geom_errorbar(aes(ymin=Mean_Mg_ha_CI_low,ymax=Mean_Mg_ha_CI_high),
                    width=.15,linewidth=.7)+facet_wrap(~Crop,scales='free_y',nrow=1)+
      scale_color_manual(values=v9_colors)+
      labs(title='Standardized mean grain productivity, 1981–2024',
           subtitle='Pointwise 95% conditional bootstrap intervals',
           x=NULL,y='Adjusted mean yield (Mg/ha)')+v9_theme()+
      theme(legend.position='none')
    v9_png(p,'MAIN_FIGURE_2_mean_productivity_95CI.png',10,4.3)
    p2<-ggplot(tab,aes(System,Detrended_total_SD_Mg_ha,color=System))+geom_point(size=3)+
      geom_errorbar(aes(ymin=Detrended_total_SD_Mg_ha_CI_low,
                        ymax=Detrended_total_SD_Mg_ha_CI_high),width=.15,linewidth=.7)+
      facet_wrap(~Crop,scales='free_y',nrow=1)+scale_color_manual(values=v9_colors)+
      labs(title='Detrended weather-standardized total temporal SD',
        subtitle='Includes modeled weather variability, excludes fitted productivity trend',
        x=NULL,y='SD (Mg/ha)')+v9_theme()+theme(legend.position='none')
    v9_png(p2,'MAIN_FIGURE_3_detrended_SD_95CI.png',10,4.3)
    cv<-ggplot(tab,aes(System,Detrended_total_CV_percent,color=System))+geom_point(size=3)+
      geom_errorbar(aes(ymin=Detrended_total_CV_percent_CI_low,
                        ymax=Detrended_total_CV_percent_CI_high),width=.15,linewidth=.7)+
      facet_wrap(~Crop,scales='free_y',nrow=1)+scale_color_manual(values=v9_colors)+
      labs(title='Relative detrended yield variability (CV)',x=NULL,y='CV (%)')+
      v9_theme()+theme(legend.position='none')
    v9_png(cv,'SUPP_FIGURE_CV_95CI.png',10,4.3)
  }
  pairs<-v9_read_result('MAIN_05_pairwise_differences_ratios_simultaneous_95CI.csv')
  if(nrow(pairs)){
    z<-pairs[pairs$Metric=='Detrended_total_SD_Mg_ha',,drop=FALSE]
    if(nrow(z)) {
      z$Comparison<-factor(z$Comparison,levels=c('LEG_vs_CNV','MNR_vs_CNV','MNR_vs_LEG'))
      p<-ggplot(z,aes(Estimate,Comparison))+geom_vline(xintercept=1,linetype=2,color='grey50')+
        geom_errorbarh(aes(xmin=Simultaneous_maxT_CI_low,
                            xmax=Simultaneous_maxT_CI_high),height=.19,linewidth=.65)+
        geom_point(size=2.6,color='#0F766E')+facet_wrap(~Crop,nrow=1,scales='free_x')+
        labs(title='Comparisons of detrended yield stability',
          subtitle='Three-comparison bootstrap max-t simultaneous 95% intervals per crop',
          x='Temporal SD ratio (1 = equal SD)',y=NULL)+v9_theme()
      v9_png(p,'MAIN_FIGURE_4_pairwise_SD_ratios_familywise_95CI.png',11,4.4)
    }
  }
  coeff<-v9_read_result('MAIN_02_primary_fixed_effects_weather_coefficients.csv')
  if(nrow(coeff)){
    z<-coeff[coeff$Term%in%c('rain_z','heat_z'),,drop=FALSE]
    z$Predictor<-ifelse(z$Term=='rain_z','Seasonal precipitation (+1 SD)',
      ifelse(z$Crop=='wheat','Apr–Jun temperature (+1 SD)','Summer hot days (+1 SD)'))
    p<-ggplot(z,aes(Predictor,Estimate))+geom_hline(yintercept=0,linetype=2,color='grey50')+
      geom_point(size=2.6,color='#0F766E')+
      geom_errorbar(aes(ymin=Wald_95_Lower,ymax=Wald_95_Upper),width=.12)+
      facet_wrap(~Crop,scales='free_y',nrow=1)+coord_flip()+
      labs(title='Weather associations with modeled grain yield',
        subtitle='Conditional Wald intervals; observational associations, not causal effects',
        x=NULL,y='Yield change (Mg/ha per 1 SD increase)')+v9_theme()
    v9_png(p,'MAIN_FIGURE_5_weather_yield_associations.png',12,4.5)
  }
  till<-v9_read_result('MAIN_07_post2008_FT_RT_mean_and_stability.csv')
  tillci<-v9_read_result('MAIN_08_post2008_FT_RT_bootstrap_CIs.csv')
  if(nrow(till)){
    till$Plot_low<-till$RT_minus_FT_Wald_CI_low
    till$Plot_high<-till$RT_minus_FT_Wald_CI_high
    if(nrow(tillci)) {
      q<-tillci[tillci$Metric=='RT_minus_FT_Mg_ha',,drop=FALSE]
      ix<-match(paste(till$Crop,till$System),paste(q$Crop,q$System))
      till$Plot_low<-q$CI_low[ix];till$Plot_high<-q$CI_high[ix]
    }
    p<-ggplot(till,aes(RT_minus_FT_Mg_ha,System,color=System))+
      geom_vline(xintercept=0,linetype=2,color='grey50')+
      geom_point(size=2.8)+
      geom_errorbarh(aes(xmin=Plot_low,xmax=Plot_high),height=.15,linewidth=.7)+
      facet_wrap(~Crop,scales='free_x',nrow=1)+scale_color_manual(values=v9_colors)+
      labs(title='Post-2008 matched-year tillage comparison',
        subtitle='Reduced/no-till minus full tillage; intervals are model-conditional',
        x='Yield difference (Mg/ha)',y=NULL)+v9_theme()+theme(legend.position='none')
    v9_png(p,'MAIN_FIGURE_6_post2008_tillage_yield_differences.png',11,4.5)
    if(nrow(tillci)) {
      r<-tillci[tillci$Metric=='RT_over_FT_SD_ratio',,drop=FALSE]
      if(nrow(r)){
        pr<-ggplot(r,aes(Estimate,System,color=System))+
          geom_vline(xintercept=1,linetype=2,color='grey50')+
          geom_point(size=2.8)+
          geom_errorbarh(aes(xmin=CI_low,xmax=CI_high),height=.14,linewidth=.7)+
          facet_wrap(~Crop,scales='free_x',nrow=1)+scale_color_manual(values=v9_colors)+
          labs(title='Post-2008 FT/RT temporal SD ratios (matched growing years)',
            subtitle='RT / FT; conditional bootstrap pointwise 95% intervals',
            x='SD ratio (1 = equal SD)',y=NULL)+v9_theme()+theme(legend.position='none')
        v9_png(pr,'SUPP_FIGURE_post2008_tillage_SD_ratios.png',11,4.5)
      }
    }
  }
  # Model-comparison summaries and residual diagnostic figures.
  tr<-v9_read_result('S2_02_ML_trend_comparisons_same_serial_weather.csv')
  if(nrow(tr)){
    z<-tr[tr$Covariance=='B',,drop=FALSE]
    p<-ggplot(z,aes(Trend,Delta_AIC_within_crop_cov,fill=Trend))+geom_col(width=.75)+
      facet_wrap(~Crop,nrow=1,scales='free_y')+
      labs(title='Trend-model sensitivity: ML AIC differences (B covariance)',
        subtitle='Same crop, weather predictors, and AR specification within each comparison',
        y='Δ AIC (smaller is better)',x=NULL)+v9_theme()+
      theme(legend.position='none',axis.text.x=element_text(angle=35,hjust=1))
    v9_png(p,'SUPP_FIGURE_trend_model_AIC.png',11,4.4)
  }
  ac<-v9_read_result('S4_03_REML_covariance_autocorrelation_comparisons.csv')
  if(nrow(ac)){
    ac$Specification<-paste(ac$Covariance,ac$Serial)
    ac<-ac[is.finite(ac$Delta_AIC_within_crop),,drop=FALSE]
    p<-ggplot(ac,aes(Specification,Delta_AIC_within_crop))+geom_col(fill='#64748B')+
      facet_wrap(~Crop,nrow=1,scales='free_y')+
      labs(title='Covariance and AR(1) sensitivity',y='REML Δ AIC',x=NULL)+v9_theme()+
      theme(axis.text.x=element_text(angle=45,hjust=1))
    v9_png(p,'SUPP_FIGURE_covariance_AR1_AIC.png',11,4.4)
  }
  qq<-v9_bind(lapply(CROPS,function(cr) v9_read_result(
    paste0('S4_whitened_residuals_',cr,'.csv'))))
  if(nrow(qq)){
    p<-ggplot(qq,aes(sample=Approx_whitened_residual))+stat_qq(alpha=.55,size=1.1)+
      stat_qq_line(color='#0F766E')+facet_wrap(~Crop,scales='free',nrow=1)+
      labs(title='Approximate whitened-residual Q–Q diagnostics')+v9_theme()
    v9_png(p,'SUPP_FIGURE_whitened_residual_QQ.png',11,4.1)
    p2<-ggplot(qq,aes(Year,Standardized_residual,color=System))+geom_hline(yintercept=0,linetype=2)+
      geom_point(size=1.4,alpha=.65)+geom_line(alpha=.36)+
      facet_wrap(~Crop,ncol=1)+scale_color_manual(values=v9_colors)+
      labs(title='Annual residuals from selected models',y='Standardized residual',x='Year')+v9_theme()
    v9_png(p2,'SUPP_FIGURE_residuals_vs_calendar_year.png',10,8)
  }
  if(nrow(desc)){
    w<-data$annual
    p<-ggplot(w,aes(System,Rain_primary,color=System))+geom_boxplot(outlier.shape=NA)+
      geom_jitter(width=.12,size=.8,alpha=.5)+facet_wrap(~Crop,scales='free_y',nrow=1)+
      scale_color_manual(values=v9_colors)+
      labs(title='Weather exposure in observed crop-growing years',
        x=NULL,y='Prespecified rainfall index (mm)')+v9_theme()+theme(legend.position='none')
    v9_png(p,'SUPP_FIGURE_observed_weather_exposure.png',10,4.3)
  }
  invisible(tab)
}
v9_finish <- function(){
  tabs<-list.files(file.path(V9_OUT,'tables'),pattern='\\.csv$',full.names=FALSE)
  figs<-list.files(file.path(V9_OUT,'figures'),pattern='\\.png$',full.names=FALSE)
  boots<-list.files(file.path(V9_OUT,'bootstrap'),pattern='\\.csv$',full.names=FALSE)
  audit<-data.frame(Artifact_type=c(rep('table',length(tabs)),rep('figure',length(figs)),
                                    rep('bootstrap',length(boots))),
                    File=c(tabs,figs,boots))
  v9_save(audit,file.path(V9_OUT,'00_OUTPUT_MANIFEST.csv'))
  info<-c(paste('FST independent V9 analysis',Sys.time()),
    paste('RUN_MODE:',RUN_MODE),paste('Yield file:',v9_file(YIELD_FILE)),
    paste('PRISM file:',v9_file(WEATHER_FILE)),
    paste('MD5 input checksums:',paste(v9_md5(),collapse=', ')),
    'Primary: corn M1_STEP B AR1; soybean M4_PIECEWISE B NONE; wheat M2_LINEAR B NONE',
    'Primary mean: full 1981–2024 standardized expected productivity at common reference weather.',
    'Primary stability SD: sqrt(latent annual variance + variance of modeled weather response)',
    '  where the fitted long-term expected yield trend does NOT contribute to that SD.',
    'Trend-including SD is a separate sensitivity and includes the 44-year expected trajectory.',
    'AR1 time differences are actual calendar-year distances, not row order.',
    'Weather standardization uses union of observed crop-years by default.',
    'All model and bootstrap inferences are conditional on chosen structure and fixed PRISM data.',
    'WARNING: CNV wheat is sparse; full-period standardized trends require extrapolation.',
    'WARNING: Full-tillage vs reduced tillage comparisons are matched-year, 2008 onward.',
    if(RUN_MODE=='SMOKE') 'ALL BOOTSTRAP CONFIDENCE INTERVALS PRELIMINARY / NOT FOR PUBLICATION' else
        'FULL bootstrap requested: verify success and model-boundary diagnostics before reporting.',
    paste('Tables produced:',length(tabs)),paste('Figures produced:',length(figs)),
    paste('Bootstrap files produced:',length(boots)))
  writeLines(info,file.path(V9_OUT,'00_READ_FIRST.txt'))
  capture.output(sessionInfo(),file=file.path(V9_OUT,'00_R_sessionInfo.txt'))
  cat('\nV9 pipeline complete: ',V9_OUT,'\n',sep='')
  cat('Tables: ',length(tabs),' | Figures: ',length(figs),
      ' | Bootstrap CSV: ',length(boots),'\n',sep='')
}
