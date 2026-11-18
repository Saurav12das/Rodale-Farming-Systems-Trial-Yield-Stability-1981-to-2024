# Publication-style figures for the FST 1981--2024 yield stability manuscript.
# R >=4.2, ggplot2 >=3.4; run via 00_RUN_ALL_FIGURES.R or:
# OUTPUT_DIR <- "path/to/FST_YieldStability_R_Figure_Supplement"
# source("02_MAKE_FIGURES.R"); fst_draw_figures()
# Figures contain observed-descriptive quantities except the mean--SD panel,
# which uses estimates from the primary REML model; model fitting and confidence
# intervals are generated in the accompanying statistical-analysis code.
# No inferential bootstrap claims are generated in this figure script.

FST_COLORS <- c(CNV="#465971", LEG="#C77B36", MNR="#15887C")
FST_CROPS <- c("corn","soybean","wheat")
FST_SYSTEMS <- c("CNV","LEG","MNR")
fst_d <- function(filename) {
  path<-file.path(OUTPUT_DIR,"Derived_data",filename)
  if(!file.exists(path)) stop("Missing derived data: ",path,
    "\nRun 00_RUN_ALL_FIGURES.R using raw input data or supply the included Derived_data CSVs")
  read.csv(path, stringsAsFactors=FALSE, check.names=FALSE)
}
fst_figure_theme <- function() {
  ggplot2::theme_classic(base_size=10.5) +
    ggplot2::theme(strip.background=ggplot2::element_rect(fill="#F4F6F7",color="#D6DDE3",linewidth=.35),
      strip.text=ggplot2::element_text(face="bold",size=10),
      panel.grid.major.y=ggplot2::element_line(color="#E5EAEE",linewidth=.23),
      panel.grid.major.x=ggplot2::element_blank(),
      legend.position="bottom",legend.title=ggplot2::element_blank(),
      legend.key=ggplot2::element_blank(),
      axis.title=ggplot2::element_text(size=10),
      axis.text=ggplot2::element_text(size=9),
      plot.margin=ggplot2::margin(10,14,9,12))
}
fst_save_plot <- function(plot, stem, width=11, height=5) {
  out<-file.path(OUTPUT_DIR,"Figures")
  dir.create(out,recursive=TRUE,showWarnings=FALSE)
  # Vector PDF for journals and 600 dpi LZW TIFF for systems that require raster.
  ggplot2::ggsave(file.path(out,paste0(stem,".pdf")),plot=plot,
                  width=width,height=height,units="in",device="pdf",useDingbats=FALSE)
  ggplot2::ggsave(file.path(out,paste0(stem,".tiff")),plot=plot,
                  width=width,height=height,units="in",device="tiff",dpi=600,compression="lzw")
  ggplot2::ggsave(file.path(out,paste0(stem,".png")),plot=plot,
                  width=width,height=height,units="in",device="png",dpi=300,bg="white")
  invisible(stem)
}

fst_draw_figures <- function() {
  if(!requireNamespace("ggplot2",quietly=TRUE)) stop("Install ggplot2 first")
  library(ggplot2)
  theme_set(fst_figure_theme())
  crops<-c(corn="Corn",soybean="Soybean",wheat="Wheat")

  # F1. Actual harvest-year support, distinguishing unpaired post-2008 arms.
  coverage<-fst_d("01_harvest_year_coverage.csv")
  coverage$Crop<-factor(coverage$Crop,levels=FST_CROPS,labels=crops)
  coverage$System<-factor(coverage$System,levels=rev(FST_SYSTEMS))
  coverage$Status<-factor(coverage$Status,
    levels=c("Included: FT","Included: FT+RT","One arm only: excluded","No harvest"))
  statcols<-c("Included: FT"="#7188A2","Included: FT+RT"="#1B8B7D",
              "One arm only: excluded"="#DDA15E","No harvest"="#E6EBEF")
  f1<-ggplot(coverage,aes(Year,System,fill=Status))+
    geom_tile(width=.94,height=.79,color="white",linewidth=.18)+
    geom_vline(xintercept=2007.5,linetype="dashed",color="#697784",linewidth=.45)+
    facet_grid(Crop~.,switch="y")+
    scale_fill_manual(values=statcols,drop=FALSE)+
    scale_x_continuous(breaks=c(1981,1985,1990,1995,2000,2005,2010,2015,2020,2024))+
    labs(x="Calendar harvest year",y=NULL)+theme(panel.grid=element_blank(),
      strip.placement="outside",legend.position="bottom")
  fst_save_plot(f1,"Figure_1_crop_year_coverage",12.3,5.5)

  # Figure 3. Yield gaps within overlapping harvest-years only (CNV minus organic).
  gap<-fst_d("03_system_gap_by_year.csv")
  sm<-fst_d("02_system_gap_summaries.csv")
  gap$Crop<-factor(gap$Crop,levels=FST_CROPS,labels=crops)
  sm$Crop<-factor(sm$Crop,levels=FST_CROPS,labels=crops)
  gap$System<-factor(gap$System,levels=c("LEG","MNR"))
  sm$System<-factor(sm$System,levels=c("LEG","MNR"))
  sm$Gap_CNV_minus_org<- -sm$Organic_minus_CNV
  sm$Label<-sprintf("Gap %.1f%% | n=%d",sm$Gap_pct,sm$N_years)
  # Set panel-specific annotation positions based on actual data spread.
  panelmax<-aggregate(Gap_CNV_minus_organic~Crop+System,gap,max)
  sm<-merge(sm,panelmax,by=c("Crop","System"),all.x=TRUE)
  sm$LabelY<-sm$Gap_CNV_minus_organic
  f2<-ggplot(gap,aes(Year,Gap_CNV_minus_organic,color=System))+
    geom_hline(yintercept=0,linetype="dashed",color="#8B949B",linewidth=.4)+
    geom_point(size=1.65,alpha=.75)+
    geom_hline(data=sm,aes(yintercept=Gap_CNV_minus_org,color=System),linetype="dotted",linewidth=.55)+
    geom_text(data=sm,aes(x=1982,y=LabelY,label=Label,color=System),
              inherit.aes=FALSE,hjust=0,vjust=-.25,size=2.6,show.legend=FALSE)+
    facet_grid(Crop~System,scales="free_y")+
    scale_color_manual(values=FST_COLORS[c("LEG","MNR")],drop=FALSE)+
    labs(x="Shared harvest year",y="CNV minus organic grain yield (Mg ha^-1)")+
    theme(legend.position="none")

  fst_save_plot(f2,"Figure_3_shared_year_yield_gaps",12,7.7)

  # Figure 2. Model-derived productivity--stability fingerprint, with iso-CV guides.
  model<-fst_d("12_primary_model_point_estimates.csv")
  model$Crop<-factor(model$Crop,levels=FST_CROPS,labels=crops)
  model$System<-factor(model$System,levels=FST_SYSTEMS)
  guide<-do.call(rbind,lapply(levels(model$Crop),function(cr){
    z<-model[model$Crop==cr,]; x0<-max(.01,min(z$Mean_Mg_ha)*.82);x1<-max(z$Mean_Mg_ha)*1.13
    expand.grid(Crop=cr,Mean_Mg_ha=seq(x0,x1,length.out=80),CV=c(15,20,25,30,35,40))
  }))
  guide$Crop<-factor(guide$Crop,levels=levels(model$Crop))
  guide$Guide_SD<-guide$Mean_Mg_ha*guide$CV/100
  f3<-ggplot()+
    geom_line(data=guide,aes(Mean_Mg_ha,Guide_SD,group=interaction(Crop,CV)),
              color="#C2CAD1",linewidth=.42)
  if(any(is.finite(model$SD_CI_low)&is.finite(model$SD_CI_high)))
    f3<-f3+geom_errorbar(data=model,aes(x=Mean_Mg_ha,ymin=SD_CI_low,ymax=SD_CI_high,
              color=System),width=.035,linewidth=.55,alpha=.65)
  if(any(is.finite(model$Mean_CI_low)&is.finite(model$Mean_CI_high)))
    f3<-f3+geom_segment(data=model,aes(x=Mean_CI_low,xend=Mean_CI_high,
                 y=SD_Mg_ha,yend=SD_Mg_ha,color=System),linewidth=.55,alpha=.65)
  f3<-f3+geom_point(data=model,aes(Mean_Mg_ha,SD_Mg_ha,color=System),size=3.4)+
    geom_text(data=model,aes(Mean_Mg_ha,SD_Mg_ha,label=System),
              nudge_x=.07,nudge_y=.035,hjust=0,size=3.1,show.legend=FALSE)+
    facet_wrap(~Crop,nrow=1,scales="free")+
    scale_color_manual(values=FST_COLORS)+
    labs(x="Model-standardized mean yield (Mg ha^-1)",
         y="Detrended modeled interannual SD (Mg ha^-1)")+
    theme(legend.position="none")
  fst_save_plot(f3,"Figure_2_mean_variability_fingerprint",12.5,4.6)

  # Figure 4. Lower-tail distribution of within-system common-year anomalies.
  anomalies<-fst_d("05_shared_year_yield_anomalies.csv")
  ns<-fst_d("04_downside_risk_summaries.csv")
  anomalies$Crop<-factor(anomalies$Crop,levels=FST_CROPS,labels=crops)
  anomalies$System<-factor(anomalies$System,levels=FST_SYSTEMS)
  ns$Crop<-factor(ns$Crop,levels=FST_CROPS,labels=crops)
  anomsystem<-paste(anomalies$Crop,anomalies$System)
  fst_assert(!anyDuplicated(paste(anomalies$Crop,anomalies$System,anomalies$Year)),
             "Duplicate common-year anomalies")
  f4<-ggplot(anomalies,aes(System,Anomaly_pct,color=System,fill=System))+
    geom_hline(yintercept=-20,linetype="dashed",color="#AD7657",linewidth=.7)+
    geom_boxplot(width=.48,alpha=.18,outlier.shape=NA,linewidth=.55)+
    geom_jitter(width=.15,height=0,size=1.8,alpha=.72)+
    facet_wrap(~Crop,nrow=1)+
    scale_color_manual(values=FST_COLORS)+scale_fill_manual(values=FST_COLORS)+
    labs(x=NULL,y="Yield anomaly vs. own common-year mean (%)")+
    theme(legend.position="none")
  fst_save_plot(f4,"Figure_4_common_year_downside_distributions",12,4.8)

  # Figure 5. Exploratory SPI-3--yield anomaly associations, no regression inference.
  wx<-fst_d("10_composite_annual_yields_with_weather.csv")
  wx$Crop<-factor(wx$Crop,levels=FST_CROPS,labels=crops)
  wx$System<-factor(wx$System,levels=FST_SYSTEMS)
  fst_assert(all(is.finite(wx$SPI3)),"SPI3 values missing in annual weather dataset")
  f5<-ggplot(wx,aes(SPI3,Anomaly_pct,color=System))+
    annotate("rect",xmin=-Inf,xmax=-1,ymin=-Inf,ymax=Inf,fill="#ECDAC8",alpha=.4)+
    annotate("rect",xmin=1,xmax=Inf,ymin=-Inf,ymax=Inf,fill="#E4EDF5",alpha=.45)+
    geom_hline(yintercept=0,color="#939CA2",linewidth=.45)+
    geom_vline(xintercept=c(-1,1),linetype="dotted",color="#99A4AD",linewidth=.4)+
    geom_point(size=1.8,alpha=.75)+
    facet_wrap(~Crop,nrow=1,scales="free_y")+
    scale_color_manual(values=FST_COLORS)+
    labs(x="Meteorological drought/wetness severity (SPI-3)",
         y="Observed yield anomaly vs. own system mean (%)")
  fst_save_plot(f5,"Figure_5_drought_yield_anomalies_exploratory",12.7,4.8)

  # Figure 6. Direct RT-FT paired-year yield difference; small conventional wheat omitted.
  till<-fst_d("06_tillage_paired_summaries.csv")
  till<-till[till$N_years>=4, ,drop=FALSE]
  till$Crop<-factor(till$Crop,levels=FST_CROPS,labels=crops)
  till$System<-factor(till$System,levels=rev(FST_SYSTEMS))
  till$Label<-paste0("n=",till$N_years)
  f6<-ggplot(till,aes(RT_minus_FT,System,color=System))+
    geom_vline(xintercept=0,linetype="dashed",color="#919CA5",linewidth=.55)+
    geom_segment(aes(x=CI_low,xend=CI_high,y=System,yend=System),linewidth=.85)+
    geom_point(size=2.8)+
    geom_text(aes(label=Label),nudge_y=.22,size=2.8,show.legend=FALSE)+
    facet_wrap(~Crop,nrow=1)+scale_color_manual(values=FST_COLORS,drop=FALSE)+
    labs(x="Reduced tillage minus full tillage yield (Mg ha^-1)",y=NULL)+
    theme(legend.position="none")
  fst_save_plot(f6,"Figure_6_paired_tillage_differences",11.8,4.3)

  # Figure S1. Supplemental difference-in-differences: organic minus CNV tillage response.
  did<-fst_d("07_system_x_tillage_interaction.csv")
  did<-did[did$N_shared_years>=4, ,drop=FALSE]
  did$Crop<-factor(did$Crop,levels=c("corn","soybean"),labels=c("Corn","Soybean"))
  did$System<-factor(did$System,levels=c("MNR","LEG"))
  did$Label<-paste0("n=",did$N_shared_years)
  f7<-ggplot(did,aes(Interaction_DiD,System,color=System))+
    geom_vline(xintercept=0,linetype="dashed",color="#919CA5",linewidth=.5)+
    geom_segment(aes(x=CI_low,xend=CI_high,y=System,yend=System),linewidth=.9)+
    geom_point(size=2.9)+
    geom_text(aes(label=Label),nudge_y=.17,size=2.8,show.legend=FALSE)+
    facet_wrap(~Crop,nrow=1)+scale_color_manual(values=FST_COLORS[c("LEG","MNR")])+
    labs(x="(RT-FT) organic minus (RT-FT) CNV (Mg ha^-1)",y=NULL)+
    theme(legend.position="none")
  fst_save_plot(f7,"Figure_S1_system_by_tillage_interactions",9.8,3.7)
  message("Created six main plot sets and one supplemental plot in Figures/",
          " (PDF, TIFF, PNG). Figure 2 uses fitted primary-model estimates.")
  invisible(TRUE)
}
