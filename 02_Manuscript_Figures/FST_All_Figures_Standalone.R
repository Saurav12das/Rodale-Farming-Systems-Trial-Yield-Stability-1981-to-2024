# FST yield-stability figures: standalone R script (no source() calls needed).
# Required: R >= 4.2, packages readxl and ggplot2.
# Put this script with the raw Excel yield workbook and PRISM daily CSV,
# OR one folder below those two raw files.
# Run: Rscript FST_All_Figures_Standalone.R
# Figures appear under Figures/ as PDF, 600-dpi TIFF, and PNG.
# Figure 2 uses archived primary-model point estimates and bootstrap CIs,
# unless you have re-fitted the primary model and put its outputs in the expected directory.
# When only cached figure-input tables exist, --plots-only redraws figures.
arg <- grep("^--file=", commandArgs(FALSE), value=TRUE)
CODE_DIR <- if(length(arg)) dirname(normalizePath(sub("^--file=", "", arg[1]))) else getwd()
OUTPUT_DIR <- CODE_DIR
candidates <- c(CODE_DIR, dirname(CODE_DIR))
has_raw <- function(p) {
   any(file.exists(file.path(p, c("FST_1981_2024_3_crops_BG(1).xlsx", "FST_1981_2024_3_crops_BG.xlsx")))) &&
   length(list.files(p, pattern="^PRISM_ppt_tmin_tmean_tmax_tdmean_soltotal_.*\\.csv$")) == 1
}
okdirs <- candidates[vapply(candidates, has_raw, logical(1))]
ROOT_DIR <- if(length(okdirs)) okdirs[1] else dirname(CODE_DIR)

# Reproducible figure-data construction for Rodale FST, 1981--2024
# Supplement to the yield-stability manuscript. R >= 4.2; packages: readxl.
# Source: plot-level yield Excel and PRISM 800-m daily weather CSV.
# This script uses replicate-first weighting, explicitly excludes unpaired
# post-2008 arms from the PRIMARY system series, and never turns a nonharvest
# year into zero yield. All comparisons below are DESCRIPTIVE / exploratory.

fst_num <- function(x) suppressWarnings(as.numeric(gsub(",", "", trimws(as.character(x)), fixed=TRUE)))
fst_write <- function(x, filename) {
  dir.create(file.path(OUTPUT_DIR, "Derived_data"), recursive=TRUE, showWarnings=FALSE)
  utils::write.csv(x, file.path(OUTPUT_DIR, "Derived_data", filename), row.names=FALSE, na="")
}
fst_assert <- function(cond, msg) if(!isTRUE(cond)) stop(msg, call.=FALSE)
SYSTEMS <- c("CNV", "LEG", "MNR")
CROPS <- c("corn", "soybean", "wheat")

fst_load_yields <- function(path) {
  if(!requireNamespace("readxl", quietly=TRUE)) stop("Package 'readxl' is required")
  sheet_names <- c(corn="corn_81_24", soybean="soy_81_24", wheat="wheat_81_24")
  pieces <- lapply(names(sheet_names), function(crop) {
    x <- readxl::read_excel(path, sheet=sheet_names[[crop]], col_types="text", .name_repair="minimal")
    fst_assert(ncol(x)>=8, paste("Bad Excel worksheet", crop))
    x <- x[, 1:8, drop=FALSE]
    x <- x[!is.na(x[[1]]) & nzchar(trimws(as.character(x[[1]]))), ,drop=FALSE]
    data.frame(Crop=crop, Year=fst_num(x[[1]]), Plot=fst_num(x[[2]]),
      Rep=fst_num(x[[3]]), EP=fst_num(x[[4]]),
      System=toupper(trimws(as.character(x[[5]]))),
      Tillage=tolower(trimws(as.character(x[[7]]))),
      Yield_Mg_ha=fst_num(x[[8]])/1000, stringsAsFactors=FALSE)
  })
  dat <- do.call(rbind, pieces); rownames(dat) <- NULL
  fst_assert(nrow(dat)==2251, sprintf("Expected 2,251 source rows, got %s", nrow(dat)))
  invalid <- !is.finite(dat$Year) | !is.finite(dat$Plot) | !is.finite(dat$Rep) |
    !is.finite(dat$EP) | !(dat$System %in% SYSTEMS) |
    !(dat$Tillage %in% c("till", "no-till")) |
    !(dat$Year %in% 1981:2024) | (!is.na(dat$Yield_Mg_ha) & dat$Yield_Mg_ha<0)
  fst_assert(!any(invalid), "Invalid source identifiers, negative yields, or out-of-range years")
  valid <- dat[is.finite(dat$Yield_Mg_ha), ,drop=FALSE]
  fst_assert(nrow(valid)==2250, "Expected 2,250 valid plot yields")
  key <- paste(valid$Crop, valid$Year, valid$Plot, sep="|")
  fst_assert(!anyDuplicated(key), "Duplicate Crop-Year-Plot identifier")
  reps <- stats::aggregate(Yield_Mg_ha ~ Crop+Year+System+Tillage+Rep,
                            data=valid, FUN=mean)
  names(reps)[names(reps)=="Yield_Mg_ha"] <- "RepMean"
  gs <- split(reps, interaction(reps$Crop, reps$Year, reps$System,
                                 reps$Tillage, drop=TRUE))
  arms <- do.call(rbind, lapply(gs, function(z) {
    fst_assert(nrow(z)>=2, paste("Insufficient replication", z$Crop[1],z$Year[1],z$System[1]))
    data.frame(Crop=z$Crop[1], Year=z$Year[1], System=z$System[1],
               Tillage=z$Tillage[1], Yield_Mg_ha=mean(z$RepMean),
               vi=stats::var(z$RepMean)/nrow(z), Nrep=nrow(z))
  }))
  rownames(arms)<-NULL
  fst_assert(all(is.finite(arms$vi)) && all(arms$vi>0), "Nonpositive annual mean sampling variance")
  pre <- arms[arms$Year<2008 & arms$Tillage=="till", ,drop=FALSE]
  pre_comp <- data.frame(Crop=pre$Crop, Year=pre$Year, System=pre$System,
              Yield_Mg_ha=pre$Yield_Mg_ha, vi=pre$vi, Inclusion="FT pre2008")
  ft <- arms[arms$Year>=2008 & arms$Tillage=="till", ,drop=FALSE]
  rt <- arms[arms$Year>=2008 & arms$Tillage=="no-till", ,drop=FALSE]
  matched <- merge(ft,rt,by=c("Crop","Year","System"), suffixes=c("_FT","_RT"))
  post_comp <- data.frame(Crop=matched$Crop, Year=matched$Year, System=matched$System,
    Yield_Mg_ha=(matched$Yield_Mg_ha_FT+matched$Yield_Mg_ha_RT)/2,
    vi=(matched$vi_FT+matched$vi_RT)/4, Inclusion="paired FT+RT post2008")
  comp <- rbind(pre_comp, post_comp)
  comp <- comp[order(match(comp$Crop,CROPS),comp$Year,match(comp$System,SYSTEMS)),]
  rownames(comp)<-NULL
  fst_assert(nrow(comp)==241, paste("Expected 241 annual composite means; got",nrow(comp)))
  fst_assert(!anyDuplicated(paste(comp$Crop,comp$Year,comp$System)), "Duplicate composite year")
  list(raw=dat, valid=valid, reps=reps, arms=arms, matched=matched, composite=comp)
}

fst_pair <- function(a,b) { # a minus b; same ordered calendar years
  n <- length(a); fst_assert(length(b)==n && n>0,"Unmatched or empty pair")
  dd <- a-b; dm <- mean(dd); se <- if(n>1) stats::sd(dd)/sqrt(n) else NA_real_
  ci <- if(n>1) dm+c(-1,1)*stats::qt(0.975,n-1)*se else c(NA_real_,NA_real_)
  c(N=n, Difference=dm, CI_low=ci[1], CI_high=ci[2],
    p_unadjusted=if(n>=3 && stats::sd(dd)>0) 2*stats::pt(-abs(dm/se),n-1) else NA_real_,
    A_mean=mean(a), B_mean=mean(b))
}

fst_coverage <- function(y) {
  grid <- expand.grid(Crop=CROPS, System=SYSTEMS, Year=1981:2024,
                       stringsAsFactors=FALSE)
  included <- paste(y$composite$Crop,y$composite$System,y$composite$Year)
  observed <- paste(y$arms$Crop,y$arms$System,y$arms$Year)
  key <- paste(grid$Crop,grid$System,grid$Year)
  grid$Status <- ifelse(key %in% included,
                        ifelse(grid$Year<2008,"Included: FT","Included: FT+RT"),
                        ifelse(key %in% observed,"One arm only: excluded","No harvest"))
  grid
}

fst_gaps <- function(y) {
  x <- y$composite
  long <- list(); rows <- list(); cnt<-0L
  for(cr in CROPS) for(s in c("LEG","MNR")) {
    a <- x[x$Crop==cr & x$System==s, c("Year","Yield_Mg_ha")]
    b <- x[x$Crop==cr & x$System=="CNV", c("Year","Yield_Mg_ha")]
    z <- merge(a,b,by="Year",suffixes=c("_ORG","_CNV"));z<-z[order(z$Year),]
    if(!nrow(z)) next
    p <- fst_pair(z$Yield_Mg_ha_ORG,z$Yield_Mg_ha_CNV)
    cnt<-cnt+1L
    rows[[cnt]]<-data.frame(Crop=cr,System=s,N_years=nrow(z),
       Mean_CNV=p["B_mean"],Mean_org=p["A_mean"],
       Organic_minus_CNV=p["Difference"],CI_low=p["CI_low"],CI_high=p["CI_high"],
       p_unadjusted=p["p_unadjusted"],Gap_pct=100*(p["B_mean"]-p["A_mean"])/p["B_mean"])
    long[[cnt]]<-data.frame(Crop=cr,System=s,Year=z$Year,
                             Gap_CNV_minus_organic=z$Yield_Mg_ha_CNV-z$Yield_Mg_ha_ORG)
  }
  list(summary=do.call(rbind,rows),yearly=do.call(rbind,long))
}

fst_tillage <- function(y) {
  m<-y$matched; a<-y$arms
  out<-list(); dd<-list(); ix<-list(); by_till<-list();n<-0L
  for(cr in CROPS) for(s in SYSTEMS) {
    z<-m[m$Crop==cr & m$System==s, ,drop=FALSE]
    if(!nrow(z)) next
    z<-z[order(z$Year),]
    p<-fst_pair(z$Yield_Mg_ha_RT,z$Yield_Mg_ha_FT)
    n<-n+1L
    out[[n]]<-data.frame(Crop=cr,System=s,N_years=nrow(z),RT_minus_FT=p["Difference"],
           CI_low=p["CI_low"],CI_high=p["CI_high"],p_unadjusted=p["p_unadjusted"],
           FT_mean=p["B_mean"],RT_mean=p["A_mean"],
           RT_penalty_pct=100*(p["B_mean"]-p["A_mean"])/p["B_mean"])
  }
  for(cr in CROPS) for(s in c("LEG","MNR")) {
    m1<-m[m$Crop==cr & m$System==s, ,drop=FALSE]
    m0<-m[m$Crop==cr & m$System=="CNV", ,drop=FALSE]
    z<-merge(m1,m0,by="Year",suffixes=c("_ORG","_CNV"))
    if(!nrow(z)) next
    org <- z$Yield_Mg_ha_RT_ORG-z$Yield_Mg_ha_FT_ORG
    cnv <- z$Yield_Mg_ha_RT_CNV-z$Yield_Mg_ha_FT_CNV
    p <- fst_pair(org,cnv)
    dd[[length(dd)+1L]]<-data.frame(Crop=cr,System=s,N_shared_years=nrow(z),
      Interaction_DiD=p["Difference"],CI_low=p["CI_low"],CI_high=p["CI_high"],
      p_unadjusted=p["p_unadjusted"])
    for(tl in c("till","no-till")) {
      za<-a[a$Crop==cr&a$System==s&a$Tillage==tl&a$Year>=2008, c("Year","Yield_Mg_ha")]
      zb<-a[a$Crop==cr&a$System=="CNV"&a$Tillage==tl&a$Year>=2008,c("Year","Yield_Mg_ha")]
      zz<-merge(za,zb,by="Year",suffixes=c("_ORG","_CNV"))
      if(!nrow(zz)) next
      by_till[[length(by_till)+1L]]<-data.frame(Crop=cr,System=s,Tillage=if(tl=="till") "FT" else "RT",
        N_years=nrow(zz),Gap_percent=100*(mean(zz$Yield_Mg_ha_CNV)-mean(zz$Yield_Mg_ha_ORG))/mean(zz$Yield_Mg_ha_CNV))
    }
  }
  list(summary=do.call(rbind,out),interaction=do.call(rbind,dd),by_tillage=do.call(rbind,by_till))
}

fst_downside <- function(y) {
  x<-y$composite; an<-list();summ<-list()
  for(cr in CROPS) {
    years<-sort(unique(x$Year[x$Crop==cr]))
    years<-years[vapply(years,function(yr) all(SYSTEMS %in% x$System[x$Crop==cr & x$Year==yr]), logical(1))]
    for(s in SYSTEMS) {
      z<-x[x$Crop==cr & x$System==s & x$Year %in% years,,drop=FALSE]
      z<-z[order(z$Year),]; fst_assert(nrow(z)==length(years),"Downside years misaligned")
      ave<-mean(z$Yield_Mg_ha)
      vals<-100*(z$Yield_Mg_ha/ave-1)
      an[[length(an)+1L]]<-data.frame(Crop=cr,System=s,Year=z$Year,Anomaly_pct=vals)
      summ[[length(summ)+1L]]<-data.frame(Crop=cr,System=s,N_common=nrow(z),
        Mean_Mg_ha=ave,SD_Mg_ha=sd(z$Yield_Mg_ha), CV_pct=100*sd(z$Yield_Mg_ha)/ave,
        N_below_80pct_mean=sum(z$Yield_Mg_ha < .8*ave))
    }
  }
  list(anomalies=do.call(rbind,an),summary=do.call(rbind,summ))
}

fst_read_prism <- function(path) {
  raw <- utils::read.csv(path, check.names=FALSE,stringsAsFactors=FALSE)
  req <- c("Date","ppt (inches)","tmin (degrees F)","tmean (degrees F)","tmax (degrees F)")
  fst_assert(all(req %in% names(raw)), "PRISM headers unexpected")
  dates <- as.Date(raw$Date, format="%m/%d/%Y")
  keep<-dates>=as.Date("1981-01-01") & dates<=as.Date("2024-12-31")
  w<-data.frame(Date=dates[keep], P_mm=fst_num(raw[["ppt (inches)"]][keep])*25.4,
    Tmin=(fst_num(raw[["tmin (degrees F)"]][keep])-32)*5/9,
    Tmean=(fst_num(raw[["tmean (degrees F)"]][keep])-32)*5/9,
    Tmax=(fst_num(raw[["tmax (degrees F)"]][keep])-32)*5/9)
  w<-w[order(w$Date),]
  expect<-seq(as.Date("1981-01-01"),as.Date("2024-12-31"),by="day")
  fst_assert(identical(as.integer(w$Date),as.integer(expect)),"Missing or repeated PRISM daily dates")
  fst_assert(all(is.finite(w$P_mm)) && all(is.finite(w$Tmean)) &&
             all(is.finite(w$Tmax)) && all(is.finite(w$Tmin)),"Missing PRISM numeric values")
  w$Year<-as.integer(format(w$Date,"%Y"));w$Month<-as.integer(format(w$Date,"%m"))
  # Hargreaves-Samani temperature-based reference evapotranspiration proxy.
  latitude <- 40.5509*pi/180
  doy <- as.integer(format(w$Date,"%j"))
  dr <- 1 + 0.033*cos(2*pi*doy/365)
  delta <- 0.409*sin(2*pi*doy/365-1.39)
  omega <- acos(pmax(-1,pmin(1,-tan(latitude)*tan(delta))))
  Ra <- (24*60/pi)*0.082*dr*(omega*sin(latitude)*sin(delta)+cos(latitude)*cos(delta)*sin(omega))
  w$ET0_HS_mm <- pmax(0,0.0023*(w$Tmean+17.8)*sqrt(pmax(w$Tmax-w$Tmin,0))*Ra/2.45)
  w
}

fst_spi3 <- function(w) {
  mm <- stats::aggregate(P_mm ~ Year+Month, data=w, FUN=sum)
  mm<-mm[order(mm$Year,mm$Month),]
  mm$P3 <- NA_real_
  for(i in 3:nrow(mm)) mm$P3[i]<-sum(mm$P_mm[(i-2):i])
  mm$SPI3<-NA_real_
  for(mo in 1:12) {
    ii <- which(mm$Month==mo & mm$Year>=1991 & mm$Year<=2020)
    obs <- mm$P3[ii];fst_assert(length(obs)==30 && all(is.finite(obs)),"SPI reference incomplete")
    p0<-mean(obs<=0);positive<-obs[obs>0]
    fst_assert(length(positive)>=10,"Insufficient positive reference precipitation")
    mean0<-mean(positive);var0<-stats::var(positive)
    initial<-log(c(max(.01,mean0^2/var0), max(.01,var0/mean0)))
    negll<-function(u) -sum(stats::dgamma(positive,shape=exp(u[1]),scale=exp(u[2]),log=TRUE))
    fit<-stats::optim(initial,negll,method="BFGS",control=list(maxit=1000))
    fst_assert(fit$convergence==0,paste("SPI gamma MLE did not converge month",mo))
    shape<-exp(fit$par[1]);scale<-exp(fit$par[2])
    jj<-which(mm$Month==mo & is.finite(mm$P3))
    xx<-mm$P3[jj]
    cdf<-ifelse(xx>0,p0+(1-p0)*stats::pgamma(xx,shape=shape,scale=scale),p0/2)
    mm$SPI3[jj]<-stats::qnorm(pmin(1-1e-5,pmax(1e-5,cdf)))
  }
  mm
}

fst_crop_weather <- function(w,mm) {
  out<-list();ix<-0L
  for(cr in CROPS) for(yr in 1981:2024) {
    if(cr=="wheat") {start<-as.Date(sprintf("%s-04-01",yr));end<-as.Date(sprintf("%s-06-30",yr));spi_m<-5}
    else {start<-as.Date(sprintf("%s-05-01",yr));end<-as.Date(sprintf("%s-08-31",yr));spi_m<-8}
    z<-w[w$Date>=start&w$Date<=end,,drop=FALSE]
    fst_assert(nrow(z)>0,"Missing seasonal weather")
    k<-which(mm$Year==yr & mm$Month==spi_m)
    d<-as.integer(z$P_mm<1)
    rr<-rle(d);max_dry<-if(any(rr$values==1)) max(rr$lengths[rr$values==1]) else 0L
    hot<-z$Tmax>=32
    if(cr!="wheat") hot<-hot & (z$Month %in% c(7,8))
    ix<-ix+1L
    out[[ix]]<-data.frame(Crop=cr,Year=yr,SPI3=mm$SPI3[k],
      Season_P_mm=sum(z$P_mm),Season_ET0_HS_mm=sum(z$ET0_HS_mm),
      P_minus_ET0_mm=sum(z$P_mm-z$ET0_HS_mm),MaxDrySpell_days=max_dry,
      HeavyRain20_days=sum(z$P_mm>=20),HotDays32C=sum(hot))
  }
  do.call(rbind,out)
}


# Use a primary-model output table if available. Otherwise keep the bundled
# previously exported model estimates. This does not re-estimate the REML model.
fst_model_export <- function() {
  candidates <- c(
    file.path(OUTPUT_DIR,"FST_V9_FINAL_ANALYSIS_outputs","tables","MAIN_01_primary_point_estimates_ALL_definitions.csv"),
    file.path(ROOT_DIR,"FST_V9_FINAL_ANALYSIS_outputs","tables","MAIN_01_primary_point_estimates_ALL_definitions.csv"),
    file.path(ROOT_DIR,"Model_Analysis_V9","FST_V9_FINAL_ANALYSIS_outputs","tables","MAIN_01_primary_point_estimates_ALL_definitions.csv")
  )
  candidates <- candidates[file.exists(candidates)]
  if(length(candidates)) {
    m <- read.csv(candidates[1], stringsAsFactors=FALSE)
    cols <- c("Crop","System","Mean_Mg_ha","Detrended_total_SD_Mg_ha")
    fst_assert(all(cols %in% names(m)),"Primary model output missing expected columns")
    out <- data.frame(Crop=m$Crop,System=m$System,Mean_Mg_ha=m$Mean_Mg_ha,
      SD_Mg_ha=m$Detrended_total_SD_Mg_ha,
      Mean_CI_low=NA_real_,Mean_CI_high=NA_real_,
      SD_CI_low=NA_real_,SD_CI_high=NA_real_,
      Source="Primary trend-weather REML model; refitted in R")
    # Where bootstrap CI outputs are available, join them by metric.
    cibase <- file.path(dirname(candidates[1]),"MAIN_04_bootstrap_mean_SD_CV_95CI.csv")
    if(file.exists(cibase)) {
      ci <- read.csv(cibase,stringsAsFactors=FALSE)
      for(metric in c("Mean_Mg_ha","Detrended_total_SD_Mg_ha")) {
        z <- ci[ci$Metric==metric,,drop=FALSE]
        j <- match(paste(out$Crop,out$System),paste(z$Crop,z$System))
        prefix <- if(metric=="Mean_Mg_ha") "Mean" else "SD"
        out[[paste0(prefix,"_CI_low")]] <- z$CI_low[j]
        out[[paste0(prefix,"_CI_high")]] <- z$CI_high[j]
      }
    }
    fst_write(out,"12_primary_model_point_estimates.csv")
    message("Updated Figure 2 point estimates from the primary REML output")
  } else {
    fn <- file.path(OUTPUT_DIR,"Derived_data","12_primary_model_point_estimates.csv")
    fst_assert(file.exists(fn),"Need primary model point estimates for Figure 2")
    message("Using archived primary-model point estimates for Figure 2; rerun model pipeline to refresh")
  }
  invisible(TRUE)
}

# main function: ROOT_DIR + OUTPUT_DIR set in 00_RUN_ALL_FIGURES.R
fst_prepare <- function() {
  yf <- c(file.path(ROOT_DIR,"FST_1981_2024_3_crops_BG(1).xlsx"),
          file.path(ROOT_DIR,"FST_1981_2024_3_crops_BG.xlsx"))
  yf <- yf[file.exists(yf)]
  fst_assert(length(yf)>0,"Original yield Excel not found under ROOT_DIR")
  wxfiles <- list.files(ROOT_DIR,pattern="^PRISM_ppt_tmin_tmean_tmax_tdmean_soltotal_.*\\.csv$",full.names=TRUE)
  fst_assert(length(wxfiles)==1,"Expected one original PRISM daily CSV in ROOT_DIR")
  y<-fst_load_yields(yf[1])
  gaps<-fst_gaps(y)
  till<-fst_tillage(y)
  down<-fst_downside(y)
  w<-fst_read_prism(wxfiles[1]);spi<-fst_spi3(w)
  cropweather<-fst_crop_weather(w,spi)
  wx_composite <- merge(y$composite, cropweather, by=c("Crop","Year"),all.x=TRUE)
  fst_assert(nrow(wx_composite)==241 && all(is.finite(wx_composite$SPI3)),"Weather merge failed")
  wx_composite$Anomaly_pct<-NA_real_
  for(cr in CROPS) for(s in SYSTEMS) {
    i<-which(wx_composite$Crop==cr & wx_composite$System==s)
    wx_composite$Anomaly_pct[i]<-100*(wx_composite$Yield_Mg_ha[i]/mean(wx_composite$Yield_Mg_ha[i])-1)
  }
  cov<-fst_coverage(y)
  fst_write(cov, "01_harvest_year_coverage.csv")
  fst_write(gaps$summary,"02_system_gap_summaries.csv")
  fst_write(gaps$yearly,"03_system_gap_by_year.csv")
  fst_write(down$summary,"04_downside_risk_summaries.csv")
  fst_write(down$anomalies,"05_shared_year_yield_anomalies.csv")
  fst_write(till$summary,"06_tillage_paired_summaries.csv")
  fst_write(till$interaction,"07_system_x_tillage_interaction.csv")
  fst_write(till$by_tillage,"08_system_gaps_by_tillage.csv")
  fst_write(cropweather,"09_PRISM_crop_year_indices.csv")
  fst_write(wx_composite,"10_composite_annual_yields_with_weather.csv")
  fst_write(y$composite,"11_primary_annual_composites.csv")
  fst_model_export()
  cat("Prepared data: ",nrow(y$valid)," plots; ",nrow(y$composite)," crop-system years; ",nrow(w)," PRISM days.\n",sep="")
  invisible(TRUE)
}

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


FST_CACHED_MODEL_CSV <- "Crop,System,Mean_Mg_ha,SD_Mg_ha,Source,Mean_CI_low,Mean_CI_high,SD_CI_low,SD_CI_high\ncorn,CNV,7.65,2.72,Primary REML model; 95% conditional bootstrap intervals,6.33,8.92,2.15,3.43\ncorn,LEG,6.79,2.37,Primary REML model; 95% conditional bootstrap intervals,5.7,7.89,1.86,3.02\ncorn,MNR,7.5,2.24,Primary REML model; 95% conditional bootstrap intervals,6.45,8.53,1.72,2.91\nsoybean,CNV,3.47,0.61,Primary REML model; 95% conditional bootstrap intervals,3.26,3.68,0.47,0.77\nsoybean,LEG,2.49,0.57,Primary REML model; 95% conditional bootstrap intervals,2.29,2.69,0.44,0.74\nsoybean,MNR,2.85,0.61,Primary REML model; 95% conditional bootstrap intervals,2.6,3.09,0.45,0.79\nwheat,CNV,4.08,0.84,Primary REML model; 95% conditional bootstrap intervals,3.6,4.56,0.63,1.09\nwheat,LEG,3.02,0.79,Primary REML model; 95% conditional bootstrap intervals,2.78,3.26,0.6,1.03\nwheat,MNR,3.43,0.9,Primary REML model; 95% conditional bootstrap intervals,3.13,3.73,0.71,1.13"

# Cached estimates from the primary Gaussian REML analysis for Figure 2.
# These are NOT re-estimated by the plotting script. Derived descriptive
# tables are recalculated from the raw workbook and daily weather file.
if(!file.exists(file.path(OUTPUT_DIR, "Derived_data", "12_primary_model_point_estimates.csv"))){
  dir.create(file.path(OUTPUT_DIR,"Derived_data"),showWarnings=FALSE,recursive=TRUE)
  .model_cached <- read.csv(text=FST_CACHED_MODEL_CSV,stringsAsFactors=FALSE)
  utils::write.csv(.model_cached,
    file.path(OUTPUT_DIR,"Derived_data","12_primary_model_point_estimates.csv"),row.names=FALSE)
}
if(!("--plots-only" %in% commandArgs(TRUE))) fst_prepare()
fst_draw_figures()
