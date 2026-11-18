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
