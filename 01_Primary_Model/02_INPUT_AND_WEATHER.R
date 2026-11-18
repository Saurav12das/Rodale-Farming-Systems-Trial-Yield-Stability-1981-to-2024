# Plot-yield processing, annual treatment means, and seasonal weather indices.
v9_num <- function(x) suppressWarnings(as.numeric(gsub(',', '', trimws(as.character(x)))))
v9_key <- function(z) paste(z,collapse=' | ')
v9_save <- function(x,path) {
  dir.create(dirname(path),recursive=TRUE,showWarnings=FALSE)
  utils::write.csv(x,path,row.names=FALSE,na='')
  invisible(path)
}
v9_bind <- function(xs) {
  xs<-Filter(function(x) is.data.frame(x)&&nrow(x)>0,xs)
  if(length(xs)) do.call(rbind,xs) else data.frame()
}
v9_file <- function(name) file.path(V9_ROOT,name)
v9_output <- function(file,folder='tables') file.path(V9_OUT,folder,file)
v9_read_sheet <- function(sheet,crop) {
  x <- readxl::read_excel(v9_file(YIELD_FILE),sheet=sheet,col_types='text',.name_repair='minimal')
  if(ncol(x)<8L) stop('Yield sheet ',sheet,' has fewer than 8 columns')
  x<-x[,1:8,drop=FALSE]
  keep<-!is.na(x[[1]]) & nzchar(trimws(as.character(x[[1]])))
  x<-x[keep,,drop=FALSE]
  out<-data.frame(Crop=crop,Excel_row=which(keep)+1L,
    Year=v9_num(x[[1]]),Plot=v9_num(x[[2]]),Rep=v9_num(x[[3]]),
    EP=v9_num(x[[4]]),System=toupper(trimws(x[[5]])),
    Tillage=tolower(trimws(x[[7]])),Yield_kg_ha=v9_num(x[[8]]),
    stringsAsFactors=FALSE)
  out$Yield_Mg_ha<-out$Yield_kg_ha/1000
  out$QC_status<-ifelse(is.na(out$Yield_kg_ha),'EXCLUDED_MISSING_YIELD','VALID')
  out
}
v9_process_yield <- function() {
  sheets<-list(corn='corn_81_24',soybean='soy_81_24',wheat='wheat_81_24')
  all<-v9_bind(lapply(names(sheets),function(cr) v9_read_sheet(sheets[[cr]],cr)))
  bad <- !is.finite(all$Year)|!is.finite(all$Plot)|!is.finite(all$Rep)|
    !is.finite(all$EP)|!(all$System %in% SYSTEMS)|
    !(all$Tillage %in% c('till','no-till'))|
    (is.finite(all$Year)&!(all$Year %in% YEAR_RANGE))|
    (!is.na(all$Yield_kg_ha)&all$Yield_kg_ha<0)|
    (all$Year<ERA_BREAK&all$Tillage=='no-till')
  v9_save(all,v9_output('S1_01_plot_level_input_audit.csv'))
  if(any(bad)) {
    v9_save(all[bad,],v9_output('STOP_INVALID_PLOT_ROWS.csv'))
    stop('Invalid yield identifiers/values: see STOP_INVALID_PLOT_ROWS.csv')
  }
  valid<-all[all$QC_status=='VALID',,drop=FALSE]
  if(anyDuplicated(valid[c('Crop','Year','Plot')])) stop('Duplicate crop/year/physical plot rows')
  if(nrow(valid)!=2250L) warning('Current source has ',nrow(valid),' usable yields (historically 2250).')
  v9_save(valid,v9_output('S1_02_valid_plot_yields.csv'))
  groups<-split(valid,interaction(valid$Crop,valid$Year,valid$System,valid$Tillage,valid$Rep,drop=TRUE))
  rep_rows<-lapply(groups,function(z) data.frame(
    Crop=z$Crop[1],Year=z$Year[1],System=z$System[1],Tillage=z$Tillage[1],
    Rep=z$Rep[1],RepMean_Mg_ha=mean(z$Yield_Mg_ha),N_entry_points=length(unique(z$EP)),
    N_plot_rows=nrow(z)))
  reps<-v9_bind(rep_rows)
  v9_save(reps,v9_output('S1_03_replicate_means.csv'))
  agroups<-split(reps,interaction(reps$Crop,reps$Year,reps$System,reps$Tillage,drop=TRUE))
  arms<-v9_bind(lapply(agroups,function(z) data.frame(Crop=z$Crop[1],Year=z$Year[1],
    System=z$System[1],Tillage=z$Tillage[1],
    Yield_Mg_ha=mean(z$RepMean_Mg_ha),
    vi=if(nrow(z)>1) var(z$RepMean_Mg_ha)/nrow(z) else NA_real_,
    SE_Mg_ha=if(nrow(z)>1) sd(z$RepMean_Mg_ha)/sqrt(nrow(z)) else NA_real_,
    N_replicates=nrow(z))))
  if(any(!is.finite(arms$vi)|arms$vi<=0)) stop('Nonpositive/undefined annual variance: inspect S1_04.')
  v9_save(arms,v9_output('S1_04_arm_annual_means_and_sampling_variances.csv'))
  pre<-arms[arms$Year<ERA_BREAK&arms$Tillage=='till',,drop=FALSE]
  pre$Actual_arms<-'FT';pre$Post2008<-0L
  ft<-arms[arms$Year>=ERA_BREAK&arms$Tillage=='till',,drop=FALSE]
  rt<-arms[arms$Year>=ERA_BREAK&arms$Tillage=='no-till',,drop=FALSE]
  matched<-merge(ft,rt,by=c('Crop','Year','System'),suffixes=c('_FT','_RT'))
  post<-data.frame(Crop=matched$Crop,Year=matched$Year,System=matched$System,
    Yield_Mg_ha=(matched$Yield_Mg_ha_FT+matched$Yield_Mg_ha_RT)/2,
    vi=(matched$vi_FT+matched$vi_RT)/4,
    SE_Mg_ha=sqrt((matched$vi_FT+matched$vi_RT)/4),
    N_replicates=matched$N_replicates_FT+matched$N_replicates_RT,
    Actual_arms='FT+RT (50:50)',Post2008=1L)
  comp<-rbind(pre[,names(post)],post)
  comp<-comp[order(match(comp$Crop,CROPS),comp$Year,match(comp$System,SYSTEMS)),,drop=FALSE]
  rownames(comp)<-NULL
  if(anyDuplicated(comp[c('Crop','Year','System')])) stop('Duplicate annual composite keys')
  if(nrow(comp)!=241L) warning('Expected 241 annual composite observations; got ',nrow(comp))
  v9_save(comp,v9_output('S1_05_composite_annual_means_PRIMARY.csv'))
  v9_save(matched,v9_output('S1_06_post2008_paired_FT_RT.csv'))
  # The composite sampling variance assumes independent FT and RT replicate groups.
  # Replicate identifiers are checked for overlap after the tillage split.
  overlap_rows<-lapply(seq_len(nrow(matched)),function(i){
    z<-matched[i,,drop=FALSE]
    rp<-reps[reps$Crop==z$Crop & reps$Year==z$Year &
      reps$System==z$System,,drop=FALSE]
    ftids<-rp$Rep[rp$Tillage=='till'];rtids<-rp$Rep[rp$Tillage=='no-till']
    data.frame(Crop=z$Crop,Year=z$Year,System=z$System,
      N_FT_replicates=length(ftids),N_RT_replicates=length(rtids),
      Overlapping_replicate_IDs=paste(intersect(ftids,rtids),collapse=';'),
      N_overlapping_replicate_IDs=length(intersect(ftids,rtids)))
  })
  overlap<-v9_bind(overlap_rows)
  v9_save(overlap,v9_output('S1_06b_FT_RT_replicate_overlap_audit.csv'))
  if(nrow(overlap)&&any(overlap$N_overlapping_replicate_IDs>0))
    warning('Some FT and RT replicate IDs overlap; revisit composite vi independence assumption.')
  unpaired<-arms[arms$Year>=ERA_BREAK &
    !(paste(arms$Crop,arms$Year,arms$System) %in%
      paste(matched$Crop,matched$Year,matched$System)),,drop=FALSE]
  v9_save(unpaired,v9_output('S1_07_unpaired_arms_excluded_from_composite.csv'))
  valid$Tillage_group<-ifelse(valid$Tillage=='till','FT','RT')
  coverage<-as.data.frame(table(valid$Crop,valid$System,valid$Tillage_group))
  names(coverage)<-c('Crop','System','Tillage','N_valid_plot_yields')
  v9_save(coverage,v9_output('S1_08_plot_coverage.csv'))
  list(plots=valid,replicate=reps,arms=arms,matched=matched,composite=comp,unpaired=unpaired)
}
v9_read_weather <- function() {
  fp<-v9_file(WEATHER_FILE)
  hdr<-readLines(fp,n=35,warn=FALSE)
  at<-grep('^[[:space:]]*"?Date"?[[:space:]]*,',sub('^\ufeff','',hdr))
  if(!length(at)) stop('No PRISM Date header found in first 35 lines')
  raw<-read.csv(fp,skip=at[1]-1,check.names=FALSE,stringsAsFactors=FALSE)
  req<-c('Date','ppt (inches)','tmin (degrees F)','tmean (degrees F)','tmax (degrees F)')
  if(!all(req%in%names(raw))) stop('PRISM header does not contain expected Fahrenheit/inches columns')
  dd<-as.Date(raw[['Date']],format='%m/%d/%Y')
  if(anyNA(dd)) stop('PRISM date parsing failed; expected month/day/year')
  wx<-data.frame(Date=dd,ppt_mm=v9_num(raw[['ppt (inches)']])*25.4,
    tmin_C=(v9_num(raw[['tmin (degrees F)']])-32)*5/9,
    tmean_C=(v9_num(raw[['tmean (degrees F)']])-32)*5/9,
    tmax_C=(v9_num(raw[['tmax (degrees F)']])-32)*5/9)
  wx<-wx[wx$Date>=as.Date('1981-01-01')&wx$Date<=as.Date('2024-12-31'),,drop=FALSE]
  wx<-wx[order(wx$Date),,drop=FALSE]
  expected<-seq.Date(as.Date('1981-01-01'),as.Date('2024-12-31'),by='day')
  if(nrow(wx)!=length(expected)||anyDuplicated(as.integer(wx$Date))||
     !identical(as.integer(wx$Date),as.integer(expected)))
     stop('PRISM contains missing, repeated or unexpected daily dates during 1981-2024')
  if(any(!is.finite(as.matrix(wx[,2:5])))||any(wx$ppt_mm<0)||
     any(wx$tmax_C<wx$tmin_C)) stop('Invalid daily PRISM weather variables')
  v9_save(wx,v9_output('S1_09_PRISM_daily_SI_1981_2024.csv'))
  season<-function(a,b) {
    a<-as.Date(a);b<-as.Date(b)
    idx<-which(wx$Date>=a&wx$Date<=b)
    dd<-seq.Date(a,b,by='day')
    if(length(idx)!=length(dd)||!all(as.integer(wx$Date[idx])==as.integer(dd))) return(NULL)
    wx[idx,,drop=FALSE]
  }
  rows<-list()
  for(cr in CROPS) for(yr in YEAR_RANGE) {
    a<-season(sprintf('%d-05-01',yr),sprintf('%d-09-30',yr))
    j<-season(sprintf('%d-07-01',yr),sprintf('%d-08-31',yr))
    k<-season(sprintf('%d-06-01',yr),sprintf('%d-08-31',yr))
    q<-season(sprintf('%d-04-01',yr),sprintf('%d-06-30',yr))
    if(cr=='wheat'){
      p<-season(sprintf('%d-10-01',yr-1),sprintf('%d-07-31',yr))
      complete<-!is.null(p)&&!is.null(q)
      r1<-if(complete) sum(p$ppt_mm) else NA_real_
      h1<-if(complete) mean(q$tmean_C) else NA_real_
      r2<-if(complete) sum(q$ppt_mm) else NA_real_
      h2<-h1
    } else {
      complete<-!is.null(a)&&!is.null(j)&&!is.null(k)
      r1<-if(complete) sum(a$ppt_mm) else NA_real_
      h1<-if(complete) sum(j$tmax_C>=HOT_DAY_C) else NA_real_
      r2<-if(complete) sum(j$ppt_mm) else NA_real_
      h2<-if(complete) sum(k$tmax_C>=HOT_DAY_C) else NA_real_
    }
    rows[[length(rows)+1L]]<-data.frame(Crop=cr,Year=yr,Post2008=as.integer(yr>=ERA_BREAK),
      WeatherWindowComplete=complete,Rain_primary=r1,Heat_primary=h1,
      Rain_alternative=r2,Heat_alternative=h2)
  }
  yearly<-v9_bind(rows)
  v9_save(yearly,v9_output('S1_10_crop_year_weather_indices.csv'))
  list(daily=wx,yearly=yearly)
}
v9_join <- function(y,weather){
  annual<-merge(y$composite,weather$yearly,by=c('Crop','Year','Post2008'),all.x=TRUE,sort=FALSE)
  annual<-annual[order(match(annual$Crop,CROPS),annual$Year,match(annual$System,SYSTEMS)),,drop=FALSE]
  if(anyNA(annual$WeatherWindowComplete)||any(!annual$WeatherWindowComplete))
     stop('Incomplete PRISM windows for annual yield observations: inspect data')
  if(anyDuplicated(annual[c('Crop','Year','System')])) stop('Duplicate merged year/system/crop')
  v9_save(annual,v9_output('S1_11_annual_yield_weather_PRIMARY.csv'))
  cover<-v9_bind(lapply(split(annual,interaction(annual$Crop,annual$System,drop=TRUE)),function(z)
    data.frame(Crop=z$Crop[1],System=z$System[1],N_years=nrow(z),
      First_year=min(z$Year),Last_year=max(z$Year),
      N_pre2008=sum(z$Year<ERA_BREAK),N_post2008=sum(z$Year>=ERA_BREAK))))
  v9_save(cover,v9_output('S1_12_system_crop_year_coverage.csv'))
  list(annual=annual,weather=weather$yearly,arms=y$arms,matched=y$matched,plots=y$plots)
}
