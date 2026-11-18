# FST grain yield and temporal stability analysis, 1981–2024.
# Pipeline entry point; analysis settings are defined in 01_SETTINGS.R.
# Annual yield records and daily weather are read from the parent data directory.
.v9_ofile <- vapply(sys.frames(),function(fr) {
   if(is.null(fr$ofile)) '' else as.character(fr$ofile)[1]
},character(1))
.v9_ofile <- .v9_ofile[nzchar(.v9_ofile)]
V9_CODE_DIR <- if(length(.v9_ofile))
  dirname(normalizePath(tail(.v9_ofile,1),mustWork=TRUE)) else getwd()
source(file.path(V9_CODE_DIR,'01_SETTINGS.R'),local=FALSE)
if(!RUN_MODE %in% c('SMOKE','FULL')) stop('RUN_MODE must be SMOKE or FULL')
required_packages<-c('readxl','ggplot2')
missing_packages<-required_packages[!vapply(required_packages,requireNamespace,logical(1),quietly=TRUE)]
if(length(missing_packages)) stop('Install packages before running: install.packages(c(',
    paste(sprintf('"%s"',missing_packages),collapse=', '),'))')
if(is.null(ROOT_OVERRIDE)||!nzchar(ROOT_OVERRIDE)) {
  if(file.exists(file.path(V9_CODE_DIR,YIELD_FILE))) V9_ROOT<-V9_CODE_DIR
  else V9_ROOT<-dirname(V9_CODE_DIR)
} else V9_ROOT<-ROOT_OVERRIDE
V9_ROOT<-normalizePath(V9_ROOT,mustWork=TRUE)
for(fn in c(YIELD_FILE,WEATHER_FILE)) {
  if(!file.exists(file.path(V9_ROOT,fn)))
    stop('Missing input file: ',file.path(V9_ROOT,fn),
         '\nKeep both original inputs in BG_DATA_UPDATE_V9, above the code folder.')
}
V9_OUT<-file.path(V9_ROOT,OUTPUT_FOLDER)
for(dir in c(V9_OUT,file.path(V9_OUT,'tables'),file.path(V9_OUT,'figures'),
             file.path(V9_OUT,'bootstrap')))
  dir.create(dir,showWarnings=FALSE,recursive=TRUE)
logfile<-file.path(V9_OUT,'00_PIPELINE_STATUS.txt')
writeLines(paste('RUNNING',Sys.time(),'Mode',RUN_MODE),logfile)
for(f in c('02_INPUT_AND_WEATHER.R','03_STATISTICAL_MODEL_CORE.R',
          '04_PRIMARY_AND_MODEL_CHECKS.R','05_BOOTSTRAP_AND_INFERENCE.R',
          '06_POST2008_TILLAGE.R','07_FIGURES_AND_MANUSCRIPT_TABLES.R'))
  source(file.path(V9_CODE_DIR,f),local=FALSE)
cat('\nV9 ROOT: ',V9_ROOT,'\nV9 OUTPUT: ',V9_OUT,'\nMODE: ',RUN_MODE,'\n',sep='')
tryCatch({
  cat('\n===== STAGE 1: VALIDATE ORIGINAL EXCEL / PRISM / ANNUAL COMPOSITE =====\n')
  .v9_yield<-v9_process_yield()
  .v9_weather<-v9_read_weather()
  .v9_data<-v9_join(.v9_yield,.v9_weather)
  v9_descriptives(.v9_data)
  cat('\n===== STAGE 2: FIT PRIMARY WEATHER-ADJUSTED TREND/COVARIANCE MODELS =====\n')
  .v9_models<-v9_primary_and_trends(.v9_data)
  if(RUN_SENSITIVITY){
    cat('\n===== STAGE 3: COMPARE TRENDS / B-UN / AR1 / WEATHER / WINDOWS =====\n')
    v9_compare_models(.v9_data,.v9_models)
    cat('\n===== STAGE 4: LEAVE-ONE-CALENDAR-YEAR-OUT DIAGNOSTICS =====\n')
    v9_influence(.v9_data,.v9_models)
  }
  if(RUN_BOOTSTRAP){
    cat('\n===== STAGE 5: MAIN CONDITIONAL BOOTSTRAP & ALTERNATIVE CORN UN =====\n')
    v9_run_selected_bootstrap(.v9_models)
    cat('\n===== STAGE 6: NULL PARAMETRIC-BOOTSTRAP OMNIBUS TESTS =====\n')
    v9_run_omnibus(.v9_models)
  }
  cat('\n===== STAGE 7: POST-2008 MATCHED FT VS RT TILLAGE =====\n')
  v9_run_tillage(.v9_data)
  cat('\n===== STAGE 8: PRODUCE MAIN TABLES, FIGURES, MANIFEST =====\n')
  v9_make_figures(.v9_data)
  v9_finish()
  writeLines(paste('COMPLETED',Sys.time(),'Mode',RUN_MODE),logfile)
},error=function(e){
  x<-paste('STOPPED',Sys.time(),'Mode',RUN_MODE,'ERROR:',conditionMessage(e))
  writeLines(x,logfile)
  message('\n',x,'\nSee ',logfile)
  stop(e)
})
