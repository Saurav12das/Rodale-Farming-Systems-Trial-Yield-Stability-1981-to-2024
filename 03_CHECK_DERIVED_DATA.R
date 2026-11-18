# Minimal figure-data integrity checks; run from supplementary directory.
root <- getwd(); d <- file.path(root,"Derived_data")
read <- function(x) read.csv(file.path(d,x),stringsAsFactors=FALSE)
a <- read("11_primary_annual_composites.csv")
stopifnot(nrow(a)==241L,all(c("corn","soybean","wheat") %in% a$Crop))
stopifnot(identical(as.integer(table(a$Crop)[c("corn","soybean","wheat")]),c(99L,89L,53L)))
g <- read("02_system_gap_summaries.csv")
leg <- g[g$Crop=="corn" & g$System=="LEG",]
stopifnot(nrow(leg)==1L,abs(leg$Gap_pct-8.87993829)<0.001)
r <- read("04_downside_risk_summaries.csv")
stopifnot(all(r$N_common[r$Crop=="wheat"]==5L))
model <- read("12_primary_model_point_estimates.csv")
stopifnot(nrow(model)==9L)
message("Figure-data shape and selected numerical spot checks passed")
