## ----setup_general, include = FALSE-------------------------------------------
knitr::opts_chunk$set(
   echo = TRUE,
  warning = FALSE,
  message = FALSE,
  collapse = TRUE,
  comment = "#>"
)

## ----setup_lib, include=FALSE-------------------------------------------------
if (!requireNamespace("htmltools", quietly = TRUE)) {
  stop("Package 'htmltools' is required to run this function. Please install it.")
}
library(growthTrendR)
library(data.table)
library(ggplot2)
library(sf)



## ----mod_data,  message = FALSE, warning = FALSE, results = 'hide'------------


# loading processed ring measurement
dt.samples_trt <- readRDS(system.file("extdata", "dt.samples_trt.rds", package = "growthTrendR"))

# climate
dt.clim <- fread(system.file("extdata", "dt.clim.csv", package = "growthTrendR"))

# merge data
dt.samples_clim <- merge(dt.samples_trt$tr_all_wide[, c("uid_site", "site_id","latitude", "longitude",  "species", "uid_tree", "uid_radius")], dt.samples_trt$tr_all_long$tr_7_ring_widths, by = "uid_radius")

# # Calculate BAI
dt.samples_clim <- calc_bai(dt.samples_clim, "uid_radius", "rw_mm")

dt.samples_clim <- merge(dt.samples_clim, dt.clim, by = c("site_id", "year"))

## ----mod_fitting,  message = FALSE, warning = FALSE, results = 'hide'---------
setorder(dt.samples_clim, uid_tree, year)

# remove the first 10 years, at least remove ageC == 1, for the log-scale model functions
m.sp <- gamm_spatial(data = dt.samples_clim[ageC > 1], resp_scale = "resp_log",
                     m.candidates =c( "bai_cm2 ~ log(ba_cm2_t_1) + s(ageC) + s(FFD)",
                                      "bai_cm2 ~ log(ba_cm2_t_1) + s(ageC) + FFD")
)

## ----mod_report_demo, eval=FALSE, message = FALSE, warning = FALSE, results = 'hide'----
# 
# generate_report(robj = m.sp)

## ----mod_report, include=FALSE, message = FALSE, warning = FALSE, results = 'hide'----

outfile_mod <- tempfile(fileext = ".html")
generate_report(robj = m.sp, output_file = outfile_mod)

## ----mod_inclu, echo = FALSE, results = 'asis'--------------------------------

cat('
<div style="margin: 20px 0; padding: 10px; 
            box-shadow: 0 2px 8px rgba(0,0,0,0.1); 
            border-radius: 8px;">
')

htmltools::includeHTML(outfile_mod)  # embed in vignette

