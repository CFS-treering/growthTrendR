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
if (!requireNamespace("magick", quietly = TRUE)) {
  stop("Package 'magick' is required to run this function. Please install it.")
}
library(growthTrendR)
library(data.table)
library(ggplot2)
library(sf)



## ----data_qa, echo= FALSE, eval = FALSE, eval=FALSE---------------------------
# 

## ----scale--------------------------------------------------------------------
# loading processed data
# otherwise need to run CFS_format() first as done in data report
dt.samples_trt <- readRDS(system.file("extdata", "dt.samples_trt.rds", package = "growthTrendR"))
# , message = FALSE, warning = FALSE, results = 'hide'
all.sites <- dt.samples_trt$tr_all_wide[,.N, by = c("species", "uid_site", "site_id")][, N:=NULL]
if (nrow(all.sites[, .N, by = .(species, site_id)][N>1]) > 0) stop("species-site_id is not unique...")
# e.g. taking the target sites
target_site <- all.sites[c(1,2), -"uid_site"]

ref.sites <- merge(dt.samples_trt$tr_all_wide[,c("species", "uid_site", "site_id", "latitude","longitude", "uid_radius")], dt.samples_trt$tr_all_long$tr_7_ring_widths, by = c("uid_radius"))

dt.scale <- CFS_scale( target_site = target_site, ref_sites = ref.sites, scale.label_data_ref = "CFS-TRenD V1.2-proj69", scale.max_dist_km = 200, scale.N_nbs = 2)

## ----scale_report_demo, eval= FALSE-------------------------------------------
# generate_report(robj = dt.scale, output_file = outfile_scale)

## ----scale_report, include = FALSE, message = FALSE, warning = FALSE, results = 'hide'----
outfile_scale <- tempfile(fileext = ".html")
generate_report(robj = dt.scale, output_file = outfile_scale)

## ----scale_inclu, echo = FALSE, results = 'asis'------------------------------
# cat('<div style="border: 5px solid #4682B4; border-radius: 8px; padding: 10px; 
#           margin: 15px 0; background-color: #f9f9f9;">')
# cat('<div style="border: 5px solid #4682B4; border-radius: 8px; padding: 10px; 
#           margin: 15px 0; background-color: #f9f9f9; 
#           max-height: 400px; overflow-y: scroll;">')


cat('
<div style="margin: 20px 0; padding: 10px; 
            box-shadow: 0 2px 8px rgba(0,0,0,0.1); 
            border-radius: 8px;">
')

htmltools::includeHTML(outfile_scale)  # embed in vignette
cat('</div>')

## ----qaa, message = FALSE, warning = FALSE, results = 'hide'------------------

# loading processed data
dt.samples_trt <- readRDS(system.file("extdata", "dt.samples_trt.rds", package = "growthTrendR"))


# data processing
dt.samples_long <- merge(dt.samples_trt$tr_all_wide[, c("uid_site", "site_id", "species", "uid_tree", "uid_sample", "sample_id", "radius_id", "uid_radius")],
                        dt.samples_trt$tr_all_long$tr_7_ring_widths, by = "uid_radius")

# rename to the reserved column name
setnames(dt.samples_long, c("sample_id", "year", "rw_mm"), c("SampleID", "Year" ,"RawRing"))

# assign treated series
# users can decide their own treated series
dt.samples_long[, RW_trt:= RawRing - shift(RawRing), by = SampleID]

# quality check on radius alignment based on the treated series
dt.qa <-CFS_qa(dt.input = dt.samples_long, qa.label_data = "CFS-TRenD V1.2-proj69", qa.label_trt = "difference", qa.min_nseries = 5)

## ----qa_report_demo, eval=FALSE, message = FALSE, warning = FALSE, results = 'hide'----
# 
# # e.g. series to check
# chk.lst <- dt.qa$dt.stats[1:2,]$SampleID
# 
# generate_report(robj = dt.qa,  qa.out_series = chk.lst)

## ----qa_report, include= FALSE, message = FALSE, warning = FALSE, results = 'hide'----
chk.lst <- dt.qa$dt.stats[1:2,]$SampleID 

outfile_qa <- tempfile(fileext = ".html")
generate_report(robj = dt.qa,  qa.out_series = chk.lst, output_file = outfile_qa)

## ----qa_inclu, echo = FALSE, results = 'asis'---------------------------------
# cat('<div style="border: 5px solid #4682B4; border-radius: 8px; padding: 10px; 
#           margin: 15px 0; background-color: #f9f9f9; 
#           max-height: 400px; overflow-y: scroll;">')

cat('
<div style="margin: 20px 0; padding: 10px; 
            box-shadow: 0 2px 8px rgba(0,0,0,0.1); 
            border-radius: 8px;">
')

htmltools::includeHTML(outfile_qa)  # embed in vignette
cat('</div>')

