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



## ----data_loading-------------------------------------------------------------
# Example: load a dataset included in the package

# ring measurement
dt.samples <- fread(system.file("extdata", "dt.samples.csv", package = "growthTrendR"))

# formatting the users' data conformed to CFS-TRenD data structure
dt.samples_trt <- CFS_format(data = list(dt.samples, 39:68), usage = 1, out.csv = NULL)
class(dt.samples_trt)

# save it to extdata for further use
# saveRDS(dt.samples_trt, "inst/extdata/dt.samples_trt.rds")

## ----data_report, message = FALSE, warning = FALSE, results = 'hide'----------

outfile_data <- tempfile(fileext = ".html")
generate_report(robj = dt.samples_trt, qa.label_data = "demo-samples ", data_report.reports_sel = c(1,2,3,4), qa.min_nseries = 50, output_file = outfile_data)


## ----data_report_inclu, echo = FALSE, results = 'asis'------------------------
# cat('<div style="border: 5px solid #4682B4; border-radius: 8px; padding: 10px; 
#           margin: 15px 0; background-color: #f9f9f9; 
#           max-height: 400px; overflow-y: scroll;">')

htmltools::includeHTML(outfile_data)  # embed in vignette

