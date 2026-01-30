params <-
list()

## ----setup_qa, include=FALSE--------------------------------------------------
# This chunk sets global options to hide all R code from appearing in the output

knitr::opts_chunk$set(echo = FALSE)

library(data.table)
library(growthTrendR)
library(gt)
library(stringr)
library(patchwork)
library(ggplot2)
# robj <- qa.trt.ccf; out.series.sel <- c("X011_104_003", "X016_104_003")
# robj <- params$robj; out.series.sel = params$parms.reports["qa.out_series"]
# print(params$parms["qa.out_series"])
# plt.lst <- plot_qa(params$robj, out.series = out.series.sel)
robj <- params$robj
qa_code <- robj$qa.parms$qa.code_desc


## ----results='asis', warning=FALSE, message=FALSE, fig.width=10, fig.height=5, fig.align='left'----
gt(qa_code) %>%
  tab_header(
    title = "Description of qa_code"
  )%>%
  tab_style(
    style = cell_text(weight = "bold"),
    locations = cells_body(columns = vars(qa_code))
  ) %>%
  opt_css(css = "td:nth-child(1) { padding-right: 20px; }")


## ----results='asis', warning=FALSE, message=FALSE, fig.width=10, fig.height=5, fig.align='left'----
plt.lst <- plot_qa(robj, qa.out_series = params$qa.out_series)



    # Loop through the plot lists and arrange them on each page
   for (i in 1:length(plt.lst$plot.trt.series)) {
  # plots <- list(plot.lst$plot.raw.series[[i]], plot.lst$plot.trt.series [[i]], plot.lst$plot.raw.ccf[[i]], plot.lst$plot.trt.ccf[[i]])
  # grid.arrange(grobs = plots, ncol = 2, nrow = 2)
  p.i<-  (plt.lst$plot.raw.series[[i]] | plt.lst$plot.trt.series [[i]]) / (plt.lst$plot.raw.ccf[[i]] | plt.lst$plot.trt.ccf[[i]]) +
    plot_annotation(
      # title = paste0("SampleID: ", samples.all[idx.lst[1]]) ,
      caption = paste0("Data source: ", robj$qa.parms$qa.label_data),
      tag_levels = 'a',
      tag_suffix = ")") &
    theme(
      plot.title = element_text(face = "bold"),
      plot.tag = element_text(face = "bold"),
      plot.caption = element_text(hjust = 0.2, face = "italic" )
      
      # plot.margin = margin(t = 10, r = 10, b = 30, unit = "pt") # Adjust plot margins
    )
  print(p.i)
}


