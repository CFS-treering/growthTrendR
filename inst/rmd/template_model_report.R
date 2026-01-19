params <-
list()

## ----setup_model, include=FALSE-----------------------------------------------
knitr::opts_chunk$set(echo = FALSE)
library(data.table)
library(mgcv)
library(ggplot2)
library(dplyr)
library(stringr)
library(growthTrendR)
library(ggeffects)
# library(cowplot)
library(patchwork)


## ----warning=FALSE, message=FALSE---------------------------------------------

robj <- params$robj

## -----------------------------------------------------------------------------
if (is.list(robj$model) && all(c("gam", "lme") %in% names(robj$model))) gam_model <- robj$model$gam else gam_model <- robj$model
summary(gam_model)

if (nrow(robj$stable) > 0){
term_important <- sterm_imp( gam_model)

method <- unique(term_important$method)
}else method <- ""


## ----echo=FALSE, results='asis', warning=FALSE, message=FALSE, fig.width=12, fig.height=8, fig.align='left'----
if (nrow(robj$stable) > 0){
# term_important <- sterm_imp( gam_model)
# 
# method <- unique(term_important$method)
tab <- term_important[, c(1,2)]
colnames(tab) <- c("Term", "Score (%)")
knitr::kable(tab, caption = "Relative importance of smooth terms",
      align = c("l", "r"))
}else{
  cat("\n**no smooth term in the model...**\n")
  cat("<br>")
}

## ----echo=FALSE, results='asis', warning=FALSE, message=FALSE, fig.width=12, fig.height=8, fig.align='left'----
if (nrow(robj$stable) > 0){
if ("bam" %in% class(robj$model)) {
  cat("\n**Predicting smooth terms on the response scale is not yet supported for bam models. The plots shown here are based on mgcv::plot.gam()**\n")
  cat("<br>")
  plot.gam(robj$model, pages = 1)
}else{
plot_list <- plot_resp(robj = robj)

# combined <- plot_grid(
#   plotlist = plot_list,
#   ncol = 2,
#   labels = LETTERS[1:length(plot_list)]
# )

combined <- wrap_plots(plot_list, ncol = 2) + 
  plot_annotation(tag_levels = 'a',
      tag_suffix = ")")

print(combined)
}
}else{
  cat("\n**no smooth term in the model...**\n")
  cat("<br>")
}

## ----echo=FALSE,   warning=FALSE, message=FALSE, fig.width=12, fig.height=8, fig.align='left'----

 # Check smoothness, residuals, K-index

 # 2x2 layout on one page
par(mfrow = c(2, 2))
  if (is.list(robj$model) && all(c("gam", "lme") %in% names(robj$model))) gam_model <- robj$model$gam else gam_model <- robj$model
gam.check(gam_model, rep = 100, verbose = TRUE)

# reset plotting layout to default
par(mfrow = c(1, 1))
 


