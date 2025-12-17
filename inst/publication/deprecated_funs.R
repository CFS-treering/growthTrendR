
#' tree growth model with gamm and model selection
#' @description
#' A wrapper function for the Generalized Additive Mixed Model (GAMM) to model annual tree growth, such as basal area increment.
#' It allows for modeling the log-transformed scale of tree growth.
#'
#' @param data data containing all necessary columns to run the model
#' @param resp_scale the scale of response variable. default is "log" for log-scale, otherwise modeling the response variable as it is
#' @param m.candidates the list of formulas.
#' @param out.csv directory of output csv files. default is NULL for no csv output.
#'
#'
#'
#' @import mgcv
#' @import nlme
#' @import stats
#' @import utils
#' @import data.table
#'
#'
#' @return list including model, fitting statistics, ptable, stable and prediction table
#' @details
#' The function implements tree identity as a random effect and incorporates a first-order autoregressive (AR1) correlation structure in the model.
#' To address the skewness of basal area increment (BAI), a Gamma family distribution with a log link function is employed.
#' When applied on the log scale of BAI, a normal distribution with an identity link function is employed,
#' assuming that the log-transformed BAI follows a normal distribution.
#'
#' If users specify multiple candidate models through the m.candiates argument, the function will fit each candidate model using the maximum likelihood (ML) method.
#' The Akaike Information Criterion (AIC) will then be compared to determine the best-fitting model. Once the optimal model is identified,
#' it will be refitted using the restricted maximum likelihood (REML) method.
#'
#' If users specify only 1 candidate model through the m.candiates argument, the model is fitted with "REML" method.

#' @export detrend_site
#'

detrend_site <- function(data, resp_scale = "log", m.candidates,  out.csv = NULL){
  if (length(m.candidates) == 0) stop("must assign the equation(s) to m.candidates")
  if (resp_scale == "log") {

    famil = gaussian("identity")
  }else {
    famil = Gamma("log")
  }
  # for comparing and selecting model on AIC
  if (length(m.candidates) > 1){
for (i in 1:length(m.candidates)){
  formul <- as.formula(m.candidates[i])
  if (resp_scale == "log") formul <- update(formul, log(.) ~ .)

  m.tmp <- gamm(formul,
                   random = list(uid_tree.fac=~1), correlation = corCAR1(value = 0.5),
                   method = "ML",family = famil, data = data)


  aic.tmp <- data.table(i = i, form = gsub("\\\\", "", paste(deparse(formul), collapse = " ")), aic =  AIC(m.tmp$lme), R2 = summary(m.tmp$gam)$r.sq, methd = "ML")


  if (i == 1) {
    aic.mn <- aic.tmp$aic
    m.sel <- m.tmp
    i.sel <- i
    aic.all <- aic.tmp
  } else{
    aic.all <- rbind(aic.all, aic.tmp)
    if (aic.mn > aic.tmp$aic){
      aic.mn <- aic.tmp$aic
      m.sel <- m.tmp
      i.sel <- i
    }
  }
# if (length(m.candidates) > 1) saveRDS(m.tmp, compress = TRUE, file =   paste0(clim.test,"/","m", i, ".",mtd,".rds"))
 rm(aic.tmp, m.tmp)
}
# aic.all$clim.test <- clim.test
aic.all[aic == min(aic), selected := "*"]
form.sel <- as.formula(aic.all[aic == min(aic)]$form)
# saveRDS(m.sel, compress = TRUE, file =   paste0(clim.test,"/", "sel.mod", " ",mtd,".rds"))
rm(m.sel)
# saveRDS(aic.all, compress = TRUE, file =   paste0(clim.test,"/","fitting ", mtd ,".rds"))
if (!is.null(out.csv)){

  if (!(dir.exists(out.csv))) dir.create(out.csv, recursive = TRUE)
write.csv(aic.all, file =  file.path(out.csv, paste0("fitting ML.csv")), row.names = FALSE, na = "")
  }


  }else{
  # for 1 equation only
  form.sel <- as.formula(m.candidates)
  if (resp_scale == "log") form.sel <- update(form.sel, log(.) ~ .)
}
  # fitting REML for prediction
  m.sel <- gamm(form.sel,
                random = list(uid_tree.fac=~1), correlation = corCAR1(value = 0.5),
                method = "REML",family = famil, data = data)


 pred.terms  <-as.data.frame( predict(m.sel$gam, type="terms",se.fit=TRUE))
  setnames(pred.terms, c("fit.s.ageC.", "se.fit.s.ageC."),c("fit.s.ageC", "se.fit.s.ageC"))
  fit.bai <- as.data.frame(predict(m.sel$gam, type = "response", se.fit = TRUE))
  names(fit.bai) <- c("fit.bai", "se.fit.bai")
  tmp.bai <- data.table(data, pred.terms, fit.bai)
  tmp.bai$res.bai <- residuals(m.sel$gam, type = "response")
  tmp.bai$res.bai_normalized <- residuals(m.sel$gam, type = "normalized")

  tmp.bai[, c("lci_bai", "uci_bai") := .(
    fit.bai - qnorm(0.975) * se.fit.bai,
    fit.bai + qnorm(0.975) * se.fit.bai)]

  # names(tmp.bai)
  # fit.s.age <- "fit.s.ageC..species.fac"



    # if (any(unique(str_detect(names(tmp.bai), fit.s.age)) == TRUE)){
  #     tmp.bai[, fit.s.ageC:=eval(parse(text= paste0(fit.s.age, "LARILAR"))) + eval(parse(text= paste0(fit.s.age, "PICEMAR")))]
  #     tmp.bai[, se.fit.s.ageC:=eval(parse(text= paste0("se.", fit.s.age, "LARILAR"))) + eval(parse(text= paste0("se.",fit.s.age, "PICEMAR")))]
  # # CI
      tmp.bai[, c("lci_ageC", "uci_ageC") := .(
        fit.s.ageC - qnorm(0.975) * se.fit.s.ageC,
        fit.s.ageC + qnorm(0.975) * se.fit.s.ageC)]
    # }
      if (resp_scale == "log") setnames(tmp.bai, c("fit.bai", "se.fit.bai", "res.bai", "lci_bai", "uci_bai"), c("fit.logbai", "se.fit.logbai", "res.logbai", "lci_logbai", "uci_logbai"))

  # for (iclim.test in c("VPDsummer","VPDprevsummer" )){
  # fit.s.clim <- paste0("fit.s.", iclim.test, "..species.fac")
  # if (any(unique(str_detect(names(tmp.bai), fit.s.clim)) == TRUE)){
  #   tmp.bai[, clim.fit:=eval(parse(text= paste0(fit.s.clim, "LARILAR"))) + eval(parse(text= paste0(fit.s.clim, "PICEMAR")))]
  #   tmp.bai[, se.clim.fit:=eval(parse(text= paste0("se.", fit.s.clim, "LARILAR"))) + eval(parse(text= paste0("se.", fit.s.clim, "PICEMAR")))]
  #   # CI
  #   tmp.bai[, c("lci_clim", "uci_clim") := .(
  #     clim.fit - qnorm(0.975) * se.clim.fit,
  #     clim.fit + qnorm(0.975) * se.clim.fit)]
  #
  # setnames(tmp.bai, c("clim.fit", "se.clim.fit", "lci_clim", "uci_clim"),
  #          c(paste0("fit.s.", iclim.test),paste0("se.fit.s.", iclim.test),paste0("lci_", iclim.test),paste0("uci_", iclim.test)))
  # }
  # }

  ptable <- data.table(Parameter = row.names(summary(m.sel$gam)$p.table), summary(m.sel$gam)$p.table )
  stable <- data.table(Parameter = row.names(summary(m.sel$gam)$s.table), summary(m.sel$gam)$s.table)
  aic.reml <- data.table(form = gsub("\\\\", "", paste(deparse(form.sel), collapse = " ")), aic =  AIC(m.sel$lme), R2 = summary(m.sel$gam)$r.sq, methd = "REML")

  if (!is.null(out.csv)){

    if (!(dir.exists(out.csv))) dir.create(out.csv, recursive = TRUE)
  write.csv(ptable, file =  file.path(out.csv, paste0("ptable ", "REML" ,".csv")), row.names = FALSE, na = "")
  write.csv(stable, file =  file.path(out.csv, paste0("stable ", "REML" ,".csv")), row.names = FALSE, na = "")
  write.csv(tmp.bai, file = file.path(out.csv, paste0("prediction ", "REML" ,".csv")), row.names = FALSE, na = "")
  write.csv(aic.reml, file = file.path(out.csv, paste0("fitting ", "REML" ,".csv")), row.names = FALSE, na = "")
  }


# if (make.plot) {
#   p.bai <- ggplot(tmp.bai, aes(x = year, y = fit.bai, color = species.fac)) + geom_line(size =1.2) +
#     geom_ribbon(aes(ymin=exp(lci_bai), ymax=exp(uci_bai), fill = species.fac), alpha=0.2) +
#     labs( y = paste0("bai"), x = "year") +
#     ggtitle( paste0("bai prediction") )
#
#   p.age <- ggplot(tmp.bai, aes(x = ageC, y = exp(fit.s.ageC), color = species.fac)) + geom_line(size =1.2) +
#     geom_ribbon(aes(ymin=exp(lci_ageC), ymax=exp(uci_ageC), fill = species.fac), alpha=0.2) +
#     labs( y = paste0("fit.s(", "ageC", ")"), x = "ageC") +
#     ggtitle( paste0("term response: ", "ageC") )
#
#
#
# }

# }

# return(i.sel)
return(list(model = m.sel, fitting = aic.reml, ptable = ptable, stable = stable, pred = tmp.bai))
}



# #' Generate table for radii report
# #'
# #' This function generates the data table for data reporting at project-species-site-radii level .
# #' @param tr_6i data table of meta data for a specific species
# #' @param tr_7i data table of tree ring data for a specific species
# #' @param rep_radii logic, for user to choose if the radii table is presented in the data report.
# #' @param ... additional arguments
# #' @param ... Additional arguments passed to \code{CFS_qa()}.
# #' @return a data table.
# #' @keywords internal
# #' @noRd
# #'
# table_spc_site_radii <- function(tr_6i, tr_7i, rep_radii = TRUE, ...){
#   args <- list(...)
#   qa.label_data <- args$qa.label_data
#   qa.min_nseries <- args$qa.min_nseries
# # if (rep_radii == TRUE & nrow(tr_6i) >= qa.min_nseries) {
#   # note thant nrow(tr_6i) checking in template_data_report.rmd 2025-08-16
#   if (rep_radii == TRUE) {
#   tr_7ispc<- copy(tr_7i)
#   # setnames(tr_7ispc, c("uid_radius", "year", "rw_mm"), c("SampleID", "Year" ,"RawRing"))
#   tr_7ispc[, c("SampleID", "Year" ,"RawRing"):= .(uid_radius, year, rw_mm)]
#   # check duplication
#  if (nrow( tr_7ispc[, .N, by = .(SampleID, Year)][N>1]) > 0) stop("SampleID-Year is not unique key, please check...")
#   # include RW_trt, this is the series that qa process works on,
#   # tr_7ispc[, RW_trt:= RawRing - shift(RawRing), by = SampleID]; label_trt <- "differentiated"
#   if (!("RW_trt" %in% names(tr_7ispc))) stop("please assign the column RW_trt for running the qa process...")
#   acf.trt <- tr_7ispc[!is.na(RW_trt), .( ar1_trt = round(acf(RW_trt, plot = FALSE)$acf[2],2) ), by = .(SampleID)]
#
#   # run quality assessment procedure
#   qa.trt.ccf <-CFS_qa(dt.input = tr_7ispc, qa.label_data = qa.label_data, qa.min_nseries = qa.min_nseries)
#   dt.radii.spc <- copy(qa.trt.ccf$dt.stats)
#
#   dt.radii.spc[, uid_radius:= as.integer(SampleID)]
#   dt.radii.spc <- merge(tr_6i[, c("uid_radius", "site_id", "radius_id")], dt.radii.spc, by = "uid_radius")
#
# }else{
#   # qa.trt.ccf <- list()
#   dt.radii.spc <- tr_7i[, .(N = .N, rw.mean = mean(rw_mm), rw.sd = sd(rw_mm), rw.min = min(rw_mm), rw.max = max(rw_mm), ymin = min(year), ymax = max(year), ar1_rw = acf(rw_mm, plot = FALSE)$acf[2] ), by = .(site_id, radius_id, uid_radius)]
#   # if more than 10 series, running the correlation with mean
#   if (nrow(tr_6i) > 10){
#     # mean.rw <- tr_7i[, mean_rw:= mean(rw_mm), by = uid_radius]
#     tmp <- merge(tr_7i, tr_7i[, .(mean_rw= mean(rw_mm)), by = year], by = "year")
#     chk <- tmp[, .(ncorr_mean_rw = .N, result = cor.test(mean_rw, rw_mm)), by = .(uid_radius)]
#
#     chk.corr <- chk[, .SD[4], by = uid_radius][, corr_mean_rw:= as.numeric(unlist(result))]
#
#     chk.pvalue <- chk[, .SD[3], by = uid_radius][, pcorr_mean_rw:= as.numeric(unlist(result))]
#
#     dt.radii.spc <- dt.radii.spc[chk.corr[, c("uid_radius", "corr_mean_rw", "ncorr_mean_rw")], on = "uid_radius"][chk.pvalue[, c("uid_radius", "pcorr_mean_rw")], on = "uid_radius"]
#   }
# }
#   return(dt.radii.spc)
# }





#' tree-ring data measurement assessment
#' @description
#' Assess tree-ring measurement accuracy using a treated series based on the differences between two consecutive tree-ring measurements.
#'
#' @param dt.input tree ring data with at least 3 columns (SampleID, Year, RawRing)
#' @param qa.label_data description of dt.input
#' @param qa.label_trt description of treated series
#' @param qa.batch_size  number of pairs to run in a batch, to avoid memory issues in processing large dataset
#' @param qa.max_lag maximum lag up to which the correlation should be calculated in CCF
#' @param qa.max_iter maximum number of iterations of step 2(see Details)
#' @param qa.min_nseries minimum number of series to run this function
#' @param qa.blcrit criteria for borderline

#'
#'
#'
#' @return A list of 5 elements:
#' 1) dt.ccf:	A data table containing the CCF results for all samples, including the quality assessment code (qa_code).
#' 2) dt.chron:	The final chronologies, including both raw chronologies, the mean of ring measurement of the series with pass,
#' and the treated  chronologies, calculated as the difference of two consecutive raw chronologies.
#' 3) dt.stats:	summary statistics of radii
#' 4) dt.plots, list of tables for generating plots
#' 5) qa.parms: parameters used


#'
#' @details Assess tree-ring measurement accuracy using a treated series based on the differences between two consecutive tree-ring measurements.
#' The algorithm consists of two main steps:
#'
#'  Step 1. Perform pairwise CCF on the treated series of all possible combinations of the samples.
#' raw ring chronologies is calculated as the average of the ring measurements of samples that achieve the maximum correlation at lag 0 with at least one other sample in the cross-correlation function (CCF).
#'  The difference series of the raw ring chronologies is used as the initial treated chronologies for next step.
#'
#'  Step 2. Perform CCF between the treated series of each sample and the initial treated chronologies.
#'  Samples that do not meet the criteria will be removed from the recalculation of the treated chronologies.
#'  This step is repeated until all remaining samples in the chronologies meet the criteria.
#'
#' This process results in five categories classification for all samples (pass, borderline, pm1, highpeak, fail)
#'
#'
#' To enhance efficiency and mitigate potential memory issues,
#' The function supports both parallel (multi-session) and sequential modes, and also offers a batch processing option for users.

# #' @export
#'
#' @keywords internal
#' @noRd
CFS_qa_old <- function(dt.input, qa.label_data = "", qa.label_trt = "",
                       qa.batch_size = 10000, qa.max_lag = 10, qa.max_iter = 100, qa.min_nseries = 100, qa.blcrit = 0.1){

  check_optional_deps()
  if (length(setdiff(c("species", "SampleID", "Year","RawRing", "RW_trt" ), names(dt.input))) > 0) stop("at least one of the mandatory columns (species, SampleID, Year, RawRing, RW_trt) doesn't exist, please check...")
  if ( qa.label_data == "") stop("please specify qa.label_data...")
  if ( "RW_trt" %in% names(dt.input) & qa.label_trt == "") stop("please specify qa.label_trt...")
  if (length(unique(dt.input$species)) > 1) stop("only 1 species is allowed in the dataset...")
  if (length(unique(dt.input$SampleID)) < qa.min_nseries) stop(paste0("please increase the sample size at minumum: ", qa.min_nseries))

  if (nrow(dt.input[, .N, by = .(SampleID, Year)][N > 1]) > 0) stop("SampleID-Year is not unique key, please check...")

  dt.rw_long <- dt.input[, c("SampleID", "Year","RawRing", "RW_trt")]
  setorder(dt.rw_long, SampleID,Year)
  # the series to be used for ccf should be in dt.input 2024-12-10
  # dt.rw_long[, RW_trt:= RawRing - shift(RawRing), by = SampleID]

  # SampleID.chr is the key sampleID, we use it in the functions to represent a sample.
  # starting with a character to ensure the validity as column name and value of a column

  dt.rw_long[, SampleID.chr:= paste0("d_", SampleID)]
  setorder(dt.rw_long, SampleID.chr)
  sample.lst <- sort(unique(dt.rw_long$SampleID.chr))

  dt.rw_wide <- dcast(dt.rw_long[!is.na(RW_trt)], Year ~ SampleID.chr, value.var = "RawRing")
  dt.trt_wide <- dcast(dt.rw_long[!is.na(RW_trt)], Year ~ SampleID.chr, value.var = "RW_trt")

  setcolorder(dt.trt_wide, c("Year", sample.lst))



  dt.trt_wide.o <- copy(dt.trt_wide); dt.trt_wide$Year <- NULL;
  # dt.trt_wide: treated series in wide format for pair-wise ccf
  # dt.rw_long: ring width series in long format for calculating mean of chronologies

  # step 1: pair-wise ccf to find all the samples which can find at least 1 sample to reach max_ccf @ lag0


  # Generate all pairs of columns
  col_pairs <- utils::combn(names(dt.trt_wide), 2, simplify = FALSE)

  # Detect available cores for parallel processing
  available_cores <- parallel::detectCores(logical = FALSE) - 1  # Adjusted cores based on system

  # Decide if parallel processing is supported
  if (available_cores > 1) {
    future::plan(future::multisession, workers = available_cores)
  } else {
    future::plan(future::sequential)
  }

  # if batch size is not specified, run as 1 batch
  if (is.null(qa.batch_size)) {
    pair_batches <- list(col_pairs)
  }else{
    pair_batches <- split(col_pairs, ceiling(seq_along(col_pairs) / qa.batch_size))
  }

  # Run processing on batches with or without parallel
  cat("Progress pair-wise ccf...\n")
  dt.ccf.pairs <- furrr::future_map(pair_batches, process_batch, dt.wide = dt.trt_wide, qa.max_lag = 10, .progress = TRUE) %>% rbindlist()


  cat("\n")

  result_dt.sel <- dt.ccf.pairs[max_lag == 0 & !is.na(max_ccf)]

  ts.sel <- data.table(ts.sel = unlist(c(result_dt.sel$ts1, result_dt.sel$ts2)))
  ts.sel <- ts.sel[, .N, by = .(ts.sel)]
  id.candi <- unique(ts.sel$ts.sel)

  # the result of step 1 is id.candi, it serves as the initial sample list of chronologies for step 2
  # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #


  # step 2: find the sample list of chronologies satisfying the condition that max_ccf @lag0 for each sample in this list with the chronologies

  # algorithm on pass

  s2.end <- FALSE;  i.iter <- 1;
  while(!s2.end & i.iter <= qa.max_iter){

    # mean of chronologies
    # dt.s2.avg <- dt.rw_long[SampleID.chr %in%  id.candi][, .(.N, mean.rw = mean(RawRing)), by = .(Year)]
    # setorder(dt.s2.avg, Year)
    # dt.s2.avg [, mean.rw.dif:= mean.rw - shift(mean.rw)]

    # mean of treated chronology as chronologies 2024-12-10
    dt.s2.avg <- dt.rw_long[SampleID.chr %in%  id.candi][, .(.N, mean.rw = mean(RawRing), mean.rw_trt = mean(RW_trt, na.rm = TRUE)), by = .(Year)]

    dt.s2.wide <- merge(dt.s2.avg[, c("Year", "mean.rw_trt")], dt.trt_wide.o, by = "Year")
    # ccf of all samples with the chronologies
    dt.s2.ccf <-rbindlist(lapply(3:ncol(dt.s2.wide), ccf_avg, data = dt.s2.wide, blcrit = qa.blcrit, lag.max = qa.max_lag, qa_code = "Fail"))
    # valid samples for chronologies
    id.pass <- unique(dt.s2.ccf[qa_code == "pass"]$SampleID.chr)
    # id.dif <- setdiff(union(id.candi, id.pass) ,intersect (id.candi, id.pass))
    s2.end <- length(setdiff(union(id.candi, id.pass) ,intersect (id.candi, id.pass))) == 0
    print(paste0(i.iter, " N.pass: ", length(id.candi)))
    # selection list stops changing
    if ( s2.end != TRUE) {
      id.candi <- id.pass
      i.iter <- i.iter + 1

    }
  }
  setorder(dt.s2.ccf, SampleID.chr)
  dt.s2.avg[, c("success", "iteration") := .(s2.end, i.iter)]

  # for generating plots
  # pre series data for both rw and treated in wide format
  dt.raw.series <- merge(dt.s2.avg[, c("Year", "mean.rw")], dt.rw_wide, by = "Year")
  dt.trt.series <- merge(dt.s2.avg[, c("Year", "mean.rw_trt")], dt.trt_wide.o, by = "Year")

  setcolorder(dt.raw.series, c("Year", "mean.rw", sample.lst))
  setcolorder(dt.trt.series, c("Year", "mean.rw_trt", sample.lst))


  # pre data for bar plotting on ccf with chronologies for raw and treated series

  dt.trt.ccf <- copy(dt.s2.ccf)
  dt.ccf.idlabel <- dt.trt.ccf[ccf.ord==1, c("SampleID.chr", "qa_code", "lag")]
  dt.ccf.idlabel[, id.label:= paste0(str_sub(SampleID.chr, 3, -1),"$", qa_code,"$", lag)]
  dt.trt.ccf<- merge(dt.trt.ccf, dt.ccf.idlabel[, c("SampleID.chr", "id.label")],by = "SampleID.chr")
  dt.trt.ccf <- dcast(dt.trt.ccf, lag ~ id.label, value.var = "acf.trt")
  names(dt.trt.ccf)
  idlabel.lst <- sort(unique(dt.ccf.idlabel$id.label))
  # test if in the same order as others
  idlabel.lst2 <- str_split_fixed(idlabel.lst, "\\$",3)[,1]
  if (!all.equal(idlabel.lst2,str_split_fixed(sample.lst, "\\_",2)[,2]  )) print("check the order of id.label in dt.ccf.idlabel")
  setcolorder(dt.trt.ccf, c("lag", idlabel.lst))

  # input data structure for ccf_avg Year, mean.value, sampleIDs...

  dt.raw.ccf <- rbindlist(lapply(3:ncol(dt.raw.series), ccf_avg, data = dt.raw.series))
  dt.raw.ccf <- dcast(dt.raw.ccf, lag ~ SampleID.chr, value.var = "acf.trt")
  names(dt.raw.ccf)
  setcolorder(dt.raw.ccf, c("lag", sample.lst))




  # for statistics per radius
  dt.s2.ccf[, SampleID := str_split_fixed(SampleID.chr, "\\_",2)[,2] ]


  # reports on radii
  dt.radii <- dt.rw_long[, .(N = .N, rw.mean = mean(RawRing), rw.sd = sd(RawRing), rw.min = min(RawRing), rw.max = max(RawRing), ymin = min(Year), ymax = max(Year), ar1_rw = acf(RawRing, plot = FALSE)$acf[2] ), by = .(SampleID.chr)]

  acf.trt <- dt.rw_long[!is.na(RW_trt), .( ar1_trt = round(acf(RW_trt, plot = FALSE)$acf[2],2) ), by = .(SampleID.chr)]

  # correlations
  dt_wide.rw <- merge(dt.s2.avg[, c("Year", "mean.rw")], dt.rw_wide, by = "Year")
  dt_wide.trt <- merge(dt.s2.avg[, c("Year", "mean.rw_trt")], dt.trt_wide.o, by = "Year")

  dt.cor <- merge(cor_avg(dt_wide.rw), cor_avg(dt_wide.trt), by = "SampleID.chr")

  stats_radii <- merge(dt.radii, acf.trt, by = "SampleID.chr")


  stats_radii <- merge(stats_radii, dt.cor, by = "SampleID.chr")


  stats_radii <- merge(dt.s2.ccf[lag == 0, c("SampleID", "SampleID.chr", "qa_code")], stats_radii, by = "SampleID.chr")


  dt.s2.ccf[, SampleID.chr := NULL]
  setcolorder(dt.s2.ccf, "SampleID")
  stats_radii[, SampleID.chr := NULL]

  # Reset to sequential
  future::plan(future::sequential)
  dt.s2.ccf <- data.table(species = unique(dt.input$species), dt.s2.ccf)
  dt.s2.avg <- data.table(species = unique(dt.input$species), dt.s2.avg)
  stats_radii <- data.table(species = unique(dt.input$species), stats_radii)
  qa_code <- data.frame(
    qa_code = c("pass", "borderline", "pm1", "highpeak", "fail"),
    Description   = c("The maximum correlation occurs at lag 0",
                      "The correlation at lag 0 ranks as the second highest, and its difference from the maximum remains within a predefined threshold, categorizing as a quasi-pass",
                      "The maximum correlation occurs at lag 1 or -1, suggesting slight misalignment.",
                      "The maximum correlation occurs at a non-zero lag and is more than twice the second-highest value, potentially signaling an issue",
                      "All other measurements that do not fit into the aforementioned categories fall under this classification.")
  )
  result <- list(dt.ccf = dt.s2.ccf, dt.chron = dt.s2.avg, dt.stats = stats_radii,
                 dt.plots = list(dt.trt.series = dt.trt.series, dt.raw.series = dt.raw.series, dt.trt.ccf = dt.trt.ccf, dt.raw.ccf =dt.raw.ccf),
                 qa.parms = list(qa.label_data = qa.label_data, qa.label_trt = qa.label_trt, qa.batch_size = qa.batch_size, qa.max_lag = qa.max_lag,
                                 qa.max_iter = qa.max_iter, qa.min_nseries = qa.min_nseries, qa.blcrit = qa.blcrit, qa.code_desc = qa_code))
  class(result) <- "cfs_qa"
  return(result)
  # the result of step 2 is dt.s2.ccf, samples with qa_code = "Pass" to form the chronologies

}


#' Generate mapped spatial-temporal plots
#'
#' This function interpolates and visualizes tree ring data (e.g., BAI growth change) across space and time
#' using inverse distance weighting (IDW). Outputs can include maps in PNG and GIF formats, interpolated raster
#' files in GeoTIFF format, and gridded values in CSV format.
#'
#' @param data A data.table object of class \code{cfs_mapping}, typically generated by \code{\link{read_series}}.
#' @param year.span A numeric vector of length two specifying the start and end years for analysis. Default is \code{c(1801, 2017)}.
#' @param extent.lim Optional. A numeric vector specifying map extent in the form \code{c(xmin, xmax, ymin, ymax)}. If \code{NULL}, extent is computed from the data.
#' @param grid.step A numeric value specifying the resolution (in degrees) of the interpolation grid. Default is \code{0.1}.
#' @param animation_fps Frames per second for the animated GIF output. Default is \code{1.0}.

#' @param by.spc Logical. If \code{TRUE}, runs separately for each species. Default is \code{FALSE}.
#' @param parms.out Character vector specifying the types of outputs to generate. Valid values are:
#'   \code{"csv"}, \code{"tif"}, \code{"png"}, and \code{"gif"}.
#' @param dir.out A character string specifying the directory where outputs will be saved. Required if \code{parms.out} is not empty.
#' @param data.crs Character string specifying the CRS (coordinate reference system) for spatial outputs. Default is WGS84.
#' @param ... Additional parameters passed to the function.
#'
#' @return This function does not return a value but saves outputs to disk depending on the options selected in \code{parms.out}.
#'   Output includes:
#'   \itemize{
#'     \item PNG maps of interpolated tree growth
#'     \item Animated GIFs showing change over time
#'     \item Raster GeoTIFFs of interpolated values
#'     \item CSV files of gridded predicted values
#'   }
#'
#' @details
#' The function performs spatial interpolation using IDW and masks outputs to Canada and boreal zones using shapefiles.
#' The elevation layer and mask files are internal to the \code{growthTrendR} package and used for contextual mapping.
#'
#' PNG and GIF outputs are generated using \code{magick}, raster outputs via \code{terra}, and interpolation via \code{gstat}.
#'
#' @section File Outputs:
#' When parms.out is specified, the following directory structure is created:
#' \itemize{
#'   \item \strong{png/}: Individual PNG files for each species and year
#'   \item \strong{gif/}: Animated GIF files for each species
#'   \item \strong{csv/}: CSV files with annual data for each species
#'   \item \strong{tif/}: GeoTIFF files for each species and year
#' }
#'
#' @seealso
#' \code{\link{CFS_mapping}} for preparing input data
#'
#' @examples
#' \dontrun{
#' # Basic usage - return animated GIFs
#' gif_results <- plot_mapping(cfs_data)
#'
#' # Generate all output formats
#' plot_mapping(cfs_data,
#'               animation_fps = 2.0,
#'               parms.out = c("png", "gif", "csv"),
#'               dir.out = "output_folder")
#'
#'
#' # Custom animation speed and plot parameters
#' plot_mapping(cfs_data,
#'               animation_fps = 0.5,
#'               parms.out = "gif",
#'               dir.out = "slow_gifs")
#' }
#'

#' @export


plot_mapping_old <- function(data, year.span = c(1801,2017),
                             # cols.meta = c("uid_tree", "uid_site", "longitude", "latitude", "species", "year"),
                             extent.lim = NULL, grid.step = 0.1 , animation_fps = 1.0,
                             # png.text = list(text1= "BAI Annual Growth Change -", text_bott = "Created by Martin P. Girardin, Canadian Forest Service (2025)" , text_side = "Growth Change (%)"),
                             by.spc = FALSE, parms.out = c("", "csv","tif","png", "gif"), dir.out= NULL,
                             data.crs = "+proj=longlat +datum=WGS84 +no_defs", ...
){

  check_optional_deps()

  if (!inherits(data, "cfs_mapping")) stop("please check the input of data, make sure it's the result of read_series() function")
  if (length(setdiff(parms.out, c("csv","tif","png", "gif"))) > 0) stop('please check parms.out, only "csv","tif","png", "gif" are supported' )
  if (is.null(dir.out) & !is.null(parms.out)) stop("please specify the output directory (dir.out) ")
  if (by.spc == TRUE) dir.out <- file.path(dir.out, "by_spc")
  if (!dir.exists(dir.out)) dir.create(dir.out, recursive = TRUE)
  data <- data$dt.w

  if (by.spc == TRUE)  data[, species.inuse:= species] else
    data$species.inuse <- "all.spp"
  year_cols <- grep("^\\d+$", str_sub(names(data), 2, -1), value = TRUE)

  # Get non-year columns
  non_year_cols <- setdiff(names(data), paste0("X", year_cols))
  val.years <- intersect(year_cols, as.character(year.span[1]:year.span[2]))
  data.yr <- data[ ,c(non_year_cols, paste0("X", val.years)), with = FALSE]


  if (is.null(extent.lim)){
    # Create a grid covering Canada
    canada_ext <- extent(min(data$longitude), max(data$longitude),
                         min(data$latitude), max(data$latitude))}else
                           canada_ext <- extent(extent.lim)


  grid <- expand.grid(
    longitude = seq(from = canada_ext@xmin, to = canada_ext@xmax, by = grid.step),
    latitude = seq(from = canada_ext@ymin, to = canada_ext@ymax, by = grid.step)
  )

  grid$grid_cell <- 1:nrow(grid)
  sp::coordinates(grid) <- ~longitude+latitude
  sp::gridded(grid) <- TRUE
  final_df <- as.data.frame(grid)

  ncells_total <- nrow(final_df)

  elevation_file <- system.file("extdata/Mapping/Canada120s.tif", package = "growthTrendR")

  # Read elevation and shapefiles
  elevation <- rast(elevation_file)


  # using the figshare that martin provided to avoid big extdata
  # the reference is here https://doi.org/10.6084/m9.figshare.30005473.v1 (Martin 2025-10-01) then redirect to
  # https://figshare.com/articles/dataset/Spatially_detailed_tree-ring_analysis_throughout_Canada/30005473
  # the true url is by right clicking download button to copy link address for the file
  # url <- "https://figshare.com/ndownloader/files/57487996"


  boreal_mask <- get_boreal_mask("https://figshare.com/ndownloader/files/57487996")
  canada_mask <- sf::st_read(system.file("extdata/Mapping/province.shp", package = "growthTrendR"))

  spc.lst <- sort(unique(data$species.inuse))
  # run over all species
  results_list.spc <- lapply(spc.lst, function(spc.i){
    dt.spc.i <- data.yr[species.inuse == spc.i]
    # remove the columns with NA only
    setDF(dt.spc.i)
    dt.spc.i <- dt.spc.i[, colSums(!is.na(dt.spc.i)) > 0]
    setDT(dt.spc.i)
    # run over year
    lst <- list()
    cols <-setdiff(names(dt.spc.i), non_year_cols)
    names(cols) <- str_sub(cols, 2)

    results_list.ispc <- lapply(cols, function(col) {
      yr <- as.numeric(str_sub(col,2))

      current_data <- dt.spc.i[!is.na(dt.spc.i[[col]]), c("longitude", "latitude", col), with = FALSE]

      # Ensure spatial coordinates
      sp::coordinates(current_data) <- ~longitude + latitude
      # sp::coordinates(grid) <- ~longitude + latitude

      idw_model <- gstat::idw(as.formula(paste(col, "~1")),
                              locations = current_data,
                              newdata = grid,
                              idp = 2,
                              nmax = 100)

      df.idw <- cbind( year = yr, as.data.table(idw_model)[, -"var1.var"])
      r <- raster(idw_model)

      # # mask
      # r_crop <- crop(r, canada_mask_proj)
      # r_mask <- mask(r_crop, canada_mask_proj)
      if (any(str_detect(parms.out, "tif")) == TRUE) {

        dir.out1<- file.path(dir.out, "tif")
        if (!dir.exists(dir.out1)) dir.create(dir.out1)
        writeRaster(r,
                    filename = file.path(dir.out1,
                                         paste0(spc.i, " ","tree_rings_", yr, ".tif")),
                    format = "GTiff",
                    overwrite = TRUE)
      }
      crs(r) <- data.crs
      # boreal_mask_proj <- sf::st_transform(boreal_mask, data.crs)
      # canada_mask_proj <- sf::st_transform(canada_mask, data.crs)
      # r_boreal <- mask(rast(r), vect(boreal_mask_proj))
      # r_masked <- mask(r_boreal, vect(canada_mask_proj))  # Second masking step with Canada

      # --- Mask rasters using terra ---
      # Convert RasterLayer to SpatRaster if needed
      if (inherits(r, "RasterLayer")) r <- terra::rast(r)

      # Transform boreal and Canada masks to the target CRS and convert to SpatVector
      boreal_mask_proj <- terra::vect(sf::st_transform(boreal_mask, data.crs))
      canada_mask_proj <- terra::vect(sf::st_transform(canada_mask, data.crs))

      # Apply masks sequentially
      r_boreal <- terra::mask(r, boreal_mask_proj)
      r_masked <- terra::mask(r_boreal, canada_mask_proj)

      # # Now r_masked is a SpatRaster compatible with terra::values()
      # values_vector <- terra::values(r_masked)


      # Prepare output path if PNG is requested
      output_png <- NULL
      if (any(str_detect(parms.out, "png"))) {
        dir.png <- file.path(dir.out, "png")
        dir.create(dir.png, showWarnings = FALSE, recursive = TRUE)
        output_png <- file.path(dir.png, paste0(spc.i, " tree_rings_", yr, ".png"))
      }

      # Call the function once, with or without output path
      p <- plot_tree_ring_map(r_masked, elevation, yr, out.png = output_png, ...)

      # Save output if GIF is requested
      if (any(str_detect(parms.out, "gif"))) {
        lst$png <- p
      }



      # Return result with year attached in csv per species with all year
      if (any(str_detect(parms.out, "csv"))) {
        cat("Adding data for year:", col, "\n")
        values_vector <- terra::values(r_masked)
        ncells_total <- ncell(r_masked)
        coords <- xyFromCell(r_masked, 1:ncells_total)
        final_df <- data.frame(grid_cell = 1:ncells_total, X = coords[, "x"], Y = coords[, "y"])
        if (length(values_vector) == ncells_total) {
          final_df.ispc <- data.frame(final_df, year = yr, var1.pred = values_vector)
        } else {
          temp_vec <- rep(NA, ncells_total)
          non_na_indices <- which(!is.na(values_vector))
          if (length(non_na_indices) > 0) {
            valid_indices <- non_na_indices[non_na_indices <= ncells_total]
            temp_vec[valid_indices] <- values_vector[valid_indices]
          }
          final_df.ispc <- data.frame(final_df, year = yr, var1.pred = temp_vec)
        }
        final_df.clean <- final_df.ispc[!is.na(final_df.ispc[,"var1.pred"]),]


        lst$df.idw <- final_df.clean

      }

      if (length(lst) > 0) return(lst)
    }
    ) # year loop end



    if (any(str_detect(parms.out, "csv")) == TRUE) {
      # Return result with year attached in csv per species with all year
      df.idw.ispc <- rbindlist(lapply(results_list.ispc, `[[`, "df.idw"), fill = TRUE)
      df.test <- rbindlist(lapply(results_list.ispc, `[[`, "df.idw.comp"), fill = TRUE)

      df.idw.ispc.all <- rbindlist(lapply(results_list.ispc, `[[`, "df.idw.all"), fill = TRUE)
      # Optional: reshape to wide format
      df.idw.ispc.wide <- dcast(df.idw.ispc, grid_cell ~ year, value.var = "var1.pred")
      df.idw.ispc.wide<- merge(final_df, df.idw.ispc.wide, by = "grid_cell", all.y = TRUE)
      dir.csv<- file.path(dir.out, "csv")
      if (!dir.exists(dir.csv)) dir.create(dir.csv)
      utils::write.csv(df.idw.ispc.wide, file = file.path(dir.csv, paste0(spc.i, " ", " by degree-", grid.step," y",  year.span[1], "-", year.span[2], ".csv")), row.names = FALSE, na = "")
      # utils::write.csv(df.idw.ispc.all, file = file.path(dir.csv, paste0(spc.i, " ", " by degree-", grid.step," y",  "idw.model", ".csv")), row.names = FALSE, na = "")

    }
    years <- names(results_list.ispc)  # or however you store years
    img_list <- unlist(lapply(results_list.ispc[order(as.numeric(years))], `[[`, "png"))

    gif_imgs <- magick::image_read(img_list)
    gif_animated <- magick::image_animate(gif_imgs, fps = animation_fps)

    if (any(str_detect(parms.out, "gif")) == TRUE) {

      dir.out1<- file.path(dir.out, "gif")
      if (!dir.exists(dir.out1)) dir.create(dir.out1)
      # img_list <- unlist(lapply(results_list.ispc, `[[`, "png"))
      # img_list <- sort(img_list)


      gif_animated %>%
        magick::image_write(file.path(dir.out1, paste0(spc.i, " ","tree_rings_animation y", year.span[1],"-", year.span[2], ".gif")))
    }
    return(gif_animated)

  })# species lapply end
  names(results_list.spc) <- spc.lst
  return(results_list.spc)
} # function end



#' Plot Tree-Ring Growth Map with Elevation
#'
#' Creates a static PNG map overlaying tree-ring growth changes on top of a grayscale elevation raster,
#' including a customized color legend and annotation.
#'
#' @param tree_ring_raster A raster object containing the tree-ring growth change values (e.g., BAI ratio).
#' @param elevation_raster A raster object representing elevation data for background shading.
#' @param year Integer. The year to display in the map title.
#' @param out.png Optional. File path for saving the PNG map. If `NULL`, a temporary file is created and returned.
#' @param png.text A named list of character strings used to annotate PNG maps. Elements include:
#'   \itemize{
#'     \item \code{text_top}: Main top title (e.g., "BAI Annual Growth Change -")
#'     \item \code{text_bott}: Bottom caption (e.g., author or data source)
#'     \item \code{text_side}: Legend or axis label (e.g., "Growth Change (%)")
#'   }
#'
#' @return The file path of the created PNG image if `out.png` is not provided.
#' @keywords internal
#' @noRd

plot_tree_ring_map <- function(tree_ring_raster, elevation_raster, year, out.png = NULL,
                               png.text = NULL, tick_positions = 0:10, tick_labels ) {

  # png.text <- utils::modifyList(list(text_top= "BAI Annual Growth Change -", text_bott = "Created by Martin P. Girardin, Canadian Forest Service (2025)" , text_side = "Growth Change (%)"), png.text)
  default_text <- list(
    text_top = "BAI Annual Growth Change -",
    text_bott = "Created by Martin P. Girardin, Canadian Forest Service (2025)",
    text_side = "Growth Change (%)"
  )

  # Only merge if png.text is provided and is a list
  if (is.list(png.text))  png.text <- utils::modifyList(default_text, png.text)
  if (is.null(png.text))  png.text <- default_text
  if (!is.list(png.text) & !is.null(png.text))
    stop("png.text must be a list with elements: text_top, text_bott, text_side")

  # Determine output file
  if (!is.null(out.png)) {
    output_file <- out.png
  } else {
    output_file <- tempfile(fileext = ".png")
  }

  png(output_file, width = 2500, height = 1300, res = 300)
  par(mfrow = c(1, 1))
  layout(matrix(c(1, 2), ncol = 2), widths = c(4, 1))
  par(mar = c(4, 4, 1, 3))

  # if (!identical(terra::crs(elevation_raster), terra::crs(tree_ring_raster))) {
  #   elevation_raster <- terra::project(elevation_raster, tree_ring_raster)
  # }
  # elevation_raster <- raster::crop(elevation_raster, tree_ring_raster)
  # elevation_raster[elevation_raster < 1] <- NA
  # Ensure both are SpatRaster
  elevation_raster   <- terra::rast(elevation_raster)
  # tree_ring_raster   <- terra::rast(tree_ring_raster)

  # Reproject elevation if CRS differs
  if (!terra::same.crs(elevation_raster, tree_ring_raster)) {
    elevation_raster <- terra::project(elevation_raster, tree_ring_raster)
  }

  # Crop elevation to tree-ring extent
  elevation_raster <- terra::crop(elevation_raster, tree_ring_raster)

  # Mask invalid elevation values
  # elevation_raster[elevation_raster < 1] <- NA


  #
  #   terra::plot(elevation_raster, col = gray.colors(100, start = 0.8, end = 0.3), legend = FALSE,
  #        main = paste(png.text[["text_top"]], year), axes = TRUE)
  # Only plot elevation if it has values
  if (!is.null(elevation_raster)  && terra::hasValues(elevation_raster)) {
    elevation_raster <- terra::ifel(elevation_raster < 1, NA, elevation_raster)
    terra::plot(elevation_raster,
                col = gray.colors(100, start = 0.8, end = 0.3),
                legend = FALSE,
                main = paste(png.text[["text_top"]], year),
                axes = TRUE)
    add_overlay <- TRUE
  } else {
    # Create empty plot as base
    plot.new()
    plot.window(xlim = c(0, 1), ylim = c(0, 1))
    title(main = paste(png.text[["text_top"]], year))
    add_overlay <- TRUE
  }


  r_clipped <- raster::clamp(tree_ring_raster, lower = 0.5, upper = 1.5, values = TRUE)
  terra::plot(r_clipped,
              col = colorRampPalette(c("red", "yellow", "green", "blue"))(100),
              alpha = 0.6,
              add = TRUE,
              legend = FALSE)
  grid()
  mtext(png.text[["text_bott"]],
        side = 1, line = 3, adj = 0, cex = 0.8)

  par(mar = c(5, 0, 3, 3))
  raster::plot(c(0, 1), c(0.3, 1.7), type = "n", axes = FALSE, xlab = "", ylab = "")

  cols <- colorRampPalette(c("darkred", "red", "yellow", "green", "blue"))(100)
  breaks <- seq(0, 10, length.out = length(cols) + 1)

  for(i in 1:length(cols)) {
    rect(0.1, breaks[i], 0.4, breaks[i + 1], col = cols[i], border = NA)
  }

  # tick_positions <- seq(0.5, 1.5, by = 0.1)
  # tick_labels <- paste0((tick_positions - 1) * 100)
  axis(4, at = tick_positions, labels = tick_labels, las = 1, cex.axis = 0.8, pos = 0.5)
  mtext(png.text[["text_side"]], side = 4, line = 1.5, cex = 1.2)

  dev.off()

  # Return the output file path
  if (is.null(out.png)) return(output_file)
}

