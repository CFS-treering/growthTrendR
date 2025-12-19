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



## ----workflow, echo = FALSE, results = 'asis', warning=FALSE, message=FALSE, out.width="80%", out.height="auto"----


# Path to the PNG inside the package

knitr::include_graphics("images/workflow_clean.png")
# Check existence

# Embed the image


