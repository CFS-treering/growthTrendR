\name{cran-comments}
\alias{cran-comments}
\title{Comments to CRAN maintainers}
\description{
  This file summarizes changes and clarifications made to address CRAN
  comments for the \pkg{growthTrendR} package.
}

\details{

\strong{Package Title and Description}

- Updated \code{Title} to be under 65 characters:
  \emph{Toolkit for Data Processing, Quality, and Statistical Models}.
- All regular words (e.g., detrending, modeling) are unquoted.
  WinBuilder may give notes for words it does not recognize, but these are not
  package names, software names, or APIs.
- Updated \code{Description} to include references in CRAN-compliant format:
  - Wood (2017, ISBN:978-1498728331) for generalized additive models.
  - 'mgcv' package (<https://CRAN.R-project.org/package=mgcv>) for model implementation.
  - Girardin et al. (2021, <https://doi.org/10.1139/er-2020-0099>) for the CFS-TRenD repository.

\strong{Function documentation (.Rd files)}

- Added \code{\@return} sections to all exported functions where missing,
  describing the structure, class, and meaning of outputs.
- Specifically updated the following functions and their .Rd files:
  - \code{generate_report()} (\code{generate_report.Rd}): invisibly returns
    the file path of the generated report; the primary purpose is the side effect
    of creating the report.
  - \code{plot_freq()} (\code{plot_freq.Rd}): returns a list of \code{ggplot}
    objects corresponding to selected species.
  - \code{plot_qa()} (\code{plot_qa.Rd}): returns a list of four \code{ggplot}
    objects showing raw and treated series and their cross-correlations.
  - \code{plot_scale()} (\code{plot_scale.Rd}): returns a list of \code{ggplot}
    objects displaying scaled series and spatial contrasts.

\strong{Internal function \code{prepare_samples_clim()} used in examples}

- The function has been changed to exported to allow its use in examples.
- Corresponding examples have been added to its documentation.
- Unit tests have been created to verify correct behavior and output structure.

\strong{Example execution in documentation}

- All examples previously wrapped in dontrun have been unwrapped, except for CFS_mapping, which is now marked as donttest due to its execution time exceeding 5 seconds, following CRAN guidelines.

\strong{CRAN-compliant messages in source files}

- All instances of \code{print()} and \code{cat()} were replaced with
  \code{message()} or with \code{message()} combined with
  \code{utils::capture.output()} for tables.
- Affected R scripts include:
  \itemize{
    \item \code{R/quality_assessment.R}
    \item \code{R/utils-dev.R}
    \item \code{R/CFSformat.R}
    \item \code{R/misc.R}
  }

\strong{Preserving user settings in functions}

- Functions that temporarily modify graphics parameters (\code{par()}) now save
  the current settings first using \code{oldpar <- par(no.readonly = TRUE)}.
- An immediate \code{on.exit()} call ensures the original settings are restored
  when the function exits:
  \code{on.exit(par(oldpar), add = TRUE)}.
- Affected R script: \code{R/mapping.R}.
- Ensures user graphics settings are restored, following CRAN best practices.

}

