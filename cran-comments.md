## CRAN comments — growthTrendR

This file summarizes changes made in response to CRAN feedback 0.1.1.

### Package Title and Description
- Updated the package Title to be under 65 characters.
- Removed quotes around regular words in the Description.
- Added CRAN-compliant references with DOI/ISBN/CRAN links
  (Wood 2017; mgcv package; Girardin et al. 2021).

### Function documentation
- Added @return sections to all exported functions, describing
  output structure and classes.


### Functions used in examples
- Exported prepare_samples_clim() to allow use in examples.
- Added corresponding examples and basic unit tests.

### Examples wrapped in \dontrun{}
- Fast examples were unwrapped; longer-running examples were
  rewrapped in \donttest{} following CRAN guidelines.

### User-facing messages
- Replaced print() and cat() with message() throughout the package.

### Preserving user settings
- using an immediate on.exit() call before change the setting by par().
- Affected code is located in R/mapping.R.

### VECTOR_ELT fix (issue on macos only)
- Fixed a potential VECTOR_ELT issue in grouped data.table
  assignments by replacing length-dependent RHS expressions
  in := + by operations with a group-invariant approach.

