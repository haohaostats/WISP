# Implementation provenance

The core algorithms were adapted from the authors' research scripts `wisp_functions.R` and `ipcw_functions.R`. Research files remain unchanged outside this package.

Package changes remove script-level library attachment and code comments, define namespace imports, add explicit input checks and Rd help, and add returned estimated dropout probabilities, model fit and stacked covariance. IPCW log probabilities use negative sums of hazards, and their score derivatives use cumulative hazards directly. These are algebraic equivalents that avoid cancellation and zero-times-infinity in unused assessment periods; weights are not clipped.

`tests/research-reference.rds` contains only 600 simulated subjects (seed 42) and outputs computed using the original research scripts. It contains no participant data. `tests/regression.R` verifies estimates, HC0 SE, stacked SE and coefficient influence functions against these frozen outputs, alongside independent mathematical checks.

Real application data are loaded by optional example scripts from icensBKL and are not redistributed.

R package structure reference: https://cran.r-project.org/doc/manuals/r-release/R-exts.html
