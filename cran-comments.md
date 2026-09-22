## Test environments

* Mac aarch64-apple-darwin23, macOS Tahoe 26.5.2, R 4.6.0

## R CMD check results

0 errors ✔ | 0 warnings ✔ | 0 notes ✔

## Downstream dependencies

`BayesfMRI`, `fMRIscrub`, `templateICAr`, and `BayesBrainMap` do not suffer 
  additional warnings or errors with this new version of `fMRItools`.

## Tests

Passes all the tests in `tests/run_fMRItools_tests.R`

## Previous submission

> We've included additional updates since our most recent submission on Sep 8.

Changes to worse in reverse depends:

Package: BayesBrainMap
Check: examples
New result: ERROR

  Found the following significant warnings:
    Note: possible error in 'dual_reg(BOLD, prior$mean, ': argument 4 matches multiple formal arguments

> This error is expected because the arguments to functions imported by `BayesBrainMap` have changed. A new version of `BayesBrainMap` compatible with the new version of `fMRItools` has been submitted to CRAN.