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

  Base package in Suggests/Enhances imported in NAMESPACE:
    'graphics'

> `graphics` is no longer imported. Now, it's just a Suggests.

Changes to worse in reverse depends:

Package: BayesBrainMap
Check: examples
New result: ERROR

  Found the following significant warnings:
    Note: possible error in 'dual_reg(BOLD, prior$mean, ': argument 4 matches multiple formal arguments

> This error is expected: `BayesBrainMap` calls `fMRItools` functions whose arguments have changed in this release. An updated version of `BayesBrainMap`, compatible with the new `fMRItools`, will be submitted to CRAN shortly after this submission is accepted.