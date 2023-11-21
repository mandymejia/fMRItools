## Test environments

* Windows x86_64-w64-mingw32/x64, R 4.2.2
* Mac x86_64-apple-darwin17.0, R 4.3.1

## R CMD check results

0 errors ✔ | 0 warnings ✔ | 0 notes ✔

## Downstream dependencies

BayesfMRI 0.3.5 passes checks.
fMRIscrub 0.14.5 passes checks.
templateICAr 0.6.2 passes checks.

## Tests

Passes all the tests in `tests/run_fMRItools_tests.R`

## Previous submission

0.4.1 Failed an automatic check because of brackets without '\code' in roxygen
of `make_mask`. This has been fixed.