## Test environments

* Mac x86_64-apple-darwin17.0, R 4.4.0

## R CMD check results

0 errors ✔ | 0 warnings ✔ | 0 notes ✔

## Downstream dependencies

BayesfMRI 0.3.11 passes checks.
fMRIscrub 0.14.5 passes checks.
templateICAr 0.9.1 passes checks.

## Tests

Passes all the tests in `tests/run_fMRItools_tests.R`

## Previous submission

  Mismatches for apparent methods not registered:
  image:
    function(x, ...)
  image.scale:
    function(z, zlim, col, breaks, axis.pos, add.axis, ...)
  See section 'Registering S3 methods' in the 'Writing R Extensions'
  manual.

This function has been renamed to `image_scale`.