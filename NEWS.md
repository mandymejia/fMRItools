# 0.4.0

* Add `all_binary`.

# 0.3.2

* Add workaround for when PESEL estimates zero components.

# 0.3.1

* Add back `dim_reduce` for now, but warn of its deprecation (moved to `templateICAr`).
* Add draft of first step of `harmonize`

# 0.3.0

* Do not transpose matrix at the end of `scale_timeseries`.
* More robust PCA, to handle co-linear columns and try another routine if fail

# 0.2.0

* Export more functions