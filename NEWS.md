# 0.9.0

* Refactor aCompCor into noise-PC calculation and data-cleaning steps. 
* Add `prelim_nuisance` support before PC calculation.
* Replace `CompCor_HCP()` and `CompCor_HCP_cii()` with `CompCor_HCP_calc()` and `CompCor_HCP_apply()`. Rename `CompCor.noise_comps()` to `noise_PCs()` and standardize the noise-PC result format.
* Correct `fsl_bptf()` and use efficient matrix multiplication to apply the filter.
* Revamp `despike_3D()`.
* Preserve the intercept when scaling design matrices, and fix edge cases in
  `scale_timeseries()` and `all_binary()`.
* More efficient CompCor and nuisance regression code.

# 0.8.0

* `norm_BOLD` overhaul: big regression; update scaling.
* Replace `dice_overlap` with `dice_coef`. Draft `match_nets`.

# 0.7.0

* Remove `dim_reduce`
* Fixes and updates to `plot_FC_gg`
* Add `dice_overlap`
* Center BOLD data for `dual_reg_parc` 

# 0.6.0

* Update `UT2mat`. Default: symmatric matrix. 
* Add `plot_FC_gg` from `templateICAr`
* enable `lines="all"` in `plot_FC`.

# 0.5.0

* Add `expand_RPs`
* Add highpass filter to `fsl_bptf`


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