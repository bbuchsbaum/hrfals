simulate_simple_data <- function(ncond = 2,
                                 nreps = 12,
                                 TR = 1,
                                 snr = 1,
                                 hrf = fmrireg::HRF_SPMG3,
                                 seed = 123) {
  fmrireg::simulate_simple_dataset(ncond = ncond,
                                   nreps = nreps,
                                   TR = TR,
                                   snr = snr,
                                   hrf = hrf,
                                   seed = seed)
}
