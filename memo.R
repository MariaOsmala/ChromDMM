bin_signal function, how works?

  DirichletMultinomial::dmn


signals_name <- "strand_ambiguous_1000_enhancers_bin_1_window_2000_only5prime.RData"
bin_size=40
save_name="dmn-result-enhancers-withShift-withFlip_40.Rds"
region_type="ranges.signal"
mod_number=9
mod_subset=c("CTCF", "H2AZ", "H3K27ac", "H3K4me1", "H3K4me2", "H3K4me3", "H3K79me2", "H3K9ac", "RNAPOL2")
Kmax=8
eta=2
nu=1
etah=1
nuh=10
shift.reads=TRUE
flip=TRUE
repetition=1
parallel=TRUE

#compile the source code and load the function into R
sourceCpp('src/dmn.cpp')

where is disable_gsl_error_handler() function?

subset.fits <- DirichletMultinomial::dmn(count=signals.subset, K=1:Kmax,
                                         bin.width=bin_size,
                                         verbose = T, eta=eta, nu=nu, etah=etah, nuh=nuh,
                                         shift.reads = shift.reads, flip=flip,
                                         repetition=repetition, parallel=parallel)