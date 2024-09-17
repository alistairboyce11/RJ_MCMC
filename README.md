# RJ_MCMC
This program performs a simultaneous Transdimensional inversion of Surface wave dispersion curves with the reversible jump algorithm

Primary contributors:

Thomas Bodin, Dorian Soergel, Alistair Boyce

README last updated: A, Boyce 09/12/2021

----------------------------- CITATION -----------------------------

If you use this code please cite the following manuscript:

Bodin, T., Leiva, J., Romanowicz, B., Maupin, V. & Yuan, H. (2016) Imaging anisotropic layering with Bayesian inversion of multiple data types. Geophys. J. Int. 206, 605–629. doi: 10.1093/gji/ggw124.
  
----------------------------- OUTLINE  -----------------------------

1. LICENSE
2. README.md
3. data_joint.h                    : parameter file for joint inversion, input data
4. params.h                        : parameter file for internal variables
5. Makefile                        : Fortran compilation file (uses mpif90)
6. Model_PREM_SIMPLE.in            : Reference model (Currently PREM)
7. RJ_MCMC_joint.f90               : MAIN CODE
8. postprocess_binary_outputs.f90  : Postprocessing of bayesian output
9. src/                            : Supporting code & functions
10. dispersion_all.in              : input dispersion curve file
11. process_MCMC_output_joint.py   : performs the re-weighting step
12. postprocess_util.py            : supporting functions for re-weighting
13. postprocess_functions.py       : functions to be applied on each sample
14. plot_results_joint.py          : Plotting scripts in python

----------------------------- DETAILS   -----------------------------

The main code is compiled using `make joint` . The relevant parameters are mostly in data_joint.h, in particular the input/output file names, the number of temperatures and the number of iterations, thinning and burn-in, as well as the bounds for inverted parameters, in particular the bounds for Vs, Xi, Vp and the uncertainty parameters for Rayleigh and Love waves. The bounds for the number of layers is in params.h . 

The code produces unformatted output files, which contain a header containing different parameters of the inversion and for each model the number of layers of the model, the values of the different parameters for each layer as well as the corresponding dispersion curves. These binary files are then converted to text files using the `postprocess_binary_outputs.f90` code compiled with `make postprocess_binary`. 

The re-weighting step, that for each dispersion curve and each sample of the bayesian inversion determines a weight and then applies it to various functions is done by `process_MCMC_output_joint.py`. The functions to be applied on each sample are defined in `postprocess_functions.py`. This produces output hdf5 files. 

The plotting can then be done using `plot_results_joint.py`. 
