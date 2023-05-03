# statistics-bootstrap package

## Package maintainer
Andrew Penn (andy.c.penn@gmail.com)

## Package contributors
Andrew Penn

## A statistics package for Octave/Matlab providing a variety of bootstrap resampling tools

This package of functions can be used to estimate bias, uncertainty (standard errors and confidence intervals) and test hypotheses (*p*-values) using bootstrap resampling. Variations of the bootstrap are included that improve the accuracy of bootstrap statistics for small samples.

## Requirements and dependencies

Core functions in this package are known to be compatible with versions of Octave v3.8.2+ and Matlab v7.4.0+. Some features of this package have specific dependencies:

 * All parallel computing options require either the parallel package (in Octave) or the Parallel Computing Matlab Toolbox (in Matlab).  
 * The optional jackknife functionality in `bootnhst` requires the Statistics package (in Octave) or the Statistics and Machine Learning Toolbox (in Matlab).  
 * The `bootcoeff` and `bootemm` functions require the Statistics package version v1.5+ and therefore also Octave v6.1.0+.  
 
## Installation
 
To install (or test) the statistics-bootstrap package at it's existing location in either Octave or Matlab, follow these steps: 
 
 * Download the package. If it is a compressed file, decompress it.
 * Open Octave or Matlab command prompt.
 * `cd` to the package directory. (The directory contains a file called 'make.m' and 'install.m')
 * Type `make` to compile the mex files from source (or use the precompiled binaries if available. If precompiled binaries are not available for your platform, then Matlab/Octave will need access to a C++ compiler. Note that if you skip this step, the package functions will still work, but will run significantly slower.) 
 * Type `install`. The package will load now (and automatically in the future) when you start Octave/Matlab.
 
 To uninstall, `cd` to the package directory and type  `uninstall`.
 
 Alternatively, users of more recent versions of Octave can install the package automatically with the following command:
 
 `pkg install "https://github.com/gnu-octave/statistics-bootstrap/archive/refs/heads/master.zip"`
 
 The package can be loaded on demand in Octave with the following commmand:
 
 `pkg load statistics-bootstrap`
 
 In Octave, you can find out basic information about the package by typing: `pkg describe -verbose statistics-bootstrap`  

## Usage

### Functions

* `boot` returns resamples data or indices created by balanced bootstrap or bootknife resampling  
* `bootknife` performs balanced bootknife resampling and calculates bootstrap bias, standard error and confidence intervals. The interval types supported are simple percentile, bias-corrected and accelerated, or calibrated percentile. This function supports iterated and stratified resampling.
* `bootbayes` performs Bayesian bootstrap and calculates bootstrap bias, standard error and confidence intervals for coefficients from a linear regression model. The interval types supported are simple percentile.
* `bootnhst` calculates *p*-values by bootstrap null-hypothesis significance testing (two-tailed). This function can be used to compare 2 or more (independent) samples in designs with a one-way layout. This function resamples under the null hypothesis.
* `bootmode` uses bootstrap to evaluate the likely number of real modes in a distribution
* `bootci` is a function for calculating bootstrap confidence intervals. This function is a wrapper of the `bootknife` function but has the same usage as the `bootci` function from Matlab's Statistics and Machine Learning toolbox.  
* `bootstrp` is a function for calculating bootstrap statistics. This function is a wrapper of the `bootknife` function but has the same usage as the `bootstrp` function from Matlab's Statistics and Machine Learning toolbox.  
* `bootcoeff` (Octave only) is a function for calculating bootstrap confidence intervals for the regression coefficients of a linear model fit using `anovan` or `fitlm`. This function uses `bootbayes`.
* `bootemm` (Octave only) is a function for calculating bootstrap confidence intervals for the estimated marginal means of a linear model fit using `anovan` or `fitlm`. This function uses `bootbayes`.

At the Octave/MATLAB command prompt, type `help function-name` for more information about the function and it's input and output arguments. In Octave, you can also request demonstrations of function usage through examples by typing 'demo function-name` at the command prompt.

## Development roadmap
 
* We intend on including the following features in version 5.2.0: 1) an option in `cor` function to select Spearman's rank or Pearson's correlation coefficient; and 2) an additional function folder for fitting linear models and returning bootstrap confidence intervals for regression coefficients and estimated marginal means. The function have similar usage to `anovan`.  
* Make `bootcoeff` and `bootemm` compatible with Matlab.  

