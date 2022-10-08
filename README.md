# statistics-bootstrap package

## Package maintainer
Andrew Penn (andy.c.penn@gmail.com)

## Package contributors
Andrew Penn

## A statistics package for Octave/Matlab providing a variety of bootstrap resampling tools

This package of functions can be used to estimate uncertainty (confidence intervals) and test hypotheses (*p*-values) using bootstrap. Variations of the bootstrap are included that improve the accuracy of bootstrap statistics for small samples.

## Requirements and dependencies

This package is known to be compatible with versions of Octave v3.6.0+ and Matlab v7.4.0+. Some features of this package have specific dependencies:

 * The optional jackknife functionality in `bootnhst` requires the Statistics package (in Octave) or the Statistics and Machine Learning Toolbox (in Matlab).  
 * The `bootanovan` function requires either the Statistics package version v1.5+ (in Octave v6.1+), or the Statistics and Machine Learning Toolbox (in Matlab)
 * The `bootcoeff` and `bootemm` functions require Octave v6.1+ and the Statistics package version v1.5+. 
 * All parallel computing options require either the parallel package (in Octave) or the Parallel Computing Matlab Toolbox (in Matlab).

## Installation
 
To install (or test) the statistics-bootstrap package at it's existing location in either Octave or Matlab, follow these steps: 
 
 * Download the package. If it is a compressed file, decompress it.
 * Open Octave or Matlab command prompt.
 * `cd` to the package directory. (The directory contains a file called 'make.m' and 'install.m')
 * Type `make` to compile the mex files from source (or use the precompiled binaries if available).  
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
* `bootknife` performs balanced bootknife resampling and calculates bootstrap bias, standard error and confidence intervals. The interval types supported are simple percentile, bias-corrected and accelerated, or calibrated. This function supports iterated and stratified resampling.
* `bootnhst` calculates *p*-values by bootstrap null-hypothesis significance testing (two-tailed). This function can be used to compare 2 or more (independent) samples. This function resamples under the null hypothesis.
* `bootmode` uses bootstrap to evaluate the likely number of real modes in a distribution
* `bootci` is a function for calculating bootstrap confidence intervals. This function is a wrapper of the `bootknife` function but has the same usage as the `bootci` function from Matlab's Statistics and Machine Learning toolbox.  
* `bootstrp` is a function for calculating bootstrap statistics. This function is a wrapper of the `bootknife` function but has the same usage as the `bootstrp` function from Matlab's Statistics and Machine Learning toolbox.  
* `bootcoeff` (Octave only) is a function for semi-parametric bootstrap of the regression coefficients from a linear model fit using `anovan` or `fitlm`.  
* `bootemm` (Octave only) is a function for semi-parametric bootstrap of the estimated marginal means from a linear model fit using `anovan` or `fitlm`.  

At the Octave command prompt, type `help function-name` for more information about the function and it's usage.

## Development roadmap
 
* Create more documentation and guidance for using the functions in this package  
* Provide the option in bootanovan to print a pretty ANOVA table of the results  
