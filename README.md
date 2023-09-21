# statistics-bootstrap package

## Package maintainer
Andrew Penn (andy.c.penn@gmail.com)

## Package contributors
Andrew Penn

## A statistics package for Octave/Matlab providing a variety of bootstrap resampling tools

This package of functions can be used to estimate bias, uncertainty (standard errors and confidence intervals), prediction error, and test hypotheses (*p*-values) using bootstrap resampling. Variations of the bootstrap are included that improve the coverage and accuracy of bootstrap confidence intervals for small samples and samples with complex dependence structures.  

## Requirements and dependencies

Core functions in this package are known to be compatible with versions of Octave v4.4.0+ and Matlab v7.4.0+. Some features of this package have specific dependencies:

 * All parallel computing options require either the parallel package (in Octave) or the Parallel Computing Matlab Toolbox (in Matlab).  
 * The optional jackknife functionality in `bootnhst` requires the Statistics package (in Octave) or the Statistics and Machine Learning Toolbox (in Matlab).  
 
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
* `bootstrp` is a function for calculating bootstrap statistics. This function has the same usage as the `bootstrp` function from Matlab's Statistics and Machine Learning toolbox.  
* `bootci` is a function for calculating bootstrap confidence intervals. This function has the same usage as the `bootci` function from Matlab's Statistics and Machine Learning toolbox. Unlike Matlab's `bootci`, this function also provides the option for obtaining calibrated bootstrap confidence intervals by iterated bootstrap.
* `bootknife` performs balanced bootknife resampling and calculates bootstrap bias, standard error and confidence intervals for any user defined-function of the data. The interval types supported are simple percentile, bias-corrected and accelerated (BCa), or calibrated percentile. This function supports iterated and stratified resampling, but not cluster resampling (see `bootclust`).
* `bootclust` can perform balanced bootstrap or bootknife resampling of clustered data and calculate bootstrap bias, standard error and confidence intervals for any user defined-function of the data. The interval types supported are simple percentile, and bias-corrected and accelerated (BCa).
* `bootbayes` performs Bayesian nonparametric bootstrap and calculates posterior statistics for the mean, or regression coefficients from a linear model. Two credible interval types are supported: shortest probability intervals and percentile intervals. Cluster and block resampling is supported in cases where the distribution of residuals and errors in the model have a dependence structure. See also `bootlm`.
* `bootwild` performs wild bootstrap-t resampling and calculates confidence intervals and frequentist *p*-values for the mean, or regression coefficients from a linear model (H0 = 0). Cluster and block resampling is supported in cases where the distribution of residuals and errors in the model have a dependence structure. See also `bootlm`.
* `bootlm` is a bootstrap function for calculating confidence intervals and *p*-values for the regression coefficients from a linear model, estimated marginal means or posthoc tests. The usage of this function is similar to `anovan` and can also return bootstrapped ANOVA statistics and prediction errors. This function uses `bootwild` or `bootbayes` to improve robustness and reduce bias in the presence of heteroscedasticity and it supports cluster/block resampling (e.g. for nested multifactorial designs).
* `bootnhst` calculates *p*-values by bootstrap null-hypothesis significance testing (two-tailed) in simple designs comparing 2 or more (independent) samples in designs with a one-way layout. This function resamples under the null hypothesis (assuming exchangeability similar to a permutation test) and can be used to compare functions of the data other than the mean (e.g. robust estimators). The function also computes and returns multiple comparison tests that control the family-wise error rate.
* `randtest2` is a function for performing permutation or randomization tests to compare independent or paired distributions using the Wassertstein metric. A feature of this function is that one can define the sampling units for clustered resampling in the case of nested experimental designs. 
* `bootmode` uses bootstrap to evaluate the likely number of real modes in a distribution


At the Octave/MATLAB command prompt, type `help function-name` for more information about the function and it's input and output arguments. In Octave, you can also request demonstrations of function usage through examples by typing 'demo function-name` at the command prompt.

## Development roadmap

