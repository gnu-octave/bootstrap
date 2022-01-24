# statistics-bootstrap package

## Package maintainer
Andrew Penn (andy.c.penn@gmail.com)

## Package contributors
Andrew Penn

## Citations
If you use this package please include the following citation(s):

* Penn, Andrew Charles. (2020). iboot: Iterated bootstrap for small samples or samples with complex dependence structures. Zenodo. https://doi.org/10.5281/zenodo.3992392  

## A statistics package for Octave/Matlab providing a variety of bootstrap resampling tools

This package of functions can be used to estimate uncertainty (confidence intervals) and test hypotheses (*p*-values) using bootstrap. Variations of the bootstrap are included that improve the accuracy of bootstrap statistics for small samples and samples with complex dependence structures. The Octave statistics-bootstrap package is forked from the GitHub repository [iboot](https://github.com/acp29/iboot).

## Requirements and dependencies

This package is known to be compatible with versions of Octave v3.6.0+ and Matlab v7.4.0+. 

Most features of of the package do not depend on other packages. An exception is the parallel computing options, which require either the Octave forge parallel package or the Parallel Computing Matlab Toolbox.

## Installation
 
To install (or test) the statistics-bootstrap package at it's existing location in either Octave or Matlab, follow these steps: 
 
 * Download the package. If it is a compressed file, decompress it.
 * Open Octave or Matlab command prompt.
 * `cd` to the package directory. (The directory contains a file called 'install.m')
 * Type `install`. The package will load now and automatically when you start Octave/Matlab
 
 To uninstall, follow the above steps but use the `uninstall` command
 
 Alternatively, users of more recent versions of Octave can install the package automatically with the following command:
 
 `pkg install "https://github.com/gnu-octave/statistics-bootstrap/archive/refs/heads/master.zip"`
 
 The package can be loaded on demand in Octave with the following commmand:
 
 `pkg load statistics-bootstrap`
 
 In Octave, you can find out basic information about the package by typing: `pkg describe -verbose statistics-bootstrap`  

## Usage

### Functions

* `bootstrp` performs bootstrap resampling 
* `bootci` calculates confidence intervals using bootstrap resampling
* `ibootci` calculates confidence intervals (calibrated) by iterated bootstrap resampling 
* `ibootp` calculates a two-tailed *p*-value for hypothesized value of the statistic using bootstrap
* `bootnhst` calculates *p*-values by bootstrap null-hypothesis significance testing (two-tailed). This function can be used to compare 2 or more (independent) samples.
* `bootanovan` calculates *p*-values for N-way, fixed effects ANOVA by bootstrapping the distribution of F-statistics under the null hypothesis. This function depends on the `anovan` function from the Statistics package in Octave, or the Statistics and Machine Learning toolbox in Matlab.
* `iboottest` is a convenience function that uses `ibootci` and `ibootp` to compute confidence intervals and *p*-values for the difference between two paired samples or between one sample and a population value (two-tailed)
* `iboottest2` is a convenience function that uses `ibootci` and `ibootp` to compute confidence intervals and *p*-values for the difference between two independent (i.e. unpaired) samples. 
* `bootmode` uses bootstrap to evaluate the likely number of real modes in a distribution
* `plotboot` plots an overlay of a histogram, kernel density estimate and interval limits from bootstrap statistics

At the Octave command prompt, type `help function-name` for more information about the function and it's usage.

### Notes 

The Matlab Statistics and Machine Learning toolbox has functions also called `bootstrp` and `bootci`. The same-named functions in the statistics-bootstrap package have almost identical usage to the Matlab functions from the Statistics and Machine Learning toolbox. Be aware though that Matlab's own `bootci` function has a couple of errors, namely in the calculation of the bias for `cper` and `bca` intervals, and in the calculation of `stud` intervals. We recommend using the statistics-bootstrap package function `bootci`, or better it's function `ibootci` for the calculation of calibrated bootstrap confidence intervals. The 'i' prefix to the bootstrap function name indicates that the functions have capabilities to calculate calibrated bootstrap statistics confidence intervals or *p*-values by iterated (a.k.a. double) bootstrap.

## Development roadmap
 
* Create a simple and efficient progress bar for the functions in this package  
* Add function to param folder to calculate D statistic from KS-test (with input format suitable for `bootnhst`, `iboottest` and (`i`)`bootci`) 
* Add the ability to perform stratified jackknife resampling in `jack` 
* Allow bootstrap to accept empty argument for `bootfun` (to enhance compatible behaviours with corresponding matlab function) 
* Create usage guide with examples of how to apply the package functions to different problems and scenarios  
* Create some compiled mex and oct files on various platforms (mac, pc and linux) to speed up computations. 
* Include optional output argument in `bootanovan` that includes confidence intervals for the model coefficients
 
