## Read me

### Package maintainer
Andrew Penn (andy.c.penn@gmail.com)

### Package contributors
Andrew Penn  
(More contributors are welcome!)

### Citations
If you use this package, please include the following citation(s):

* Penn, Andrew Charles. (2020). Resampling methods for small samples or samples with complex dependence structures. *Zenodo*. [https://doi.org/10.5281/zenodo.3992392](https://doi.org/10.5281/zenodo.3992392) 

(Note that package versions 5.4.3 and below were named the 'statistics-bootstrap' package. The 'statistics-resampling' package is a more developed version of the older ['iboot'](https://github.com/acp29/iboot) package). 

### Description

The statistics-resampling package is an Octave package and Matlab toolbox that can be used to overcome a wide variety of statistics problems using non-parametric resampling methods. In particular, the functions included can be used to estimate bias, uncertainty (standard errors and confidence intervals), prediction error, and test hypotheses (*p*-values). Variations of the resampling methods are included that improve the accuracy of the statistics for small samples and samples with complex dependence structures. 

### Requirements and dependencies

Core functions in this package are known to be compatible with versions of Octave 4.4.0+ and Matlab R2007a 7.4.0+. Some features of this package have specific dependencies:

 * All parallel computing options require either the parallel package (in Octave) or the Parallel Computing Matlab Toolbox (in Matlab).  
 * The optional jackknife functionality in `boot1way` requires the Statistics package (in Octave) or the Statistics and Machine Learning Toolbox (in Matlab).  
 
### Installation
 
To install (or test) the statistics-resampling package at a location of your choice, follow these steps: 
 
 * Download the latest package release from [here](https://github.com/gnu-octave/statistics-resampling/releases/). Extract (not just browse) the contents of the compressed file (.zip or .tar.gz), and move the package directory to the desired location.
 * Open Octave or Matlab (command prompt).
 * Change directory (cd) into the package folder. (The directory contains a file called 'make.m' and 'install.m', among others)
 * Type `make` to compile the MEX files from source (or use the precompiled binaries if available. If suitable precompiled binaries are not available for your platform, then Matlab/Octave will need access to a C++11 compiler. Note that if you skip the make step, then the package functions will still work, but some will run slower. This step is interactive so check the command window.) 
 * Type `install`. The package will load now (and automatically in the future) when you start Octave/Matlab.
 
 If you want or need to uninstall the package/toolbox, change directory (cd) into the package folder and type uninstall.
 
 Alternatively, users of more recent versions of Octave can install the package automatically with the following command:
 
 `pkg install "https://github.com/gnu-octave/statistics-resampling/archive/refs/heads/master.zip"`
 
 The package can be loaded on demand in Octave with the following commmand:
 
 `pkg load statistics-resampling`
 
 MATLAB users can conveniently install the package functions as a toolbox by double-clicking the 'statistics-resampling.mltbx' file in the matlab subdirectory. The toolbox installed in this way can be disabled or uninstalled via MATLAB's Add-On manager. Currently, MEX files are included with the toolbox installation for Windows (32- or 64-bit), MacOS (Intel 64-bit) and Linux (64-bit). Without the MEX files, all functionality of the package is available, but some of the functions run slower.  
 
 N.B. The package does not yet include any MEX files (for Octave or Matlab) precompiled for macOS with Apple silicon processors, since the package developers do not have access to this computer platform. If you used this package on macOS with an Apple silicon processor (M1-3 chip), please consider contacting the package maintainer to contribute the MEX files to this project.  
  
### Usage

All help and demos are documented on the ['Function Reference'](https://gnu-octave.github.io/statistics-resampling/function_reference) page in the [manual](https://gnu-octave.github.io/statistics-resampling/). If you do not see the navigation pane on the manual web pages, please enable javascript in your browser. If you need help with using any of the functions in this package, please post your questions on the [discussions page](https://github.com/gnu-octave/statistics-resampling/discussions).  
  
Function help can also be requested directly from the Octave/MATLAB command prompt, by typing `help function-name` - substituting in the actual function name.
  
In Octave only, you can get a basic overview of the package and it's functions by typing: `pkg describe -verbose statistics-resampling`, or request demonstrations of function usage by typing `demo function-name`. Users can also request help with using functions and programming in Octave at the [discourse group](https://octave.discourse.group/c/help/6).  

### Quick start

Below are some links to demonstrations of how the bootstrap or randomization functions from this package can be used to perform variants of the commonly used statistical tests below, but without the Normality assumption:  
  
 * [t-test to compare two independent samples (assuming equal variances)](https://gnu-octave.github.io/statistics-resampling/function/boot1way.html#1)  

 * [t-test to compare two independent samples (assuming equal variances but robust to outliers)](https://gnu-octave.github.io/statistics-resampling/function/boot1way.html#2)  

 * [t-test to compare two independent samples (not assuming equal variances)](https://gnu-octave.github.io/statistics-resampling/function/bootlm.html#1)  

 * [t-test to compare two paired or matching samples (also robust to heteroscedasticity)](https://gnu-octave.github.io/statistics-resampling/function/bootlm.html#2)  

 * [One-way ANOVA to compare more than two independent samples (also robust to heteroscedasticity)](https://gnu-octave.github.io/statistics-resampling/function/bootlm.html#3)  

 * [Nested one-way ANOVA to compare more than two independent samples (also robust to heteroscedasticity and grouping of observations)](https://gnu-octave.github.io/statistics-resampling/function/bootlm.html#11)  
 
 * [Two-way factorial ANOVA (also robust to heteroscedasticity)](https://gnu-octave.github.io/statistics-resampling/function/bootlm.html#6)  

 * [One-way ANCOVA (also robust to heteroscedasticity)](https://gnu-octave.github.io/statistics-resampling/function/bootlm.html#7)  
  
### Issues

If you find bugs or have any suggestions, please raise an issue on GitHub [here](https://github.com/gnu-octave/statistics-resampling/issues). Please, make sure that when reporting a bug you provide as much information as possible for other users to be able to replicate it.   
  
Package: [statistics-resampling](https://gnu-octave.github.io/statistics-resampling/)  

