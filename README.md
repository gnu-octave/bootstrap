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

The statistics-resampling package is an Octave package and Matlab toolbox that can be used to perform a wide variety of statistics tasks using non-parametric resampling methods. In particular, the functions included can be used to estimate bias, uncertainty (standard errors and confidence intervals), prediction error, and calculate *p*-values for null hypothesis significance tests. Variations of the resampling methods are included that improve the accuracy of the statistics for small samples and samples with complex dependence structures. 

### Requirements and dependencies

Core functions in this package are known to be compatible with versions of Octave 4.4.0+ and Matlab R2007a 7.4.0+. Some features of this package have specific dependencies:

 * All parallel computing options require either the parallel package (in Octave) or the Parallel Computing Matlab Toolbox (in Matlab).  
 * The optional jackknife functionality in `boot1way` requires the Statistics package (in Octave) or the Statistics and Machine Learning Toolbox (in Matlab).  
 
  
### Running the statistics-resampling package online

You can try out the statistics-resampling package online: 
 
 [![Open in MATLAB Online](https://www.mathworks.com/images/responsive/global/open-in-matlab-online.svg)](https://matlab.mathworks.com/open/github/v1?repo=gnu-octave/statistics-resampling&file=README.md)  
 
  [![Binder](https://mybinder.org/badge_logo.svg)](https://mybinder.org/v2/gh/acpennlab/statistics-resampling/master?labpath=statistics-resampling.ipynb)  

  Note that the first time you load use statistics-resampling online with Binder it will take a few minutes to build a docker image.  
 
### Installation
 
To install (or test) the statistics-resampling package at a location of your choice for either **Matlab or Octave**, follow these steps: 
 
 * Download the latest package release from [here](https://github.com/gnu-octave/statistics-resampling/releases/). Extract (not just browse) the contents of the compressed file (.zip or .tar.gz), and move the package directory to the desired location.
 * Open Octave or Matlab (command prompt).
 * Change directory (cd) into the package folder. (The directory contains a file called 'make.m' and 'install.m', among others)
 * Type `make` to compile the MEX files from source (or use the precompiled binaries if available. If suitable precompiled binaries are not available for your platform, then Matlab/Octave will need access to a C++11 compiler. Note that if you skip the make step, then the package functions will still work, but some will run slower. This step is interactive so check the command window.) 
 * Type `install`. The package will load now (and automatically in the future) when you start Octave/Matlab.
 
If you want or need to uninstall the package, change directory (cd) into the package folder and type uninstall.
 
Alternatively, **Octave** users can install the latest release of the package just like any other Octave package by typing:

 `pkg install -forge statistics-resampling`
 
Or for the most recent developmental version of the package:
 
 `pkg install "https://github.com/gnu-octave/statistics-resampling/archive/refs/heads/master.zip"`
 
The package can then be loaded on demand in Octave with the following commmand:
 
 `pkg load statistics-resampling`  
 
  (Note that this isn't necessary if you used the local installation instructions first described in this section)
   
Alternatively, **MATLAB** users can conveniently install the package functions as a toolbox by double-clicking the 'statistics-resampling.mltbx' file in the matlab subdirectory. The toolbox installed in this way can be disabled or uninstalled via MATLAB's Add-On manager. MEX files are included with the toolbox installation for Windows (32- or 64-bit), MacOS (Intel or Apple Silicon 64-bit) and Linux (64-bit). 
  
### Usage

All help and demos are documented on the ['Function Reference'](https://gnu-octave.github.io/statistics-resampling/function_reference) page in the [manual](https://gnu-octave.github.io/statistics-resampling/). If you do not see the navigation pane on the manual web pages, please enable javascript in your browser. If you need further help with using any of the functions in this package, please post your questions on the [discussions page](https://github.com/gnu-octave/statistics-resampling/discussions).  
  
Function help can also be requested directly from the Octave/MATLAB command prompt, by typing `help function-name` - substituting in the actual function name.
  
In Octave only, you can get a basic overview of the package and it's functions by typing: `pkg describe -verbose statistics-resampling`, or request demonstrations of function usage by typing `demo function-name`. Users can also request help with using functions and programming in Octave at the [discourse group](https://octave.discourse.group/c/help/6).  

*TIP: Integrate use of the statistics-resampling package into your analysis workflow along with other programming languages (e.g. Python, R or Julia etc.) using [Jupyter Lab with an Octave kernel](https://jupyter.org/try).*  

### Quick start

Below are links to demonstrations of how the bootstrap or randomization functions from this package can be used to perform variants of some commonly used statistical tests, but without the Normality assumption:  
   
 * [t-test to compare two independent samples (assuming equal variances)](https://gnu-octave.github.io/statistics-resampling/function/boot1way.html#1)  

 * [t-test to compare two independent samples (assuming equal variances)](https://gnu-octave.github.io/statistics-resampling/function/boot1way.html#2) (but also robust to outliers)  

 * [t-test to compare two independent samples](https://gnu-octave.github.io/statistics-resampling/function/bootlm.html#1) (but also robust to heteroscedasticity, i.e. unequal variance)  

 * Nested t-test to compare two independent samples. (See example listed below for nested one-way ANOVA)   
 
 * [t-test to compare two paired or matching samples](https://gnu-octave.github.io/statistics-resampling/function/bootlm.html#2)

 * [One-way ANOVA to compare two or more independent samples](https://gnu-octave.github.io/statistics-resampling/function/bootlm.html#3) (but also robust to heteroscedasticity)  

 * [Nested one-way ANOVA to compare two or more independent samples](https://gnu-octave.github.io/statistics-resampling/function/bootlm.html#13) (but also robust to heteroscedasticity and grouping of observations)    
 
 * [Balanced two-way factorial ANOVA](https://gnu-octave.github.io/statistics-resampling/function/bootlm.html#6) (but also robust to heteroscedasticity)  

 * [Unbalanced two-way factorial ANOVA](https://gnu-octave.github.io/statistics-resampling/function/bootlm.html#7) (but also robust to heteroscedasticity)  
 
 * [Simple linear regression](https://gnu-octave.github.io/statistics-resampling/function/bootlm.html#8) (but also robust to heteroscedasticity)  
 
 * [One-way ANCOVA](https://gnu-octave.github.io/statistics-resampling/function/bootlm.html#9) (but also robust to heteroscedasticity)  
 
 * [One-way repeated measures ANOVA](https://gnu-octave.github.io/statistics-resampling/function/bootlm.html#4)  
 
 * [Randomization or permutation tests to compare the distributions of two independent or paired/matching samples](https://gnu-octave.github.io/statistics-resampling/function/randtest2.html#1)  
 
 * [Statistically evaluate the number of modes (i.e. peaks) in a distribution](https://gnu-octave.github.io/statistics-resampling/function/bootmode.html#1)  
 
 * [Correcting sample size calculations for a two-sample test with nested design](https://gnu-octave.github.io/statistics-resampling/function/sampszcalc.html#6)  
  
### Issues

If you find bugs or have any suggestions, please raise an issue on GitHub [here](https://github.com/gnu-octave/statistics-resampling/issues). Please, make sure that when reporting a bug you provide as much information as possible for other users to be able to replicate it.   
  
Package: [statistics-resampling](https://gnu-octave.github.io/statistics-resampling/)  

