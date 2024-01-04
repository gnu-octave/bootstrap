## statistics-resampling package

### Package maintainer
Andrew Penn (andy.c.penn@gmail.com)

### Package contributors
Andrew Penn

### Citations
If you use this package, please include the following citation(s):

* Penn, Andrew Charles. (2020). Resampling methods for small samples or samples with complex dependence structures. *Zenodo*. [https://doi.org/10.5281/zenodo.3992392](https://doi.org/10.5281/zenodo.3992392)  

### Description

The statistics-resampling package is an Octave package and Matlab toolbox that can be used to overcome a wide variety of statistics problems using non-parametric resampling methods. In particular, the functions included can be used to estimate bias, uncertainty (standard errors and confidence intervals), prediction error, and test hypotheses (*p*-values). Variations of the resampling methods are included that improve the accuracy of the statistics for small samples and samples with complex dependence structures. 

### Requirements and dependencies

Core functions in this package are known to be compatible with versions of Octave v4.4.0+ and Matlab v7.4.0+. Some features of this package have specific dependencies:

 * All parallel computing options require either the parallel package (in Octave) or the Parallel Computing Matlab Toolbox (in Matlab).  
 * The optional jackknife functionality in `boot1way` requires the Statistics package (in Octave) or the Statistics and Machine Learning Toolbox (in Matlab).  
 
### Installation
 
To install (or test) the statistics-resampling package at it's existing location in either Octave or Matlab, follow these steps: 
 
 * Download the package. If it is a compressed file (.zip or .tar.gz), extract it's contents and move the package directory to the desired location.
 * Open Octave or Matlab (command prompt).
 * Change directory (cd) into the package folder. (The directory contains a file called 'make.m' and 'install.m', among others)
 * Type `make` to compile the mex files from source (or use the precompiled binaries if available. If suitable precompiled binaries are not available for your platform, then Matlab/Octave will need access to a C++11 compiler. Note that if you skip the make step, then the package functions will still work, but some will run slower. This step is interactive so check the command window.) 
 * Type `install`. The package will load now (and automatically in the future) when you start Octave/Matlab.
 
 If you want or need to uninstall the package/toolbox, change directory (cd) into the package folder and type uninstall.
 
 Alternatively, users of more recent versions of Octave can install the package automatically with the following command:
 
 `pkg install "https://github.com/gnu-octave/statistics-resampling/archive/refs/heads/master.zip"`
 
 The package can be loaded on demand in Octave with the following commmand:
 
 `pkg load statistics-resampling`
 
 In Octave, you can find out basic information about the package by typing: `pkg describe -verbose statistics-resampling`  

### Usage

All help and demos are documented on the 'Function Reference' page in the [manual](https://gnu-octave.github.io/statistics-resampling/). If you do not see the navigation pane on the manual web pages, please enable javascript in your browser.

Function help can also be requested from the Octave/MATLAB command prompt, by typing `help function-name`. 

In Octave only, you can also request the demonstrations of function usage by typing `demo function-name` at the command prompt.  
