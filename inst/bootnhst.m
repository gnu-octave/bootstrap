%  Function File: bootnhst
%
%  Bootstrap null hypothesis significance test(s) (NHST)
%
%  bootnhst(DATA,GROUP,...);
%  bootnhst(DATA,GROUP);
%  bootnhst(DATA,GROUP);
%  bootnhst(...,'bootfun',bootfun);
%  bootnhst(...,'bootfun',{bootfun,bootfun_args});
%  bootnhst(...,'nboot',nboot);
%  bootnhst(...,'ref',ref);
%  bootnhst(...,'block',blocks);    % for bootfun 'mean' or 'robust' only
%  bootnhst(...,'nested',clusters); % for bootfun 'mean' or 'robust' only
%  bootnhst(...,'Options',paropt);
%  bootnhst(...,'alpha',alpha);
%  p = bootnhst(...,'dim',dim)
%  [p,c] = bootnhst(DATA,GROUP,...)
%  [p,c,stats] = bootnhst(DATA,GROUP,...)
%
%  This function uses non-parametric (or semi-parametric) bootstrap for 
%  null hypothesis (H0) significance testing on univariate (vector) or
%  multivatiate (matrix) DATA, to compare bootfun (default is the 'mean') 
%  evaluated on independent or dependent GROUPs (i.e. samples) [1]. This 
%  function is appropriate for post hoc comparisons among a family of 
%  hypothesis tests and comparing groups. The family-wise error rates 
%  (FWER) of pairwise comparisons, or comparisons to a reference group,
%  are controlled by the single-step maxT procedure on bootstrap resamples.
%  Thus, depending on the comparisons requested using the ref input   
%  argument, the p-value adjustments are essentially bootstrap versions of 
%  either Tukey-Kramer's or Dunnett's tests. Unlike these tests though, 
%  bootnhst does not make the normality assumption. Since DATA across the 
%  GROUPs are resampled, as for a permutation test, bootnhst assumes  
%  exchangeability among the groups under the null hypothesis. The sampling 
%  method used for bootstrap and bootknife is balanced, unless computations 
%  are accelerated by parallel processing, in which case bootstrap sampling 
%  (only) is no longer balanced. Note that this function will return an 
%  error if any GROUP label is not represented by more than one data row. 
%
%  bootnhst(DATA,GROUP) performs a k-sample bootstrap test where DATA is 
%  a column vector or matrix, and GROUP is a vector or cell array the same 
%  number of rows as DATA. GROUP contains numbers or strings which denote 
%  GROUP labels. If the number of GROUPs (k) is 2, this is a 2-sample test. 
%  If k > 2, this is test for k GROUPs (like in one-way layout). If GROUP 
%  is numeric and any GROUP is assigned NaN then their respective DATA rows 
%  will be omitted from analysis. Likewise, any rows of DATA containing an
%  NaN will be omitted. If no output arguments are requested, the test 
%  statistics and multiplicity-adjusted p-values for both the overall 
%  hypothesis test, and the post hoc tests for comparison between the 
%  GROUPs, are returned in a pretty table. The differences between GROUPs 
%  are also plot along with the symmetic 100*(1-alpha)% bootstrap-t 
%  confidence intervals (also adjusted to control the FWER). Markers and 
%  error bars are red if p < .05 or blue if p > .05. The default alpha 
%  level is 0.05, which produces 95% confidence intervals.
%
%  bootnhst(...,'bootfun',bootfun) sets the statistic calculated from
%  the bootstrap samples. This can be a function handle, string or cell 
%  array with the function handle or string in the first cell and input 
%  arguments to that function in subsequent cells. The calculation of 
%  bootfun on the DATA must return a scalar value. Note that bootfun MUST 
%  calculate a statistic representative of the finite data sample, it 
%  should NOT be an estimate of a population parameter. For example, for 
%  the variance, set bootfun to {@var,1}, not @var or {@var,0}. The default 
%  value of bootfun is 'mean'.  If empty, the default is @mean or 'mean'. 
%  If DATA is multivariate, bootfun is the grand mean, which is the mean of 
%  the means of each column (i.e. variates). If a robust statistic for 
%  central location is required, setting bootfun to 'robust' implements a 
%  smoothed version of the median (see function help for smoothmedian). 
%  Smooth functions of the data are preferable.
%    Standard errors are estimated by bootknife resampling by default [2], 
%  where nboot(2) corresponds to the number of bootknife resamples. If 
%  nboot(2) is 0 and standard errors are calculated without resampling 
%  (if bootfun is 'mean'), or using leave-one-out jackknife (or cluster-
%  jackknife (if using a 'nested' design). Note that if bootfun is not 
%  the mean, the t-statistics returned by this function will not be 
%  comparable with tabulated values.  
%
%  bootnhst(...,'nboot',nboot) is a vector of upto two positive integers
%  indicating the number of replicate samples for the first (bootstrap) 
%  and second (bootknife) levels of iterated resampling. The default
%  value of nboot is [1000,200]. If a scalar value is provided for nboot,
%  the value will set the number of first level bootstrap samples; the 
%  number of second level bootknife samples will assume the default of 
%  200. Increasing the values of nboot reduces the Monte Carlo error of  
%  the p-value (and confidence interval) estimates but the calculations  
%  take longer to complete. If nboot(2) is explicitly set to 0 (or if a 
%  hierarchical data structure is defined with 'nested') then bootnhst 
%  calculates standard errors for studentization using jackknife (or 
%  cluster-jackknife) resampling instead.
%
%  bootnhst(...,'ref',ref) also sets the GROUP to use as the reference 
%  GROUP for post hoc tests. For a one-way experimental design or family 
%  of tests, this will usually be a control GROUP. If all pairwise 
%  comparisons are desired, then set ref to 'pairwise' or leave empty. 
%  By default, pairwise comparisons are computed for post hoc tests.
%
%  If the ref input argument is empty, the resampling procedure for 
%  pairwise comparisons produces p-value adjustments analagous to the 
%  Tukey-Kramer Honest Significance Difference (HSD) test. 
%
%  If the ref input argument is specified, then the resampling procedure 
%  for treatment versus reference produces p-value adjustments analagous 
%  to Dunnett's post hoc tests (since the range of bootfun in the null 
%  distribution of test statistics is restricted to differences with 
%  the reference GROUP).
%
%  The specification of H0 for the overall hypothesis test depends on whether 
%  a reference GROUP is specified with the ref input argument.
%
%  If ref is empty: 
%    H0: GROUPs of DATA are all sampled from the same population with respect
%        to the parameter defined by bootfun (which by default is the mean). 
%
%  If ref is specified:
%    H0: GROUPs of DATA are all sampled from the same population as DATA in 
%        the GROUP ref with respect to the parameter defined by bootfun (which   
%        by default is the mean).  
%
%  bootnhst(...,'block',blocks) specifies a column vector of numeric 
%  identifiers with the same number of rows as DATA. The identifiers should 
%  indicate block membership of the data rows. Data in a data block are 
%  centered and the resampling is stratified to impose restrictions on the
%  exchangeability of data to within blocks. Since the data must be centered
%  using bootfun, this feature only supports location parameters, of which
%  bootnhst supports bootfun being either 'mean' or 'robust'. This option
%  is appropriate when the family of tests has a randomized block design or
%  one-way repeated measures layout. See end of this help for an example. The
%  'block' option here should not be confused with the block option in ibootci.
%
%  bootnhst(...,'nested',clusters) specifies a column vector of numeric 
%  identifiers with the same number of rows as DATA. The identifiers should 
%  indicate subgroups (a.k.a clusters) of the data rows nested within the 
%  GROUPs. The clusters are resampled by two-stage bootstrap resampling of 
%  residuals with shrinkage correction (see bootstrp help for more 
%  information). For multivariate data, the residuals are calculated as the 
%  difference of the data values from the column-wise mean vector. Because 
%  specifying clusters causes bootnhst to resample residuals, the bootstrap 
%  becomes semi-parametric. The clusters input argument can be used to 
%  accomodate for a single level of nesting in heirarchical data structures. 
%  Since the bootstrap becomes semiparametric, this feature only supports 
%  location parameters, of which bootnhst supports bootfun being either 
%  'mean' or 'robust'. This resampling strategy is appropriate when the 
%  family of tests has a split plot design layout, and is a bootstrap 
%  version of a nested 1-way ANOVA. See end of this help for an example of 
%  this application. This function will return an error if any GROUP is not
%  represented by more than one cluster, but there is no restriction on the 
%  number of data rows in any cluster. Note that the value in nboot(2) is 
%  ignored since specifying cluster identifiers enforces cluster-jackknife 
%  to compute standard errors. If empty, this argument is ignored. 
%
%  bootnhst(...,'Options',paropt) specifies options that govern if and 
%  how to perform bootstrap iterations using multiple processors (if the 
%  Parallel Computing Toolbox or Octave Forge package is available). If 
%  empty, calculations are performed in serial.
%
%  paropt argument is a structure with the following recognised fields:
%
%   'UseParallel' - If true, compute bootstrap iterations in parallel.
%                   Default is false for serial computation. In MATLAB,
%                   the default is true if a parallel pool has already
%                   been started.
%   'nproc'       - The number of processors to use to accelerate 
%                   computations. 
% 
%  bootnhst(...,'alpha',alpha) specifies the two-tailed significance level
%  for confidence interval coverage of 0 (in c).
%
%  bootnhst(...,'dim',dim) specifies which dimension to average over the
%  DATA first when DATA is a matrix. dim can take values of 1 or 2. Note
%  that while setting dim can affect the result when bootfun is the median,
%  both values give the same result when bootfun is the mean (i.e. for the
%  grand mean). This name-value pair is only used if bootfun is 'mean', 
%  'median', 'smoothmedian', or 'robust'.
%
%  p = bootnhst(DATA,GROUP) returns a single p-value for the overall,
%  omnibus hypothesis test and represents the multiplicity-adjusted p-value 
%  for the maximum t-statistic from the set of comparisons (relevant to the 
%  test of the overall null hypothesis, see below). Note that the p-value 
%  returned will be truncated at the resolution limit determined by the 
%  number of bootstrap replicates, specifically 1/nboot(1). 
%
%  [p, c] = bootnhst(DATA,GROUP,...) also returns a 9 column matrix that
%  summarises post hoc test results. The family-wise error rate is 
%  simultaneously controlled since the null distribution for each test 
%  represents the maximum studentized test statistic of the resamples. 
% 
%  The columns of output argument c contain:
%    column 1: reference GROUP number
%    column 2: test GROUP number
%    column 3: value of bootfun evaluated using reference GROUP DATA
%    column 4: value of bootfun evaluated using test GROUP DATA
%    column 5: columns 4 minus column 3
%    column 6: t-ratio
%    column 7: p-adj 
%    column 8: LOWER bound of the bootstrap-t confidence interval
%    column 9: UPPER bound of the bootstrap-t confidence interval
%
%  Note that the p-values returned in column 7 and the length of 
%  the confidence interval limits returned columns 8 and 9 are 
%  corrected/adjusted to control the FWER.
%
%  [p,c,stats] = bootnhst(DATA,GROUP,...) also returns a structure 
%  containing additional statistics. The stats structure contains the 
%  following fields:
%
%   gnames   - group names used in the GROUP input argument. The index of 
%              gnames corresponds to the numbers used to identify GROUPs
%              in columns 1 and 2 of the output argument c
%   ref      - reference group
%   groups   - group number and bootfun for each group with sample size 
%              and standard error.
%   Var      - weighted mean (pooled) sampling variance
%   maxT     - omnibus test statistic (maxT) 
%   df       - degrees of freedom
%   nboot    - number of bootstrap resamples
%   alpha    - two-tailed significance level for the confidence interval 
%              coverage of 0 (in c).
%   blocks   - vector of numeric identifiers indicating block membership.
%   clusters - vector of numeric identifiers indicating cluster membership
%   bootstat - test statistic computed for each bootstrap resample 
%
%  [...] = bootnhst(...,'DisplayOpt',logical) a logical value (true or 
%  false) is used to specify whether to display and graph the results in 
%  addition to creating the output arguments. The default is true.
%
%  Many examples of using bootnhst to analyse data obtained from a 
%  variety of experimental designs:
%
%
%  EXAMPLE 1A: 
%  ONE-WAY ANOVA WITH EQUAL SAMPLE SIZES: Treatment vs. Control (1)
%
%   >> y = [111.39 110.21  89.21  76.64  95.35  90.97  62.78;
%           112.93  60.36  92.29  59.54  98.93  97.03  79.65;
%            85.24 109.63  64.93  75.69  95.28  57.41  75.83;
%           111.96 103.40  75.49  76.69  77.95  93.32  78.70];
%   >> g = [1 2 3 4 5 6 7;
%           1 2 3 4 5 6 7;
%           1 2 3 4 5 6 7;
%           1 2 3 4 5 6 7];
%   >> bootnhst (y(:),g(:),'ref',1,'nboot',5000);
%
% Summary of bootstrap null hypothesis (H0) significance test(s)
% ******************************************************************************
% Overall hypothesis test from single-step maxT procedure
% H0: Groups of data are all sampled from the same population as data in ref
% 
% Maximum t(21) = 3.24, p-adj = .018 
% ------------------------------------------------------------------------------
% 
% POST HOC TESTS with control of the FWER by the single-step maxT procedure
% ------------------------------------------------------------------------------
% | Comparison |  Reference # |       Test # |  Difference |    t(df)|   p-adj |
% |------------|--------------|--------------|-------------|---------|---------|
% |          1 |            1 |            2 |   -9.48e+00 |    0.92 |    .864 |
% |          2 |            1 |            3 |   -2.49e+01 |    2.42 |    .102 |
% |          3 |            1 |            4 |   -3.32e+01 |    3.24 |    .018 |*
% |          4 |            1 |            5 |   -1.35e+01 |    1.31 |    .613 |
% |          5 |            1 |            6 |   -2.07e+01 |    2.01 |    .224 |
% |          6 |            1 |            7 |   -3.11e+01 |    3.03 |    .031 |*
% 
% Where degrees of freedom (df) = 21
% 
% ------------------------------------------------------------------------------
% |    GROUP # |                                                   GROUP label |
% |------------|---------------------------------------------------------------|
% |          1 |                                                             1 |
% |          2 |                                                             2 |
% |          3 |                                                             3 |
% |          4 |                                                             4 |
% |          5 |                                                             5 |
% |          6 |                                                             6 |
% |          7 |                                                             7 |
%
%  EXAMPLE 1B: 
%  ROBUST ONE-WAY ANOVA WITH EQUAL SAMPLE SIZES: Treatment vs. Control (1)
%
%   >> y = [111.39 110.21  89.21  76.64  95.35  90.97  62.78;
%           112.93  60.36  92.29  59.54  98.93  97.03  79.65;
%            85.24 109.63  64.93  75.69  95.28  57.41  75.83;
%           111.96 103.40  75.49  76.69  77.95  93.32  78.70];
%   >> g = [1 2 3 4 5 6 7;
%           1 2 3 4 5 6 7;
%           1 2 3 4 5 6 7;
%           1 2 3 4 5 6 7];
%   >> bootnhst (y(:),g(:),'ref',1,'nboot',5000,'bootfun','robust');
%
% Summary of bootstrap null hypothesis (H0) significance test(s)
% ******************************************************************************
% Overall hypothesis test from single-step maxT procedure
% H0: Groups of data are all sampled from the same population as data in ref
% 
% Maximum t(21) = 2.74, p-adj = .024 
% ------------------------------------------------------------------------------
% 
% POST HOC TESTS with control of the FWER by the single-step maxT procedure
% ------------------------------------------------------------------------------
% | Comparison |  Reference # |       Test # |  Difference |    t(df)|   p-adj |
% |------------|--------------|--------------|-------------|---------|---------|
% |          1 |            1 |            2 |   -6.66e+00 |    0.51 |    .986 |
% |          2 |            1 |            3 |   -3.02e+01 |    2.32 |    .069 |
% |          3 |            1 |            4 |   -3.56e+01 |    2.74 |    .024 |*
% |          4 |            1 |            5 |   -1.62e+01 |    1.25 |    .582 |
% |          5 |            1 |            6 |   -2.00e+01 |    1.54 |    .370 |
% |          6 |            1 |            7 |   -3.49e+01 |    2.68 |    .027 |*
% 
% Where degrees of freedom (df) = 21
% 
% ------------------------------------------------------------------------------
% |    GROUP # |                                                   GROUP label |
% |------------|---------------------------------------------------------------|
% |          1 |                                                             1 |
% |          2 |                                                             2 |
% |          3 |                                                             3 |
% |          4 |                                                             4 |
% |          5 |                                                             5 |
% |          6 |                                                             6 |
% |          7 |                                                             7 |
%
%  EXAMPLE 2A:
%  COMPARISON OF TWO INDEPENDENT GROUPS WITH UNEQUAL SAMPLE SIZES 
%  (analagous to Student's t-test)
%
%   >> y =    [54       43
%              23       34
%              45       65
%              54       77
%              45       46
%             NaN       65];
%   >> g = {'male' 'female'
%           'male' 'female'
%           'male' 'female'
%           'male' 'female'
%           'male' 'female'
%           'male' 'female'};
%   >> bootnhst (y(:),g(:),'ref','male','nboot',5000);
%
% Summary of bootstrap null hypothesis (H0) significance test(s)
% ******************************************************************************
% Overall hypothesis test from single-step maxT procedure
% H0: Groups of data are all sampled from the same population as data in ref
%    
% Maximum t(9) = 1.20, p-adj = .259 
% ------------------------------------------------------------------------------
%
%  EXAMPLE 2B:
%  ONE-WAY ANOVA WITH UNEQUAL SAMPLE SIZES: pairwise comparisons (the 'ref' default)
%
%   >> y = [54  87  45
%           23  98  39
%           45  64  51
%           54  77  49
%           45  89  50
%           47 NaN  55];
%   >> g = [ 1   2   3
%            1   2   3
%            1   2   3 
%            1   2   3
%            1   2   3
%            1   2   3];
%   >> bootnhst (y(:),g(:),'nboot',5000);
%
% Summary of bootstrap null hypothesis (H0) significance test(s)
% ******************************************************************************
% Overall hypothesis test from single-step maxT procedure
% H0: Groups of data are all sampled from the same population
%
% Maximum t(14) = 6.15, p-adj = <.001 
% ------------------------------------------------------------------------------
%
% POST HOC TESTS with control of the FWER by the single-step maxT procedure
% ------------------------------------------------------------------------------
% | Comparison |  Reference # |       Test # |  Difference |    t(df)|   p-adj |
% |------------|--------------|--------------|-------------|---------|---------|
% |          1 |            1 |            2 |   +3.83e+01 |    6.15 |   <.001 |***
% |          2 |            1 |            3 |   +3.50e+00 |    0.59 |    .831 |
% |          3 |            2 |            3 |   -3.48e+01 |    5.59 |   <.001 |***
%
% Where degrees of freedom (df) = 14
%
% ------------------------------------------------------------------------------
% |    GROUP # |                                                   GROUP label |
% |------------|---------------------------------------------------------------|
% |          1 |                                                             1 |
% |          2 |                                                             2 |
% |          3 |                                                             3 |
%
%
%  EXAMPLE 3: 
%  NESTED ONE-WAY ANOVA: pairwise comparisons
%
%   >> y =        [28   32   27   25   26   25   21   19   18
%                  26   27   25   24   28   26   19   18   20
%                  27   28   29   27   29   24   17   23   19
%                  31   29   27   23   27   23   20   20   18];
%   >> g =        [ 1    1    1    2    2    2    3    3    3
%                   1    1    1    2    2    2    3    3    3
%                   1    1    1    2    2    2    3    3    3
%                   1    1    1    2    2    2    3    3    3];
%   >> clusters = [ 1    2    3    4    5    6    7    8    9
%                   1    2    3    4    5    6    7    8    9
%                   1    2    3    4    5    6    7    8    9
%                   1    2    3    4    5    6    7    8    9];
%   >> bootnhst(y(:),g(:),'nested',clusters(:),'nboot',5000);
%
% Summary of bootstrap null hypothesis (H0) significance test(s)
% ******************************************************************************
% Overall hypothesis test from single-step maxT procedure
% H0: Groups of data are all sampled from the same population
%
% Maximum t(6) = 9.01, p-adj = .002 
% ------------------------------------------------------------------------------
%
% POST HOC TESTS with control of the FWER by the single-step maxT procedure
% ------------------------------------------------------------------------------
% | Comparison |  Reference # |       Test # |  Difference |    t(df)|   p-adj |
% |------------|--------------|--------------|-------------|---------|---------|
% |          1 |            1 |            2 |   -2.42e+00 |    2.51 |    .109 |
% |          2 |            1 |            3 |   -8.67e+00 |    9.01 |    .002 |**
% |          3 |            2 |            3 |   -6.25e+00 |    6.50 |    .005 |**
%
% Where degrees of freedom (df) = 6
% 
% ------------------------------------------------------------------------------
% |    GROUP # |                                                   GROUP label |
% |------------|---------------------------------------------------------------|
% |          1 |                                                             1 |
% |          2 |                                                             2 |
% |          3 |                                                             3 |
%
%  EXAMPLE 4A: 
%  COMPARISON OF TWO DEPENDENT GROUPS 
%  (analagous to paired t-test)
%
%   >> y =      [4.5  5.6
%                3.7  6.4
%                5.3  6.4
%                5.4  6.0
%                3.9  5.7];
%   >> g =      [  1    2
%                  1    2
%                  1    2
%                  1    2
%                  1    2];
%   >> blocks = [  1    1
%                  2    2
%                  3    3
%                  4    4
%                  5    5];
%   >> bootnhst (y(:),g(:),'block',blocks(:),'nboot',5000);   
%
% Summary of bootstrap null hypothesis (H0) significance test(s)
% ******************************************************************************
% Overall hypothesis test from single-step maxT procedure
% H0: Groups of data are all sampled from the same population
% 
% Maximum t(4) = 3.91, p-adj = .001 
% ------------------------------------------------------------------------------
%
%  EXAMPLE 4B:
%  ONE-WAY REPEATED MEASURES ANOVA: pairwise comparisons  
%
%   >> y =      [54 43 78 111
%                23 34 37 41
%                45 65 99 78
%                31 33 36 35
%                15 25 30 26];
%   >> g =      [ 1  2  3  4
%                 1  2  3  4
%                 1  2  3  4
%                 1  2  3  4
%                 1  2  3  4];
%   >> blocks = [ 1  1  1  1
%                 2  2  2  2
%                 3  3  3  3
%                 4  4  4  4 
%                 5  5  5  5];
%   >> bootnhst (y(:),g(:),'block',blocks(:),'nboot',5000);
%
% Summary of bootstrap null hypothesis (H0) significance test(s)
% ******************************************************************************
% Overall hypothesis test from single-step maxT procedure
% H0: Groups of data are all sampled from the same population
% 
% Maximum t(12) = 2.83, p-adj = .003 
% ------------------------------------------------------------------------------
% 
% POST HOC TESTS with control of the FWER by the single-step maxT procedure
% ------------------------------------------------------------------------------
% | Comparison |  Reference # |       Test # |  Difference |    t(df)|   p-adj |
% |------------|--------------|--------------|-------------|---------|---------|
% |          1 |            1 |            2 |   +6.40e+00 |    0.74 |    .846 |
% |          2 |            1 |            3 |   +2.24e+01 |    2.58 |    .010 |**
% |          3 |            1 |            4 |   +2.46e+01 |    2.83 |    .003 |**
% |          4 |            2 |            3 |   +1.60e+01 |    1.84 |    .163 |
% |          5 |            2 |            4 |   +1.82e+01 |    2.10 |    .086 |
% |          6 |            3 |            4 |   +2.20e+00 |    0.25 |    .993 |
% 
% Where degrees of freedom (df) = 12
% 
% ------------------------------------------------------------------------------
% |    GROUP # |                                                   GROUP label |
% |------------|---------------------------------------------------------------|
% |          1 |                                                             1 |
% |          2 |                                                             2 |
% |          3 |                                                             3 |
% |          4 |                                                             4 |
%
% 
%  EXAMPLE 5: 
%  MANOVA or TWO-WAY REPEATED MEASURES ANOVA: pairwise comparisons
%
%   >> y = [34  54
%           65  91
%           12  13
%           35  29
%           55  79
%           99 121];
%   >> g = [1
%           1
%           2
%           2
%           3
%           3];
%   >> bootnhst (y,g,'nboot',5000);
%
% Summary of bootstrap null hypothesis (H0) significance test(s)
% ******************************************************************************
% Overall hypothesis test from single-step maxT procedure
% H0: Groups of data are all sampled from the same population
%
% Maximum t(3) = 2.77, p-adj = .190 
% ------------------------------------------------------------------------------
%
% POST HOC TESTS with control of the FWER by the single-step maxT procedure
% ------------------------------------------------------------------------------
% | Comparison |  Reference # |       Test # |  Difference |    t(df)|   p-adj |
% |------------|--------------|--------------|-------------|---------|---------|
% |          1 |            1 |            2 |   -3.88e+01 |    1.62 |    .445 |
% |          2 |            1 |            3 |   +2.75e+01 |    1.15 |    .624 |
% |          3 |            2 |            3 |   +6.62e+01 |    2.77 |    .190 |
%
% Where degrees of freedom (df) = 3
%
% ------------------------------------------------------------------------------
% |    GROUP # |                                                   GROUP label |
% |------------|---------------------------------------------------------------|
% |          1 |                                                             1 |
% |          2 |                                                             2 |
% |          3 |                                                             3 |
%
%  EXAMPLE 6A: 
%  MANOVA or TWO-WAY REPEATED MEASURES ANOVA: comparing 2 groups
%   >> y = [34    35    78    54    42
%           65    67   111    98    89
%           39    41   167   143   136
%           65    54   211   178   146];
%   >> g = [ 1
%            1
%            2
%            2];
%   >> bootnhst (y,g,'nboot',5000);
%
% Summary of bootstrap null hypothesis (H0) significance test(s)
% ******************************************************************************
% Overall hypothesis test from single-step maxT procedure
% H0: Groups of data are all sampled from the same population
% 
% Maximum t(2) = 2.51, p-adj = .376 
% ------------------------------------------------------------------------------
% 
%
%  EXAMPLE 6B: 
%  TWO-WAY REPEATED MEASURES ANOVA: SIMPLE EFFECTS
%   >> y =      [ 34    35    78    54    42
%                 65    67   111    98    89
%                 39    41   167   143   136
%                 65    54   211   178   146];
%   >> g =      [  1     2     3     4     5
%                  1     2     3     4     5 
%                NaN   NaN   NaN   NaN   NaN
%                NaN   NaN   NaN   NaN   NaN];
%    >> blocks = [  1     1     1     1     1
%                   2     2     2     2     2
%                   3     3     3     3     3
%                   4     4     4     4     4];
%   >> bootnhst (y(:),g(:),'block',blocks(:),'ref',1,'nboot',5000);
%
% Summary of bootstrap null hypothesis (H0) significance test(s)
% ******************************************************************************
% Overall hypothesis test from single-step maxT procedure
% H0: Groups of data are all sampled from the same population as data in ref
% 
% Maximum t(4) = 8.48, p-adj = .001 
% ------------------------------------------------------------------------------
% 
% POST HOC TESTS with control of the FWER by the single-step maxT procedure
% ------------------------------------------------------------------------------
% | Comparison |  Reference # |       Test # |  Difference |    t(df)|   p-adj |
% |------------|--------------|--------------|-------------|---------|---------|
% |          1 |            1 |            2 |   +1.50e+00 |    0.28 |    .988 |
% |          2 |            1 |            3 |   +4.50e+01 |    8.48 |    .001 |**
% |          3 |            1 |            4 |   +2.65e+01 |    4.99 |    .008 |**
% |          4 |            1 |            5 |   +1.60e+01 |    3.02 |    .056 |
% 
% Where degrees of freedom (df) = 4
% 
% ------------------------------------------------------------------------------
% |    GROUP # |                                                   GROUP label |
% |------------|---------------------------------------------------------------|
% |          1 |                                                             1 |
% |          2 |                                                             2 |
% |          3 |                                                             3 |
% |          4 |                                                             4 |
% |          5 |                                                             5 |
% 
%   >> g =      [NaN   NaN   NaN   NaN   NaN
%                NaN   NaN   NaN   NaN   NaN
%                  1     2     3     4     5
%                  1     2     3     4     5 ];
%   >> bootnhst (y(:),g(:),'block',blocks(:),'ref',1,'nboot',5000);
%
% Summary of bootstrap null hypothesis (H0) significance test(s)
% ******************************************************************************
% Overall hypothesis test from single-step maxT procedure
% H0: Groups of data are all sampled from the same population as data in ref
% 
% Maximum t(4) = 13.46, p-adj = .002 
% ------------------------------------------------------------------------------
% 
% POST HOC TESTS with control of the FWER by the single-step maxT procedure
% ------------------------------------------------------------------------------
% | Comparison |  Reference # |       Test # |  Difference |    t(df)|   p-adj |
% |------------|--------------|--------------|-------------|---------|---------|
% |          1 |            1 |            2 |   -4.50e+00 |    0.44 |    .956 |
% |          2 |            1 |            3 |   +1.37e+02 |   13.46 |    .002 |**
% |          3 |            1 |            4 |   +1.08e+02 |   10.66 |    .004 |**
% |          4 |            1 |            5 |   +8.90e+01 |    8.74 |    .006 |**
% 
% Where degrees of freedom (df) = 4
% 
% ------------------------------------------------------------------------------
% |    GROUP # |                                                   GROUP label |
% |------------|---------------------------------------------------------------|
% |          1 |                                                             1 |
% |          2 |                                                             2 |
% |          3 |                                                             3 |
% |          4 |                                                             4 |
% |          5 |                                                             5 |
%
% Note that the simple (main) effects analyses here will only pool the variance   
% within each level of the main effect (not across all levels of the main effect)
%
%
%  EXAMPLE 7: 
%  TWO-WAY ANOVA: SIMPLE EFFECTS
%
%   >> y = [  34  36  41 NaN  43  98  87  95  99  88
%             23  19  26  29  25  32  29  26  33 30];
%   >> g = [   1   1   1   1   1   2   2   2   2   2 
%            NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN]; 
%   >> bootnhst (y(:),g(:),'nboot',5000);
%
% Summary of bootstrap null hypothesis (H0) significance test(s)
% ******************************************************************************
% Overall hypothesis test from single-step maxT procedure
% H0: Groups of data are all sampled from the same population
%
% Maximum t(7) = 16.15, p-adj = .002  
% ------------------------------------------------------------------------------
%    
%   >> g = [ NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN
%              1   1   1   1   1   2   2   2   2   2];
%   >> bootnhst (y(:),g(:),'nboot',5000);
%    
% Summary of bootstrap null hypothesis (H0) significance test(s)
% ******************************************************************************
% Overall hypothesis test from single-step maxT procedure
% H0: Groups of data are all sampled from the same population
%
% Maximum t(8) = 2.69, p-adj = .031 
% ------------------------------------------------------------------------------
%
%
%   Bibliography:
%   [1] Efron and Tibshirani. Chapter 16 Hypothesis testing with the
%        bootstrap in An introduction to the bootstrap (CRC Press, 1994)
%   [2] Hesterberg, Tim C. (2004), Unbiasing the Bootstrap - Bootknife-
%        Sampling vs. Smoothing, Proceedings of the Section on Statistics 
%        and the Environment, American Statistical Association, 2924-2930.
%
%  bootnhst v2.0.0.0 (28/06/2022)
%  Author: Andrew Charles Penn
%  https://www.researchgate.net/profile/Andrew_Penn/
%
%  Copyright 2019 Andrew Charles Penn
%  This program is free software: you can redistribute it and/or modify
%  it under the terms of the GNU General Public License as published by
%  the Free Software Foundation, either version 3 of the License, or
%  (at your option) any later version.
%
%  This program is distributed in the hope that it will be useful,
%  but WITHOUT ANY WARRANTY; without even the implied warranty of
%  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%  GNU General Public License for more details.
%
%  You should have received a copy of the GNU General Public License
%  along with this program.  If not, see <http://www.gnu.org/licenses/>.


function [p, c, stats] = bootnhst (data, group, varargin)

  % Check if running in Octave (else assume Matlab)
  info = ver; 
  isoctave = any (ismember ({info.Name}, 'Octave'));
  
  % Apply defaults
  bootfun = 'mean';
  nboot = [1000,200];
  ref = [];
  alpha = 0.05;
  strata = [];
  clusters = [];
  dim = 1;
  DisplayOpt = true;
  paropt = struct;
  paropt.UseParallel = false;
  if isoctave
    paropt.nproc = nproc;
  else
    paropt.nproc = feature('numcores');
  end

  % Fetch extra input arguments
  argin3 = varargin;
  narg = numel(argin3);
  if narg > 1
    while ischar(argin3{end-1})
      if strcmpi(argin3{end-1},'bootfun')
        bootfun = argin3{end};
      elseif strcmpi(argin3{end-1},'nboot')
        nboot = argin3{end};
      elseif strcmpi(argin3{end-1},'ref')
        ref = argin3{end};
      elseif any(strcmpi({'Options','Option'},argin3{end-1}))
        paropt = argin3{end};
      elseif strcmpi(argin3{end-1},'alpha')
        alpha = argin3{end};
      elseif strcmpi(argin3{end-1},'dim')
        dim = argin3{end};
      elseif strcmpi(argin3{end-1},'DisplayOpt')
        DisplayOpt = argin3{end};
      elseif any(strcmpi(argin3{end-1},{'cluster','clusters','clustered','nested','nests','nesting'}))
        clusters = argin3{end};
      elseif any(strcmpi(argin3{end-1},{'block','blocks','stratum','strata','stratified'}))
        strata = argin3{end};
      else
        error('unrecognised input argument to bootnhst')
      end
      argin3 = {argin3{1:end-2}};
      narg = numel(argin3);
      if narg < 1
        break
      end
    end
  end

  % Error checking
  % Check and process bootnhst input arguments
  nvar = size(data,2);
  if (nargin < 2)
    error('bootnhst requires atleast two input arguments');
  end
  if ischar(group)
    group = cellstr(group);
  end
  if (size(group,1)>1) && (size(data,1) ~= size(group,1))
    error('DATA and GROUP must have the same number of rows')
  end
  if iscell(group)
    if ~iscellstr(group)
      group = cell2mat(group);
    end
  end
  if ~isa(nboot,'numeric')
    error('nboot must be numeric');
  end
  if any(nboot~=abs(fix(nboot)))
    error('nboot must contain positive integers')
  end
  if numel(nboot) > 2
    error('the vector nboot cannot have length > 2')
  elseif numel(nboot) < 2
    % If nboot is scalar, set default number of bootknife samples to calculate standard errors
    % nboot(2) must be explcitly set to 0 to request jackknife standard errors
    nboot = cat(2, nboot, 200);
  end
  if nboot(1) < 1000
    error('the minimum allowable value of nboot(1) is 1000')
  end 
  if iscell(bootfun)
    func = bootfun{1};
    args = bootfun(2:end);
    bootfun = @(data) feval(func, data, args{:});
  end
  if isa(bootfun,'function_handle')
    if strcmp (func2str (bootfun), 'mean')
      if nvar > 1
        % Grand mean for multivariate data
        bootfun = @(data) mean(mean(data,dim));
      else
        bootfun = 'mean';
      end
    elseif strcmp (func2str (bootfun), 'smoothmedian')
      if nvar > 1 
        % Grand smoothed median for multivariate data
        bootfun = @(data) smoothmedian(smoothmedian(data,dim));
      else
        bootfun = 'smoothmedian';
      end
    else
      if ~isempty(strata)
        error('bootfun must be ''mean'' or ''robust'' for block or repeated measures designs.')
      end
      if ~isempty(clusters)
        error('bootfun must be ''mean'' for nested designs.')
      end
    end
  elseif isa(bootfun,'char')
    if strcmpi(bootfun,'mean') 
      if nvar > 1
        % Grand mean for multivariate data
        bootfun = @(data) mean(mean(data,dim));
      else
        bootfun = 'mean';
      end
    elseif any(strcmpi(bootfun,{'robust','smoothmedian'}))
      if nvar > 1
        % Grand smoothed median for multivariate data
        bootfun = @(data) smoothmedian(smoothmedian(data,dim));
      else
        bootfun = 'smoothmedian';
      end
    else
      if ~isempty(strata)
        error('bootfun must be ''mean'' or ''robust'' for block or repeated measures designs.')
      end
      if ~isempty(clusters)
        error('bootfun must be ''mean'' or ''robust'' for nested designs.')
      end
    end
  end
  
  % Error checking
  if ~isempty(ref) && strcmpi(ref,'pairwise')
    ref = [];
  end
  if ~isa(dim,'numeric')
    error('dim must be numeric');
  end
  if (dim ~= 1) && (dim ~= 2)
    error('dim must be either 1 or 2');
  end
  if ~isempty(clusters)
    if (size(clusters,2) > 1) || (size(clusters,1) ~= size(data,1))
      error('clusters must be a column vector with the same number of rows as DATA')
    end
    opt = struct;
  end
  if ~isempty(strata)
    if (size(strata,2) > 1) || (size(strata,1) ~= size(data,1))
      error('strata must be a column vector with the same number of rows as DATA')
    end
  end
  if ~isempty(strata) && ~isempty(clusters)
    error('block and nested options cannot be used together')
  end
  if nargout > 3
    error('bootnhst only supports up to 3 output arguments')
  end

  % Data or group exclusion using NaN 
  if isnumeric(group)
    if any(isnan(group))
      data(isnan(group),:) = [];
      if ~isempty(clusters)
        clusters(isnan(group)) = [];
      end
      if ~isempty(strata)
        strata(isnan(group)) = [];
      end
      group(isnan(group)) = [];
    end
  end
  if any(any(isnan([data]),2))
    group(any(isnan([data]),2)) = [];
    if ~isempty(clusters)
      clusters(any(isnan([data]),2)) = [];
    end
    if ~isempty(strata)
      strata(any(isnan([data]),2)) = [];
    end
    data(any(isnan([data]),2),:) = [];
  end

  % Assign non-zero numbers to group labels
  [gnames,junk,g] = unique(group);
  gk = unique(g);
  k = numel(gk);
  if ~isempty(ref)
    if isnumeric(ref)
      ref = gk(ismember(gnames,ref));
    else
      ref = gk(strcmp(gnames,ref));
    end
  end
            
  % Get data structure information
  if isempty(clusters)
    N = size(g,1);
  else
    N = numel(unique(clusters));
  end
  if isempty(strata)
    l = 1;
  else
    sid = unique (strata);      % strata ID
    l = numel (sid);            % number of strata
  end
 
  % If applicable, center each stratum on it's respective (grand) mean or smoothed median 
  % Note that the bootstrap becomes semiparametric
  if ~isempty(strata)
    S = zeros(N, l, 'logical');
    for i=1:l
      % Create strata matrix
      S(:,i) = ismember(strata, sid(i));   % strata logical indexing
      data(S(:,i),:) = data(S(:,i),:) - feval(bootfun, feval(bootfun,data(S(:,i),:),2) ,1);
    end
    ns = sum(S);                           % number of data elements per stratum
  end
  
  % If applicable, setup a parallel pool 
  if ~isoctave
    % MATLAB
    if paropt.UseParallel 
      % PARALLEL
      if (paropt.nproc > 0) 
        % MANUAL
        try 
          pool = gcp ('nocreate'); 
          if isempty (pool)
            if (paropt.nproc > 1)
              % Start parallel pool with nproc workers
              pool = parpool (paropt.nproc);
            else
              % Parallel pool is not running and nproc is 1 so run function evaluations in serial
              paropt.UseParallel = false;
            end
          else
            if (pool.NumWorkers ~= paropt.nproc)
              % Check if number of workers matches nproc and correct it accordingly if not
              delete (pool);
              if (paropt.nproc > 1)
                parpool (paropt.nproc);
              end
            end
          end
        catch
          % Parallel toolbox not installed, run function evaluations in serial
          paropt.nproc = 1;
        end
      else
        % AUTOMATIC
        try 
          pool = gcp ('nocreate'); 
          if isempty (pool)
            % Parallel pool not running, start parallel pool using all available workers
            parpool;
          else
            % Parallel pool is already running, set nproc to the number of workers
            paropt.nproc = pool.NumWorkers;
          end
        catch
          % Parallel toolbox not installed, run function evaluations in serial
          paropt.UseParallel = false;
        end
      end
    end
  end

  % Define a function to calculate maxT
  func = @(data) maxstat (data, g, nboot(2), bootfun, ref, clusters, strata, isoctave);

  % Perform resampling and calculate bootstrap statistics
  if isempty(clusters)
    % Use newer, faster and balanced (less biased) resampling function (boot)
    if ~isempty (strata)
      bootsam = zeros (N, nboot(1), 'int16');
      for i = 1:l
        bootsam(S(:,i),:) = boot (ns(i), nboot(1), false);
        rows = find (S(:,i));
        bootsam(S(:,i),:) = rows(bootsam(S(:,i), :));
      end
    else
      bootsam = boot (N, nboot(1), false);
    end
    if isoctave
      % OCTAVE
      if paropt.UseParallel
        % Set unique random seed for each parallel thread
        pararrayfun(paropt.nproc, @boot, 1, 1, false, 1, 1:paropt.nproc);
        % Evaluate maxstat on each bootstrap resample in PARALLEL 
        cellfunc = @(bootsam) feval (func, data (bootsam, :));
        Q = parcellfun(paropt.nproc, cellfunc, num2cell (bootsam, 1), 'ChunksPerProc', 100);
      else
        % Evaluate maxstat on each bootstrap resample in SERIAL
        cellfunc = @(bootsam) feval (func, data (bootsam, :));
        Q = cellfun (cellfunc, num2cell (bootsam, 1));
      end
    else
      % MATLAB
      if paropt.UseParallel
        % Set unique random seed for each parallel thread
        parfor i = 1:paropt.nproc; boot (1, 1, false, 1, i); end;
        % Evaluate maxstat on each bootstrap resample in PARALLEL 
        Q = zeros (1, nboot(1));
        parfor h = 1:nboot(1)
          Q(h) = feval (func, data (bootsam (:, h), :));
        end
      else
        % Evaluate maxstat on each bootstrap resample in SERIAL
        cellfunc = @(bootsam) feval (func, data (bootsam, :));
        Q = cellfun (cellfunc, num2cell (bootsam, 1));
      end
    end
  else
    % Use legacy bootstrp function for two-stage nonparametric bootstrap sampling 
    % with shrinkage correction for clustered data
    state = warning; 
    warning off;    % silence warnings about non-vectorized bootfun
    Q = bootstrp (nboot(1),func,data,'cluster',clusters,'Options',paropt);
    warning(state);
  end

  % Compute the estimate (theta) and it's pooled (weighted mean) sampling variance 
  theta = zeros(k,1);
  SE = zeros(k,1);
  Var = zeros(k,1);
  t = zeros(nboot(2),1);
  nk = zeros(size(gk));
  for j = 1:k
    if ~isempty(clusters)
      theta(j) = feval(bootfun,data(g==gk(j),:));
      % Compute unbiased estimate of the standard error by
      % cluster-jackknife resampling
      opt.clusters = clusters(g==gk(j));
      nk(j) = numel(unique(opt.clusters));
      SE(j) = jack(data(g==gk(j),:), bootfun, [], opt);
    elseif (nboot(2) == 0)
      nk(j) = sum(g==gk(j));
      if strcmp (bootfun, 'mean')
        theta(j) = mean(data(g==gk(j),:));
        % Quick calculation for the standard error of the mean
        SE(j) = std(data(g==gk(j),:),0) / sqrt(nk(j));
      else
        theta(j) = feval(bootfun,data(g==gk(j),:));
        % If requested, compute unbiased estimates of the standard error using jackknife resampling
        SE(j) = jack(data(g==gk(j),:), bootfun);
      end
    else
      % Compute unbiased estimate of the standard error by balanced bootknife resampling
      % Bootknife resampling involves less computation than Jackknife when sample sizes get larger
      theta(j) = feval(bootfun,data(g==gk(j),:));
      nk(j) = sum(g==gk(j));
      stats = bootknife(data(g==gk(j),:),[nboot(2),0],bootfun,[],[],0,isoctave);
      SE(j) = stats(3);
    end
    Var(j) = ((nk(j)-1)/(N-k-(l-1))) * SE(j)^2;
  end
  if any(nk <= 1)
    error('the number of observations or clusters per group must be greater than 1')
  end
  nk_bar = sum(nk.^2)./sum(nk);  % weighted mean sample size
  Var = sum(Var.*nk/nk_bar);     % pooled sampling variance weighted by sample size
  df = sum(nk)-k-(l-1);          % degrees of freedom

  % Calculate weights to correct for unequal sample size  
  % when calculating standard error of the difference
  w = nk_bar./nk;

  % Prepare to make symmetrical bootstrap-t confidence intervals
  [cdf,QS] = empcdf(Q,0);

  % Compute resolution limit of the p-values as determined by resampling with nboot(1) resamples
  res = 1/nboot(1);

  % Calculate p-values for comparisons adjusted to simultaneously control the FWER
  if isempty(ref)
    % Single-step maxT procedure for pairwise comparisons is a resampling version of Tukey-Kramer HSD test
    A = ones(k,1)*gk';
    B = tril(gk*ones(1,k),-1);
    M = [A(:) B(:)];
    ridx = (M(:,2)==0);
    M(ridx,:)=[];
    n = size(M,1);
    c = zeros(n,9);
    c(:,1:2) = M;
    for i = 1:n
      c(i,3) = theta(c(i,1));
      c(i,4) = theta(c(i,2));
      c(i,5) = c(i,4) - c(i,3);
      SED = sqrt(Var * (w(c(i,1)) + w(c(i,2))));
      c(i,6) = abs(c(i,5)) / SED;
      if (c(i,6) < QS(1))
        c(i,7) = interp1(QS,1-cdf,c(i,6),'linear',1);
      else
        c(i,7) = interp1(QS,1-cdf,c(i,6),'linear',res);
      end
      c(i,8) = c(i,5) - SED * interp1(cdf,QS,1-alpha,'linear');
      c(i,9) = c(i,5) + SED * interp1(cdf,QS,1-alpha,'linear');
    end
  else
    % Single-step maxT procedure for treatment vs control comparisons is a resampling version of Dunnett's test
    c = zeros(k,9);
    c(:,1) = ref;
    c(:,3) = theta(ref);
    for j = 1:k
      c(j,2) = gk(j);
      c(j,4) = theta(c(j,2));
      c(j,5) = c(j,4) - c(j,3); 
      SED = sqrt(Var * (w(c(j,1)) + w(c(j,2))));
      c(j,6) = abs(c(j,5)) / SED;
      if (c(j,6) < QS(1))
        c(j,7) = interp1(QS,1-cdf,c(j,6),'linear',1);
      else
        c(j,7) = interp1(QS,1-cdf,c(j,6),'linear',res);
      end
      c(j,8) = c(j,5) - SED * interp1(cdf,QS,1-alpha,'linear');
      c(j,9) = c(j,5) + SED * interp1(cdf,QS,1-alpha,'linear');
    end
    c(ref,:) = [];
  end

  % Calculate (maximum) test statistic and (minimum) p-value for the omnibus test
  maxT = max(c(:,6));
  p = min(c(:,7));

  % Prepare stats output structure
  stats = struct;
  stats.gnames = gnames;
  stats.ref = ref;
  stats.groups = zeros(k,5);
  stats.groups = zeros(k,5);
  stats.groups(:,1) = gk;
  stats.groups(:,2) = theta;
  stats.groups(:,3) = nk;
  stats.groups(:,4) = SE;
  %stats.groups(:,5) = theta - sqrt((0.5*(w+1)).*Var/2) * interp1(cdf,QS,1-alpha,'linear');
  %stats.groups(:,6) = theta + sqrt((0.5*(w+1)).*Var/2) * interp1(cdf,QS,1-alpha,'linear');
  stats.Var = Var;
  stats.maxT = maxT;
  stats.df = df;
  stats.nboot = nboot;
  stats.alpha = alpha;
  stats.blocks = strata;
  stats.clusters = clusters;
  stats.bootstat = Q;

  % Print output and plot graph with confidence intervals if no output arguments are requested
  cols = [1,2,5,6,7]; % columns in c that we want to print data for
  if (nargout == 0) || (DisplayOpt == true)
    fprintf (['\n',...
                    'Summary of bootstrap null hypothesis (H0) significance test(s)\n',...
                    '******************************************************************************\n']);
    if isempty(ref)
      fprintf (['Overall hypothesis test from single-step maxT procedure\n',...
               'H0: Groups of data are all sampled from the same population\n\n']);
    else
      fprintf (['Overall hypothesis test from single-step maxT procedure\n',...
               'H0: Groups of data are all sampled from the same population as data in ref\n\n']);
    end
    if (p <= 0.001)
      fprintf (['Maximum t(%u) = %.2f, p-adj = <.001 \n',...
                '------------------------------------------------------------------------------\n'],[df,maxT]);
    elseif (p > 0.999)
      fprintf (['Maximum t(%u) = %.2f, p-adj = 1.000 \n',...
          '------------------------------------------------------------------------------\n'],[df,maxT]);
    else
      fprintf (['Maximum t(%u) = %.2f, p-adj = .%03u \n',...
                '------------------------------------------------------------------------------\n'],[df,maxT,round(p*1000)]);
    end
    if size(c,1) > 1
      fprintf (['\n',...
                'POST HOC TESTS with control of the FWER by the single-step maxT procedure\n',...
                '------------------------------------------------------------------------------\n',...
                '| Comparison |  Reference # |       Test # |  Difference |    t(df)|   p-adj |\n',...
                '|------------|--------------|--------------|-------------|---------|---------|\n']);
      if isempty(ref)
        for i = 1:n
          tmp = num2cell(c(i,cols));
          tmp{end} = round(tmp{end} * 1000);
          if (c(i,7) <= 0.001)
            tmp(end) = [];
            fprintf ('| %10u | %12u | %12u | %+11.2e | %7.2f |   <.001 |***\n',i,tmp{:});
          elseif (c(i,7) > 0.999)
            tmp(end) = [];
            fprintf ('| %10u | %12u | %12u | %+11.2e | %7.2f |   1.000 |\n',i,tmp{:});
          else
            fprintf ('| %10u | %12u | %12u | %+11.2e | %7.2f |    .%03u |',i,tmp{:});
            if c(i,7) < 0.01
              fprintf('**\n')
            elseif c(i,7) < 0.05
              fprintf('*\n')
            else
              fprintf('\n')
            end
          end
        end
      else
        for j = 1:k-1
          tmp = num2cell(c(j,cols));
          tmp{end} = round(tmp{end} * 1000);
          if (c(j,7) <= 0.001)
            tmp(end) = [];
            fprintf ('| %10u | %12u | %12u | %+11.2e | %7.2f |   <.001 |***\n',j,tmp{:});
          elseif (c(j,7) > 0.999)
            tmp(end) = [];
            fprintf ('| %10u | %12u | %12u | %+11.2e | %7.2f |   1.000 |\n',j,tmp{:});
          else
            fprintf ('| %10u | %12u | %12u | %+11.2e | %7.2f |    .%03u |',j,tmp{:});
            if c(j,7) < 0.01
              fprintf('**\n')
            elseif c(j,7) < 0.05
              fprintf('*\n')
            else
              fprintf('\n')
            end
          end
        end
      end
      fprintf('\nWhere degrees of freedom (df) = %u\n',df)
      fprintf (['\n',...
                '------------------------------------------------------------------------------\n',...
                '|    GROUP # |                                                   GROUP label |\n',...
                '|------------|---------------------------------------------------------------|\n']);
      if ~iscellstr(gnames)
        gnames = cellstr(num2str(gnames));
      end
      for j = 1:k
        fprintf ('| %10u | %61s |\n',gk(j),gnames{j});
      end
      fprintf ('\n')
    end

    % Plot graph of the difference in bootfun for each comparison with 100*(1-alpha)% confidence intervals
    figure;
    nc = size(c,1);                               % Calculate number of comparisons to plot
    plot([0;0],[0;nc+1]','k:');                   % Plot vertical dashed line at 0 effect
    set(gca,'Ydir','reverse')                     % Flip y-axis direction
    ylim([0.5,nc+0.5]);                           % Set y-axis limits
    hold on                                       % Enable plotting new data on the same axis
    for i=1:nc
      if c(i,7) < 0.05
        plot(c(i,5),i,'or','MarkerFaceColor','r') % Plot marker for the difference in bootfun 
        plot([c(i,8),c(i,9)],i*ones(2,1),'r-')    % Plot line for each confidence interval 
      else
        plot(c(i,5),i,'ob','MarkerFaceColor','b') % Plot marker for the difference in bootfun 
        plot([c(i,8),c(i,9)],i*ones(2,1),'b-')    % Plot line for each confidence interval 
      end
    end
    hold off
    xlabel(sprintf('%g%% bootstrap-t confidence interval for the difference',100*(1-alpha)));
    ylabel('Comparison number (Test - Reference)');   

  end

end

