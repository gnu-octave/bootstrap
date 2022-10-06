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
%  bootnhst(...,'Options',paropt);
%  bootnhst(...,'alpha',alpha);
%  p = bootnhst(...,'dim',dim)
%  [p,c] = bootnhst(DATA,GROUP,...)
%  [p,c,stats] = bootnhst(DATA,GROUP,...)
%
%  This function uses non-parametric bootstrap for null hypothesis (H0) 
%  significance testing on univariate (vector) or multivatiate (matrix) 
%  DATA, to compare bootfun (default is the 'mean') evaluated on 
%  independent GROUPs (i.e. samples) of data in a one-way layout [1].  
%  This function is appropriate for post hoc comparisons among a family  
%  of hypothesis tests and comparing groups. The family-wise error rates 
%  (FWER) of pairwise comparisons, or comparisons to a reference group,
%  are controlled by the single-step maxT procedure on bootstrap resamples.
%  Thus, depending on the comparisons requested using the ref input   
%  argument, the p-value adjustments are essentially bootstrap versions of 
%  either Tukey-Kramer's or Dunnett's tests. Unlike these tests though, 
%  bootnhst does not make the normality assumption. Since DATA across the 
%  GROUPs are resampled, as for a permutation test, bootnhst assumes  
%  exchangeability among the groups under the null hypothesis. The sampling 
%  method used is balanced bootknife resampling. Note that this function 
%  will return an error if any GROUP label is not represented by more than 
%  one data row. 
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
%  error bars are red if p < .05, or blue if p > .05. The default alpha 
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
%  value of bootfun is 'mean'. If empty, the default is @mean or 'mean'. 
%  If a robust statistic for central location is required, setting bootfun
%  to 'robust' implements a smoothed version of the median (see function help
%  for smoothmedian). Smooth functions of the data are preferable for bootstrap.
%    Standard errors are estimated by bootknife resampling by default [2], 
%  where nboot(2) corresponds to the number of bootknife resamples. If 
%  nboot(2) is 0 and standard errors are calculated without resampling 
%  (if bootfun is 'mean'), or using leave-one-out jackknife. Note that if 
%  bootfun is not the mean, the t-statistics returned by this function  
%  will not be comparable with tabulated values.
%
%  bootnhst(...,'nboot',nboot) is a vector of upto two positive integers
%  indicating the number of replicate samples for the first (bootstrap) 
%  and second (bootknife) levels of iterated resampling. The default
%  value of nboot is [1000,200]. If a scalar value is provided for nboot,
%  the value will set the number of first level bootstrap samples; the 
%  number of second level bootknife samples will assume the default of 
%  200, except for the mean, where the default is 0. Increasing the values 
%  of nboot reduces the Monte Carlo error of the p-value (and confidence 
%  interval) estimates but the calculations take longer to complete. If 
%  nboot(2) is explicitly set to 0 then bootnhst calculates standard 
%  errors for studentization using jackknife resampling instead.
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
%  the confidence interval limits returned in columns 8 and 9 are 
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
%   groups   - group number and bootfun for each group with sample size,
%              standard error and confidence intervals that start to overlap
%              at a FWER-controlled p-value of approximately 0.05
%   Var      - weighted mean (pooled) sampling variance
%   maxT     - omnibus test statistic (maxT) 
%   df       - degrees of freedom
%   nboot    - number of bootstrap resamples
%   alpha    - two-tailed significance level for the confidence interval 
%              coverage of 0 (in c).
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
% Bootstrap settings:
%  Function: mean
%  Bootstrap resampling method: Balanced, bootknife resampling
%  Number of bootstrap resamples: 5000
%  Method for estimating standard errors: Calculated without resampling
%  Multiple comparison method: Single-step maxT procedure based on Dunnett
%  Reference group used for comparisons: 1
% ------------------------------------------------------------------------------
% 
% Overall hypothesis test from single-step maxT procedure
% H0: Groups of data are all sampled from the same population as data in ref
% 
% Maximum t(21) = 3.24, p-adj = .018
% ------------------------------------------------------------------------------
% POST HOC TESTS with control of the FWER by the single-step maxT procedure
% ------------------------------------------------------------------------------
% | Comparison |  Reference # |       Test # |  Difference |    t(df)|   p-adj |
% |------------|--------------|--------------|-------------|---------|---------|
% |          1 |            1 |            2 |   -9.48e+00 |    0.93 |    .855 |
% |          2 |            1 |            3 |   -2.49e+01 |    2.43 |    .101 |
% |          3 |            1 |            4 |   -3.32e+01 |    3.24 |    .018 |*
% |          4 |            1 |            5 |   -1.35e+01 |    1.32 |    .600 |
% |          5 |            1 |            6 |   -2.07e+01 |    2.02 |    .216 |
% |          6 |            1 |            7 |   -3.11e+01 |    3.04 |    .031 |*
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
% Bootstrap settings:
%  Function: smoothmedian
%  Bootstrap resampling method: Balanced, bootknife resampling
%  Number of bootstrap resamples: 5000
%  Method for estimating standard errors: Balanced, bootknife resampling
%  Number of bootknife resamples used to estimate standard errors: 200
%  Multiple comparison method: Single-step maxT procedure based on Dunnett
%  Reference group used for comparisons: 1
% ------------------------------------------------------------------------------
% 
% Overall hypothesis test from single-step maxT procedure
% H0: Groups of data are all sampled from the same population as data in ref
% 
% Maximum t(21) = 2.79, p-adj = .022
% ------------------------------------------------------------------------------
% POST HOC TESTS with control of the FWER by the single-step maxT procedure
% ------------------------------------------------------------------------------
% | Comparison |  Reference # |       Test # |  Difference |    t(df)|   p-adj |
% |------------|--------------|--------------|-------------|---------|---------|
% |          1 |            1 |            2 |   -6.66e+00 |    0.52 |    .982 |
% |          2 |            1 |            3 |   -3.02e+01 |    2.36 |    .063 |
% |          3 |            1 |            4 |   -3.56e+01 |    2.79 |    .022 |*
% |          4 |            1 |            5 |   -1.62e+01 |    1.27 |    .552 |
% |          5 |            1 |            6 |   -2.00e+01 |    1.56 |    .343 |
% |          6 |            1 |            7 |   -3.49e+01 |    2.73 |    .026 |*
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
% Bootstrap settings:
%  Function: mean
%  Bootstrap resampling method: Balanced, bootknife resampling
%  Number of bootstrap resamples: 5000
%  Method for estimating standard errors: Calculated without resampling
%  Multiple comparison method: Single-step maxT procedure based on Dunnett
%  Reference group used for comparisons: male
% ------------------------------------------------------------------------------
% 
% Overall hypothesis test from single-step maxT procedure
% H0: Groups of data are all sampled from the same population as data in ref
% 
% Maximum t(9) = 1.20, p-adj = .256
% ------------------------------------------------------------------------------
%
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
% Bootstrap settings:
%  Function: mean
%  Bootstrap resampling method: Balanced, bootknife resampling
%  Number of bootstrap resamples: 5000
%  Method for estimating standard errors: Calculated without resampling
%  Multiple comparison method: Single-step maxT procedure based on Tukey-Kramer
% ------------------------------------------------------------------------------
% 
% Overall hypothesis test from single-step maxT procedure
% H0: Groups of data are all sampled from the same population
% 
% Maximum t(14) = 6.17, p-adj = <.001
% ------------------------------------------------------------------------------
% POST HOC TESTS with control of the FWER by the single-step maxT procedure
% ------------------------------------------------------------------------------
% | Comparison |  Reference # |       Test # |  Difference |    t(df)|   p-adj |
% |------------|--------------|--------------|-------------|---------|---------|
% |          1 |            1 |            2 |   +3.83e+01 |    6.17 |   <.001 |***
% |          2 |            1 |            3 |   +3.50e+00 |    0.59 |    .830 |
% |          3 |            2 |            3 |   -3.48e+01 |    5.60 |   <.001 |***
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
%
%
%   Bibliography:
%   [1] Efron and Tibshirani. Chapter 16 Hypothesis testing with the
%        bootstrap in An introduction to the bootstrap (CRC Press, 1994)
%   [2] Hesterberg, Tim C. (2004), Unbiasing the Bootstrap - Bootknife-
%        Sampling vs. Smoothing, Proceedings of the Section on Statistics 
%        and the Environment, American Statistical Association, 2924-2930.
%
%  bootnhst v2.1.0.0 (24/07/2022)
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

  % Evaluate the number of function arguments
  if (nargin < 2)
    error('bootnhst usage: ''bootnhst (data, group, varargin)''; atleast 2 input arguments required');
  end

  % Store local functions in a stucture for parallel processes
  localfunc = struct ('maxstat',@maxstat,...
                      'empcdf',@empcdf);

  % Check if running in Octave (else assume Matlab)
  info = ver; 
  ISOCTAVE = any (ismember ({info.Name}, 'Octave'));
  
  % Apply defaults
  bootfun = 'mean';
  nboot = [1000,200];
  ref = [];
  alpha = 0.05;
  dim = 1;
  DisplayOpt = true;
  paropt = struct;
  paropt.UseParallel = false;
  if ISOCTAVE
    paropt.nproc = nproc;
  else
    paropt.nproc = feature('numcores');
  end

  % Fetch extra input arguments
  argin3 = varargin;
  narg = numel(argin3);
  if (narg > 1)
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
      else
        error('bootnhst: unrecognised input argument to bootnhst')
      end
      argin3 = {argin3{1:end-2}};
      narg = numel(argin3);
      if (narg < 1)
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
    error('bootnhst: DATA and GROUP must have the same number of rows')
  end
  if iscell(group)
    if ~iscellstr(group)
      group = cell2mat(group);
    end
  end
  if iscell(bootfun)
    func = bootfun{1};
    args = bootfun(2:end);
    bootfun = @(data) func (data, args{:});
  end
  if ~isa(nboot,'numeric')
    error('bootnhst: nboot must be numeric');
  end
  if any(nboot~=abs(fix(nboot)))
    error('bootnhst: nboot must contain positive integers')
  end
  if isa(bootfun,'char')
    if any(strcmpi(bootfun,'robust'))
      bootfun = 'smoothmedian';
    end
    bootfun = str2func(bootfun);
  end
  if numel(nboot) > 2
    error('bootnhst: the vector nboot cannot have length > 2')
  elseif numel(nboot) < 2
    if strcmp(func2str(bootfun),'mean')
      % Avoid resampling when estimating standard errors of the mean
      % nboot(2) must be explcitly set to > 0 to request bootknife standard errors of the mean
      nboot = cat(2, nboot, 0);
    else
      % If nboot is scalar, set default number of bootknife samples to calculate standard errors
      % nboot(2) must be explcitly set to 0 to request jackknife standard errors
      nboot = cat(2, nboot, 200);
    end
  end
  if nboot(1) < 1000
    error('bootnhst: the minimum allowable value of nboot(1) is 1000')
  end 
  if (nboot(2) == 0) && ~strcmp(func2str(bootfun),'mean')
    if ISOCTAVE
      statspackage = ismember ({info.Name}, 'statistics');
      if (~ any (statspackage))
        error ('bootnhst: jackknife calculations require the ''Statistics'' package')
      end
    else
      if ~ismember ('Statistics and Machine Learning Toolbox', {info.Name})
        error ('bootnhst: jackknife calculations require the ''Statistics and Machine Learning Toolbox''')
      end
    end
  end
  
  % Error checking
  if ~isempty(ref) && strcmpi(ref,'pairwise')
    ref = [];
  end
  if ~isa(dim,'numeric')
    error('bootnhst: dim must be numeric');
  end
  if (dim ~= 1) && (dim ~= 2)
    error('bootnhst: dim must be either 1 or 2');
  end
  if nargout > 3
    error('bootnhst only supports up to 3 output arguments')
  end
  if ~islogical(DisplayOpt) || (numel(DisplayOpt)>1)
    error('bootnhst: the value for DisplayOpt must be a logical scalar value')
  end

  % Data or group exclusion using NaN 
  if isnumeric(group)
    if any(isnan(group))
      data(isnan(group),:) = [];
      group(isnan(group)) = [];
    end
  end
  if any(any(isnan([data]),2))
    group(any(isnan([data]),2)) = [];
    data(any(isnan([data]),2),:) = [];
  end

  % Assign non-zero numbers to group labels
  [gnames,junk,g] = unique(group,'legacy');
  clear junk;
  gk = unique(g,'legacy');
  k = numel(gk);
  if ~isempty(ref)
    if isnumeric(ref)
      ref = gk(ismember(gnames,ref));
    else
      ref = gk(strcmp(gnames,ref));
    end
  end
  N = numel(g);
  
  % If applicable, check we have parallel computing capabilities
  if paropt.UseParallel
    if ISOCTAVE  
      pat = '^parallel';
      software = pkg('list');
      names = cellfun(@(S) S.name, software, 'UniformOutput', false);
      status = cellfun(@(S) S.loaded, software, 'UniformOutput', false);
      index = find(~cellfun(@isempty,regexpi(names,pat)));
      if ~isempty(index)
        if logical(status{index})
          PARALLEL = true;
        else
          PARALLEL = false;
        end
      else
        PARALLEL = false;
      end
    else
      if ismember ('Parallel Computing Toolbox', {info.Name})
        PARALLEL = true;
      else
        PARALLEL = false;
      end
    end
  end
  
  % If applicable, setup a parallel pool (required for MATLAB)
  if ~ISOCTAVE
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
          % MATLAB Parallel Computing Toolbox is not installed
          warning('MATLAB Parallel Computing Toolbox is not installed. Falling back to serial processing.')
          paropt.UseParallel = false;
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
  else
    if paropt.UseParallel && (paropt.nproc > 1) && ~PARALLEL
      if ISOCTAVE
        % OCTAVE Parallel Computing Package is not installed or loaded
        warning('OCTAVE Parallel Computing Package is not installed and/or loaded. Falling back to serial processing.')
      else
        % MATLAB Parallel Computing Toolbox is not installed or loaded
        warning('MATLAB Parallel Computing Toolbox is not installed and/or loaded. Falling back to serial processing.')
      end
      paropt.UseParallel = false;
      paropt.nproc = 0;
    end
  end

  % Create handle to a local function for calculating the maximum test statistic
  localfunc = struct ('maxstat', @maxstat);
  func = @(data) localfunc.maxstat (data, g, nboot(2), bootfun, ref, ISOCTAVE);

  % Perform resampling and calculate bootstrap statistics to estimate sampling distribution under the null hypothesis
  boot (1, 1, true, 1, 0); % set random seed to make bootstrap resampling deterministic  
  % Use newer, faster and balanced (less biased) resampling functions (boot and bootknife)
  if paropt.UseParallel
    [null,Q] = bootknife (data,nboot(1),func,[],[],paropt.nproc,[],ISOCTAVE);
  else
    [null,Q] = bootknife (data,nboot(1),func,[],[],0,[],ISOCTAVE);
  end
  
  % Compute the estimate (theta) and it's pooled (weighted mean) sampling variance 
  theta = zeros(k,1);
  SE = zeros(k,1);
  Var = zeros(k,1);
  t = zeros(nboot(2),1);
  nk = zeros(size(gk));
  for j = 1:k
    if (nboot(2) == 0)
      nk(j) = sum(g==gk(j));
      if strcmp (func2str(bootfun), 'mean')
        theta(j) = mean(data(g==gk(j),:));
        % Quick calculation for the standard error of the mean
        SE(j) = std(data(g==gk(j),:),0)/sqrt(nk(j));
        if (j==1); se_method = 'Calculated without resampling'; end;
      else
        theta(j) = bootfun(data(g==gk(j),:));
        % If requested, compute unbiased estimates of the standard error using jackknife resampling
        jackstat = jackknife(bootfun,data(g==gk(j),:));
        SE(j) = sqrt ((nk(j)-1)/nk(j) * sum(((mean(jackstat)-jackstat)).^2));
        if (j==1); se_method = 'Leave-one-out jackknife'; end;
      end
    else
      % Compute unbiased estimate of the standard error by balanced bootknife resampling
      % Bootknife resampling involves less computation than Jackknife when sample sizes get larger
      theta(j) = bootfun(data(g==gk(j),:));
      nk(j) = sum(g==gk(j));
      stats = bootknife(data(g==gk(j),:),[nboot(2),0],bootfun,[],[],0,[],ISOCTAVE);
      SE(j) = stats.std_error;
      if (j==1); se_method = 'Balanced, bootknife resampling'; end;
    end
    Var(j) = ((nk(j)-1)/(N-k)) * SE(j)^2;
  end
  if any(SE==0)
    error('bootnhst: samples must have non-zero standard error')
  end
  if any(isnan(SE))
    error('bootnhst: evaluating bootfun on the bootknife resamples created NaN values for the standard error')
  end
  nk_bar = sum(nk.^2)./sum(nk);  % weighted mean sample size
  Var = sum(Var.*nk/nk_bar);     % pooled sampling variance weighted by sample size
  df = sum(nk)-k;                % degrees of freedom

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
  stats.groups = zeros(k,6);
  stats.groups = zeros(k,6);
  stats.groups(:,1) = gk;
  stats.groups(:,2) = theta;
  stats.groups(:,3) = nk;
  stats.groups(:,4) = SE;
  stats.groups(:,5) = theta - sqrt((0.5*(w+1)).*Var/2) * interp1(cdf,QS,1-alpha,'linear');
  stats.groups(:,6) = theta + sqrt((0.5*(w+1)).*Var/2) * interp1(cdf,QS,1-alpha,'linear');
  stats.Var = Var;
  stats.maxT = maxT;
  stats.df = df;
  stats.nboot = nboot;
  stats.alpha = alpha;
  stats.bootstat = Q;

  % Print output and plot graph with confidence intervals if no output arguments are requested
  cols = [1,2,5,6,7]; % columns in c that we want to print data for
  if (nargout == 0) || (DisplayOpt == true)
    if ~iscellstr(gnames)
      gnames = cellstr(num2str(gnames));
    end
    fprintf (['\n',...
                    'Summary of bootstrap null hypothesis (H0) significance test(s)\n',...
                    '******************************************************************************\n']);
    fprintf ('Bootstrap settings: \n');
    fprintf (' Function: %s\n',func2str(bootfun));
    fprintf (' Bootstrap resampling method: Balanced, bootknife resampling\n')
    fprintf (' Number of bootstrap resamples: %u \n', nboot(1));
    fprintf (' Method for estimating standard errors: %s\n', se_method)
    if (nboot(2) > 0)
      fprintf (' Number of bootknife resamples used to estimate standard errors: %u \n', nboot(2));
    end
    if isempty(ref)
      fprintf (' Multiple comparison method: %s \n', 'Single-step maxT procedure based on Tukey-Kramer');
    else
      fprintf (' Multiple comparison method: %s \n', 'Single-step maxT procedure based on Dunnett');
      fprintf (' Reference group used for comparisons: %s \n', gnames{ref});
    end
    fprintf ('------------------------------------------------------------------------------\n\n');
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
      fprintf (['POST HOC TESTS with control of the FWER by the single-step maxT procedure\n',...
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

%--------------------------------------------------------------------------

function maxT = maxstat (Y, g, nboot, bootfun, ref, ISOCTAVE)

  % Helper function file required for bootnhst
  % Calculate maximum test statistic
  
  % maxstat cannot be a subfunction or nested function since 
  % Octave parallel threads won't be able to find it

  % Calculate the size of the data (N) and the number (k) of unique groups
  N = size(g,1);
  gk = unique(g);
  k = numel(gk);

  % Compute the estimate (theta) and it's pooled (weighted mean) sampling variance 
  theta = zeros(k,1);
  SE = zeros(k,1);
  Var = zeros(k,1);
  nk = zeros(size(gk));
  for j = 1:k
    if (nboot == 0)
      nk(j) = sum(g==gk(j));
      if strcmp (func2str(bootfun), 'mean')
        theta(j) = mean(Y(g==gk(j),:));
        % Quick calculation for the standard error of the mean
        SE(j) = std(Y(g==gk(j),:),0) / sqrt(nk(j));
      else
        theta(j) = bootfun(Y(g==gk(j),:));
        % If requested, compute unbiased estimates of the standard error using jackknife resampling
        jackstat = jackknife(bootfun,Y(g==gk(j),:));
        SE(j) = sqrt ((nk(j)-1)/nk(j) * sum(((mean(jackstat)-jackstat)).^2));
      end
    else
      % Compute unbiased estimate of the standard error by balanced bootknife resampling
      % Bootknife resampling involves less computation than Jackknife when sample sizes get larger
      theta(j) = bootfun(Y(g==gk(j),:));
      nk(j) = sum(g==gk(j));
      stats = bootknife(Y(g==gk(j),:),[nboot,0],bootfun,[],[],0,[],ISOCTAVE);
      SE(j) = stats.std_error;
    end
    Var(j) = ((nk(j)-1)/(N-k)) * SE(j)^2;
  end
  if any(isnan(SE))
    error('maxstat: evaluating bootfun on the bootknife resamples created NaN values for the standard error')
  end
  nk_bar = sum(nk.^2)./sum(nk);  % weighted mean sample size
  Var = sum(Var.*nk/nk_bar);     % pooled sampling variance weighted by sample size

  % Calculate weights to correct for unequal sample size  
  % when calculating standard error of the difference
  w = nk_bar./nk;

  % Calculate the maximum test statistic 
  if isempty(ref)
    % Calculate Tukey-Kramer test statistic (without sqrt(2) factor)
    %
    % Bibliography:
    %  [1] https://en.wikipedia.org/wiki/Tukey%27s_range_test
    %  [2] https://cdn.graphpad.com/faq/1688/file/MulitpleComparisonAlgorithmsPrism8.pdf
    %  [3] www.graphpad.com/guides/prism/latest/statistics/stat_the_methods_of_tukey_and_dunne.htm
    idx = logical(triu(ones(k,k),1));
    i = (1:k)' * ones(1,k);
    j = ones(k,1) * (1:k);
    t = abs(theta(i(idx)) - theta(j(idx))) ./ sqrt(Var * (w(i(idx)) + w(j(idx))));
  else
    % Calculate Dunnett's test statistic 
    t = abs((theta - theta(ref))) ./ sqrt(Var * (w + w(ref)));
  end
  maxT = max(t);
  
end

%--------------------------------------------------------------------------

function [F, x] = empcdf (bootstat, c)

  % Subfunction to calculate empirical cumulative distribution function of bootstat
  %
  % Set c to:
  %  1 to have a complete distribution with F ranging from 0 to 1
  %  0 to avoid duplicate values in x
  %
  % Unlike ecdf, empcdf uses a denominator of N+1

  % Check input argument
  if ~isa(bootstat,'numeric')
    error('bootnhst:empcdf: bootstat must be numeric')
  end
  if all(size(bootstat)>1)
    error('bootnhst:empcdf: bootstat must be a vector')
  end
  if size(bootstat,2)>1
    bootstat = bootstat.';
  end

  % Create empirical CDF
  bootstat = sort(bootstat);
  N = sum(~isnan(bootstat));
  [x,F] = unique(bootstat,'rows','last','legacy');
  F = F/(N+1);

  % Apply option to complete the CDF
  if c > 0
    x = [x(1);x;x(end)];
    F = [0;F;1];
  end

  % Remove impossible values
  F(isnan(x)) = [];
  x(isnan(x)) = [];
  F(isinf(x)) = [];
  x(isinf(x)) = [];

end

%--------------------------------------------------------------------------

%!test
%! y = [111.39 110.21  89.21  76.64  95.35  90.97  62.78;
%!      112.93  60.36  92.29  59.54  98.93  97.03  79.65;
%!       85.24 109.63  64.93  75.69  95.28  57.41  75.83;
%!      111.96 103.40  75.49  76.69  77.95  93.32  78.70];
%! g = [1 2 3 4 5 6 7;
%!      1 2 3 4 5 6 7;
%!      1 2 3 4 5 6 7;
%!      1 2 3 4 5 6 7];
%! p = bootnhst (y(:),g(:),'ref',1,'nboot',[1000,0],'DisplayOpt',false);
%! assert (p, 0.01210939735473963, 1e-09);
%! p = bootnhst (y(:),g(:),'nboot',[1000,0],'DisplayOpt',false);
%! assert (p, 0.04407742932277153, 1e-09);
%! # Result from anova1 is 0.0387

%!test
%! y = [54       43
%!      23       34
%!      45       65
%!      54       77
%!      45       46
%!     NaN       65];
%! g = {'male' 'female'
%!      'male' 'female'
%!      'male' 'female'
%!      'male' 'female'
%!      'male' 'female'
%!      'male' 'female'};
%! p = bootnhst (y(:),g(:),'ref','male','nboot',[1000,0],'DisplayOpt',false);
%! assert (p, 0.2577543618442567, 1e-09);
%! p = bootnhst (y(:),g(:),'nboot',[1000,0],'DisplayOpt',false);
%! assert (p, 0.2577543618442567, 1e-09);
%! # Result from anova1 is 0.2613

%!test
%! y = [54  87  45
%!      23  98  39
%!      45  64  51
%!      54  77  49
%!      45  89  50
%!      47 NaN  55];
%! g = [ 1   2   3
%!       1   2   3
%!       1   2   3
%!       1   2   3
%!       1   2   3
%!       1   2   3];
%! p = bootnhst (y(:),g(:),'nboot',[1000,0],'DisplayOpt',false);
%! assert (p, 0.001, 1e-09); # truncated at 0.001
%! # Result from anova1 is 4.162704768129188e-05

