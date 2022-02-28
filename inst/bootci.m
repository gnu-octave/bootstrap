%  Function File: bootci
%
%  Two-sided nonparametric bootstrap confidence intervals
%
%  ci = bootci(nboot,bootfun,d)
%  ci = bootci(nboot,{bootfun,d1,...,dN},...,'alpha',alpha)
%  ci = bootci(nboot,{bootfun,...},...,'type',type)
%  ci = bootci(nboot,{bootfun,...},...,'type','student','nbootstd',nbootstd)
%  ci = bootci(nboot,{bootfun,...},...,'type','student','stderr',stderr)
%  ci = bootci(nboot,{bootfun,...},...,'Weights',weights)
%  ci = bootci(nboot,{bootfun,...},...,'Options',paropt)
%  [ci,bootstat] = bootci(...)
%
%  ci = bootci(nboot,bootfun,...) computes the 95% bootstrap confidence
%  interval of the statistic computed by bootfun. nboot is a scalar value
%  indicating the number of replicate bootstrap samples. bootfun is a
%  function handle (e.g. specified with @), or a string indicating the
%  function name. The third and later input arguments are data (column
%  vectors, or a matrix), that are used to create inputs for bootfun.
%  bootci creates each first level bootstrap by balanced resampling from
%  the rows of the data argument(s) (which must be the same size) [1-3].
%  Linear interpolation of the empirical cumulative distribution function
%  of bootstat is then used to construct two-sided confidence intervals
%  [4]. The default value for nboot is 5000.
%
%  ci = bootci(nboot,{bootfun,...},...,'alpha',alpha) computes the
%  bootstrap confidence interval of the statistic defined by the function
%  bootfun with coverage 100*(1-alpha)%, where alpha is a scalar value
%  between 0 and 1. bootfun and the data that bootci passes to it are
%  contained in a single cell array. The default value of alpha is 0.05,
%  corresponding to intervals with a coverage of 95% confidence.
%
%  ci = bootci(nboot,{bootfun,...},...,'type',type) computes the bootstrap
%  confidence interval of the statistic defined by the function bootfun.
%  type is the confidence interval type, chosen from among the following:
%    'norm' or 'normal': Normal approximated interval with
%                        bootstrapped bias and standard error.
%    'per' or 'percentile': Percentile method.
%    'cper': Bias-corrected percentile method.
%    'bca': Bias-corrected and accelerated method (Default).
%    'stud' or 'student': Studentized (bootstrap-t) confidence interval.
%    The bootstrap-t method includes an additive correction to stabilize
%    the variance when the sample size is small [5].
%
%  ci = bootci(nboot,{bootfun,...},...,'type','stud','nbootstd',nbootstd)
%  computes the Studentized bootstrap confidence interval of the statistic
%  defined by the function bootfun. The standard error of the bootstrap
%  statistics is estimated using bootstrap, with nbootstd bootstrap data
%  samples. nbootstd is a positive integer value. The default value of
%  nbootstd is 50. The nbootstd argument is ignored when the interval
%  type is set to anything other than Studentized (bootstrap-t) intervals.
%
%  ci = bootci(nboot,{bootfun,...},...,'type','stud','stderr',stderr)
%  computes the studentized bootstrap confidence interval of statistics
%  defined by the function bootfun. The standard error of the bootstrap
%  statistics is evaluated by the function stderr. stderr is a function
%  handle. stderr takes the same arguments as bootfun and returns the
%  standard error of the statistic computed by bootfun.
%
%  ci = bootci(nboot,{bootfun,...},...,'Weights',weights) specifies
%  observation weights. weights must be a vector of non-negative numbers.
%  The length of weights must be equal to first dimension of the non-
%  scalar input argument(s) to bootfun. Balanced resampling is extended
%  to resampling with weights [6], which are used as bootstrap sampling
%  probabilities. Note that weights are not currently implemented for 
%  Studentized-type intervals.
%
%  ci = bootci(nboot,{bootfun,...},...,'Options',paropt) specifies
%  options that govern if and how to perform bootstrap and jacknife
%  iterations using multiple processors (if the Parallel Computing
%  Toolbox or Octave Forge parallel package is available). This 
%  argument is a structure with the following recognised fields:
%
%   'UseParallel' - If true, compute bootstrap iterations in parallel.
%                   Default is false for serial computation. In MATLAB,
%                   the default is true if a parallel pool has already
%                   been started.
%   'nproc'       - The number of processors to use by Octave. Default
%                   is the number of available processors. If you choose
%                   In Matlab, nproc is ignored and the number of parallel
%                   workers should be predefined beforehand by starting
%                   a parallel pool, else it will use the preferred number
%                   of workers.
%
%  [ci,bootstat] = bootci(...) also returns the bootstrapped statistic
%  computed for each of the bootstrap replicate samples sets.
%
%  Bibliography:
%  [1] Efron, and Tibshirani (1993) An Introduction to the
%        Bootstrap. New York, NY: Chapman & Hall
%  [2] DiCiccio and Efron (1996) Bootstrap confidence intervals.
%        Statist. Sci. 11(3):189-228.
%  [3] Davison et al. (1986) Efficient Bootstrap Simulation.
%        Biometrika, 73: 555-66
%  [4] Efron (1981) Censored data and the bootstrap. JASA
%        76(374): 312-319
%  [5] Polansky (2000) Stabilizing bootstrap-t confidence intervals
%        for small samples. Can J Stat. 28(3):501-516
%  [6] Booth, Hall and Wood (1993) Balanced Importance Resampling
%        for the Bootstrap. The Annals of Statistics. 21(1):286-298
%
%  bootci v2.8.7.1 (23/12/2021)
%  Author: Andrew Charles Penn
%  https://www.researchgate.net/profile/Andrew_Penn/
%
%  Cite as:
%  Andrew Penn (2019). bootci (https://www.github.com/acp29/iboot), GitHub.
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


function [ci,bootstat] = bootci(argin1,argin2,varargin)

  % Evaluate the number of function arguments
  if nargin<2
    error('Too few input arguments');
  end
  if nargout>5
   error('Too many output arguments');
  end

  % Initialise nproc if it doesn't exist
  if ~exist('nproc','builtin')
    nproc = 1;
  end

  % Assign input arguments to function variables
  if ~iscell(argin2)
    % Normal usage without options
    nboot = argin1;
    bootfun = argin2;
    data = varargin;
    alpha = 0.05;
    type = 'bca';
    weights = [];
    nbootstd = [];
    stderr = [];
    paropt = struct;
    paropt.UseParallel = false;
    paropt.nproc = nproc;
  else
    % Normal usage with options
    % Evaluate option input arguments
    nboot = argin1;
    bootfun = argin2{1};
    data = {argin2{2:end}};
    options = varargin;
    alpha = 1+find(strcmpi('alpha',options));
    type = 1+find(strcmpi('type',options));
    nbootstd = 1+find(strcmpi('nbootstd',options));
    stderr = 1+find(strcmpi('stderr',options));
    weights = 1+find(strcmpi('weights',options));
    paropt = 1+find(cellfun(@(options) any(strcmpi({'options','option'},options)),options));
    if ~isempty(alpha)
      try
        alpha = options{alpha};
      catch
        alpha = 0.05;
      end
    else
      alpha = 0.05;
    end
    if ~isempty(type)
      try
        type = options{type};
      catch
        type = 'bca';
      end
    else
      type = 'bca';
    end
    if any(strcmpi(type,{'stud','student'}))
      if ~isempty(nbootstd)
        try
          nbootstd = options{nbootstd};
        catch
          nbootstd = 100;
        end
      else
        nbootstd = 100;
      end
    else
      nbootstd = [];
    end
    if ~isempty(stderr)
      try
        stderr = options{stderr};
        nbootstd = [];
      catch
        stderr = [];
      end
    else
      stderr = [];
    end
    if ~isempty(weights)
      try
        weights = options{weights};
      catch
        weights = [];
      end
    else
      weights = [];
    end
    if ~isempty(paropt)
      try
        paropt = options{paropt};
        if isstruct(paropt)
          if isfield(paropt,'UseParallel')
            if ~islogical(paropt.UseParallel)
              error('value of parallel option field ''UseParallel'' must be logical')
            end
          else
            paropt.UseParallel = false;
          end
          if isfield(paropt,'nproc') && isoctave
            if ~isa(paropt.nproc,'numeric') || numel(paropt.nproc)~=1
              error('value of parallel option field ''nproc'' must be scalar')
            end
            if paropt.nproc ~= fix(paropt.nproc)
              error('value of parallel option field ''nproc'' must be an integer')
            end
            if paropt.nproc < 1
              error('value of parallel option field ''nproc'' must greater than 0')
            end
            if paropt.nproc > nproc
              error('nproc cannot be > the number of physical cores')
            end
          else
            paropt.nproc = nproc;
          end
        else
          error('parallel options must be contained in a structure')
        end
      catch
        error(lasterr);
      end
    else
      paropt = struct;
      paropt.UseParallel = false;
      paropt.nproc = nproc;
    end
  end
  % Evaluate data input
  opt = struct;
  nvar = size(data,2);
  if (min(size(data{1}))>1)
    if (nvar == 1)
      opt.matflag = 1;   % Flag for matrix input set to 1
    end
  else
    opt.matflag = 0;
  end
  opt.blocksize = [];
  opt.clusters = [];
  opt.weights = [];
  
  % Set 'type' input variable for ibootci
  if any(strcmpi(type,{'norm','normal'}))
    % Calculate percentile intervals (we will modify the intervals later)
    inptype = 'per';
  else
    inptype = type;
  end

  % Error checking
  if ~all(size(nboot) == [1,1])
    error('nboot must be a scalar value')
  end

  % Create options structure
  options = {'alpha', alpha,...
             'type', inptype,...
             'nbootstd', nbootstd,...
             'stderr', stderr,...
             'Weights', weights,...
             'Options', paropt};

  % Parse input arguments to the function ibootci
  [ci, bootstat, S] = ibootci (nboot, {bootfun, data{:}}, options{:});

  % Apply normal approximation
  switch lower(type)
    case {'norm','normal'}
      % Normal approximation intervals
      % Davison and Hinkley (1997), p198-200
      stdnorminv = @(p) sqrt(2) * erfinv (2 * p - 1);
      za = stdnorminv(alpha/2);   % normal confidence point
      LL = S.bc_stat + S.SE * za;
      UL = S.bc_stat - S.SE * za;
      ci = [LL; UL];     
  end

end
