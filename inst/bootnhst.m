%  Function File: bootnhst
%
%  Bootstrap null hypothesis significance tests (NHST)
%
%  p = bootnhst(DATA,GROUP)
%  p = bootnhst(DATA,GROUP)
%  p = bootnhst(...,'bootfun',bootfun)
%  p = bootnhst(...,'nboot',nboot)
%  p = bootnhst(...,'ref',ref)
%  p = bootnhst(...,'cluster',clusters)
%  p = bootnhst(...,'Options',paropt)
%  p = bootnhst(...,'alpha',alpha)
%  [p,c] = bootnhst(DATA,GROUP,...)
%  [p,c,stats] = bootnhst(DATA,GROUP,...)
%  bootnhst(DATA,GROUP,...);
%
%  This non-parametric (or semi-parametric) bootstrap function can be used 
%  for null hypothesis (H0) significance testing with univariate (vector) 
%  or multivatiate (matrix) data, to compare bootfun (default is the 'mean') 
%  evaluated on independent GROUPs (i.e. samples) [1]. This function is
%  appropriate for posthoc comparisons among a family of hypothesis tests.
%
%  The tests conducted by bootnhst do not make the normality assumption 
%  of parametric statistical tests and the calculations of the weighted mean 
%  sampling variance (for studentization using a pooled standard error) 
%  accomodates for unequal sample size and gives bootnhst more statistical 
%  power (compared to not pooling the sampling variance). Since DATA across 
%  the GROUPs is resampled, as for a permutatiion test, bootnhst assumes 
%  exchangeability among the groups under the null hypothesis. Note that 
%  this function will return an error if any GROUP is not represented by 
%  more than one data row.
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
%  The method used for bootstrap is balanced resampling (unless computations 
%  are accelerated by parallel processing). 
%
%  p = bootnhst(DATA,GROUP) is a k-sample bootstrap test where GROUP is a 
%  vector or cell array the same number of rows as y and contains numbers 
%  or strings which denote GROUP labels. If the number of GROUPs (k) is 2, 
%  this is a 2-sample test. If k > 2, this is test for k GROUPs (like in 
%  one-way ANOVA layout). If GROUP is numeric and any GROUP is assigned 
%  NaN then their respective DATA rows will be excluded from analysis.
%
%  p = bootnhst(...,'bootfun',bootfun) sets the statistic calculated
%  from the bootstrap samples. This can be a function handle or string
%  corresponding to the function name. The calculation of bootfun on the
%  data must return a scalar value. If empty, the default is @mean or 
%  'mean'. If DATA is multivariate, bootfun is the grand mean, which is 
%  is the mean of the means of each column (i.e. variates). The standard 
%  errors are estimated from 200 bootknife resamples [2], or cluster-
%  jacknife if clusters are specified. If a robust statistic for central 
%  location is required, setting bootfun to 'robust' implements a smoothed 
%  version of the median (see function help for smoothmedian).
%
%  p = bootnhst(...,'nboot',nboot) sets the number of bootstrap resamples.
%  Increasing nboot reduces the monte carlo error of the p-value estimates 
%  but the calculations take longer to complete. When nboot is empty or not 
%  provided, the default (and minimum allowable nboot to compute two-tailed 
%  p-values down to 0.001) is 1000 - an error is returned if the nboot 
%  provided by the user is lower than this. 
%
%  p = bootnhst(...,'ref',ref) also sets the GROUP to use as the reference 
%  GROUP for post hoc tests. For a one-way experimental design or family of 
%  tests, this will usually be a control GROUP. If all pairwise comparisons 
%  are desired, then set ref to 'pairwise' or leave empty. By default, 
%  pairwise comparisons are computed for post hoc tests.
%
%  p = bootnhst(...,'cluster',clusters) specifies a column vector of numeric 
%  identifiers with the same number of rows as DATA. The identifiers should 
%  indicate cluster membership of the data rows. Clusters are resampled by 
%  two-stage bootstrap resampling of residuals with shrinkage correction 
%  (see bootstrp help for more information). Because specifying clusters 
%  causes bootnhst to resample residuals, the bootstrap becomes semi-
%  parametric.The clusters input argument can be used to accomodate for a 
%  single level of nesting in heirarchical data structures. This resampling 
%  stategy is appropriate when the family of tests has a split plot design 
%  layout. If left empty, this argument is ignored. This function will 
%  return an error if any GROUP is not represented by more than one cluster, 
%  but there is no restriction on the number of data rows in any cluster.
%
%  p = bootnhst(...,'Options',paropt) specifies options that govern if and 
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
%   'nproc'       - The number of processors to use by Octave. Default
%                   is the number of available processors. If you choose
%                   In Matlab, nproc is ignored and the number of parallel
%                   workers should be predefined beforehand by starting
%                   a parallel pool, else it will use the preferred number
%                   of workers.
% 
%  p = bootnhst(...,'alpha',alpha) specifies the two-tailed significance 
%  level for confidence interval coverage of 0 (in c) or interval overlap 
%  (in stats.groups).
%
%  [p, c] = bootnhst(DATA,GROUP,...) also returns a 9 column matrix that
%  summarises post hoc test results. The family-wise error rate is 
%  simultaneously controlled since the null distribution for each test 
%  represents the studentized range statistic. 
%
%  If the ref input argument is empty, the resampling procedure above 
%  for pairwise comparisons is analagous to Tukey's Honest Significance 
%  Difference (HSD) test. 
%
%  If the ref input argument is specified, then the resampling above 
%  for treatment versus reference is analagous to Dunnett's post hoc  
%  tests (since the range of bootfun in the null distribution test  
%  statistics is restricted to differences with the reference GROUP).
% 
%  The q-ratio (analagous to a t-statistic) is computed using the 
%  difference in bootfun evaluated for two GROUPs divided by the  
%  standard error of the difference (derived from the pooled, weighted 
%  mean, sampling variance). To compare the q-ratio reported here with 
%  Tukey's more traditional q-statistic, multiply it by sqrt(2). Note 
%  that because unbiased sampling variance is estimated using bootknife 
%  resampling [2], bootnhst can be used to compare a wide variety of 
%  statistics (not just the mean). 
% 
%  The columns of output argument c contain:
%    column 1: reference GROUP number
%    column 2: test GROUP number
%    column 3: value of bootfun evaluated using reference GROUP DATA
%    column 4: value of bootfun evaluated using test GROUP DATA
%    column 5: columns 4 minus column 3
%    column 6: q-ratio
%    column 7: p-value
%    column 8: LOWER bound of the bootstrap-t confidence interval
%    column 9: UPPER bound of the bootstrap-t confidence interval
%
%  [p,c] = bootnhst(data,group,...) also returns the group names used in 
%  the GROUP input argument. The index of gnames corresponds to the 
%  numbers used to identify GROUPs in columns 1 and 2 of the output 
%  argument c.
%
%  [p,c,stats] = bootnhst(DATA,GROUP,...) also returns a structure 
%  containing additional statistics. The stats structure contains the 
%  following fields:
%
%   gnames   - group names used in the GROUP input argument. The index of 
%              gnames corresponds to the numbers used to identify GROUPs
%              in columns 1 and 2 of the output argument c
%   ref      - reference group
%   groups   - bootfun for each group with sample size, standard error, 
%              and lower and upper bootstrap-t confidence intervals, which  
%              have coverage such that they overlap with eachother if the  
%              ref input argument is 'pairwise', or with the reference   
%              group, at a FWER-controlled p-value of greater than alpha.
%   Var      - weighted mean (pooled) sampling variance
%   stat     - omnibus test statistic (q) 
%   nboot    - number of bootstrap resamples
%   alpha    - two-tailed significance level for the confidence interval 
%              coverage of 0 (in c) or interval overlap (in stats.groups)
%   clusters - vector of numeric identifiers indicating cluster membership
%   bootstat - test statistic computed for each bootstrap resample 
%
%  bootnhst(DATA,GROUP,...); performs the calculations as per above but 
%  prints the columns 1, 2 and 5-7 of the results (c) in a pretty table.
%  The differences between groups are also plot along with the symmetic
%  95% bootstrap-t confidence intervals (adjusted to control the FWER).
%  Markers and error bars are red if p < 0.05 or blue if p > 0.05.
%
%   Bibliography:
%   [1] Efron and Tibshirani. Chapter 16 Hypothesis testing with the
%        bootstrap in An introduction to the bootstrap (CRC Press, 1994)
%   [2] Hesterberg, Tim C. (2004), Unbiasing the Bootstrap - Bootknife-
%        Sampling vs. Smoothing, Proceedings of the Section on Statistics 
%        and the Environment, American Statistical Association, 2924-2930.
%
%  bootnhst v1.5.0.0 (26/01/2022)
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

  % Apply defaults
  bootfun = 'mean';
  nboot = [1000,200];
  ref = [];
  alpha = 0.05;
  strata = [];
  clusters = [];
  blocksize = [];
  paropt = struct;
  paropt.UseParallel = false;
  % Initialise nproc if it doesn't exist
  if ~exist('nproc','builtin')
    nproc = 1;
  end
  paropt.nproc = nproc;

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
      elseif any(strcmpi(argin3{end-1},{'cluster','clusters'}))
        clusters = argin3{end};
      else
        error('unrecognised input argument to bootstrp')
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
  if isa(bootfun,'function_handle')
    if all(bootfun(data) == mean(data))
      if nvar > 1
        % Grand mean for multivariate data
        bootfun = @(data) mean(mean(data));
      else
        bootfun = @mean;
      end
    elseif all(bootfun(data) == smoothmedian(data))
      if nvar > 1 
        % Grand smoothed median for multivariate data
        bootfun = @(data) smoothmedian(smoothmedian(data));
      else
        bootfun = @smoothmedian;
      end
    elseif all(bootfun(data) == median(data))
      error('bootfun cannot be the median, use ''robust'' instead.')
    end
  elseif isa(bootfun,'char')
    if strcmpi(bootfun,'mean') 
      if nvar > 1
        % Grand mean for multivariate data
        bootfun = @(data) mean(mean(data));
      else
        bootfun = @mean;
      end
    elseif any(strcmpi(bootfun,{'robust','smoothmedian'}))
      if nvar > 1
        % Grand smoothed median for multivariate data
        bootfun = @(data) smoothmedian(smoothmedian(data));
      else
        bootfun = @smoothmedian;
      end
    elseif strcmpi(bootfun,'median')
      error('bootfun cannot be the median, use ''robust'' instead.')
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
    % set default number of bootknife samples;
    nboot = cat(2,nboot,200);
  end
  if nboot(1) < 1000
    error('the minimum allowable value of nboot is 1000')
  end 
  if ~isempty(ref) && strcmpi(ref,'pairwise')
    ref = [];
  end
  if ~isempty(clusters)
    if (size(clusters,2) > 1) || (size(clusters,1) ~= size(data,1))
      error('clusters must be a column vactor with the same number of rows as DATA')
    end
    opt = struct;
  end
  if nargout > 3
    error('bootnhst only supports up to 3 output arguments')
  end

  % Group exclusion using NaN 
  if isnumeric(group)
    if any(isnan(group))
      data(isnan(group),:) = [];
      group(isnan(group)) = [];
    end
  end

  % Assign non-zero numbers to group labels
  [gnames,~,g] = unique(group);
  if isempty(clusters)
    N = size(g,1);
  else
    N = numel(unique(clusters));
  end
  gk = unique(g);
  k = numel(gk);
  if ~isempty(ref)
    if isnumeric(ref)
      ref = gk(ismember(gnames,ref));
    else
      ref = gk(strcmp(gnames,ref));
    end
  end
  
  % Define function to calculate maximum difference among groups
  func = @(data) maxstat(data,g,nboot(2),bootfun,ref,clusters);

  % Perform resampling and calculate bootstrap statistics
  state = warning;
  warning off;    % silence warnings about non-vectorized bootfun
  Q = bootstrp (nboot(1),func,data,'cluster',clusters,'Options',paropt);
  warning(state);

  % Compute the estimate (theta) and it's pooled (weighted mean) sampling variance 
  theta = zeros(k,1);
  SE = zeros(k,1);
  Var = zeros(k,1);
  t = zeros(nboot(2),1);
  nk = zeros(size(gk));
  for j = 1:k
    theta(j) = feval(bootfun,data(g==gk(j),:));
    if ~isempty(clusters)
      % Compute unbiased estimate of the standard error by cluster-jackknife resampling
      opt = struct;
      opt.clusters = clusters(g==gk(j));
      nk(j) = numel(unique(opt.clusters));
      SE(j) = jack (data(g==gk(j),:), bootfun, [], opt);
    else
      % Compute unbiased estimate of the standard error by bootknife resampling
      % Bootknife resampling involves less computation than Jackknife when sample sizes get larger
      nk(j) = sum(g==gk(j));
      if nvar > 1
        t = zeros(nboot(2),1); 
        for b = 1:nboot(2)
          idx = 1+fix(rand(nk(j)-1,1)*nk(j));
          tmp = data(g==gk(j),:);
          t(b) = feval(bootfun,tmp(idx,:));
        end
      else
        % Vectorized if data is univariate
        idx = 1+fix(rand(nk(j)-1,nboot(2))*nk(j));
        tmp = data(g==gk(j),:);
        t = feval(bootfun,tmp(idx));
      end
      SE(j) = std(t);
    end
    Var(j) = ((nk(j)-1)/(N-k)) * SE(j)^2;
  end
  if any(nk <= 1)
    error('the number of observations or clusters per group must be greater than 1')
  end
  nk_bar = sum(nk.^2)./sum(nk);  % weighted mean sample size
  Var = sum(Var.*nk/nk_bar);     % pooled sampling variance weighted by sample size

  % Calculate weights to correct for unequal sample size  
  % when calculating standard error of the difference
  w = nk_bar./nk;

  % Prepare to make symmetrical bootstrap-t confidence intervals
  [cdf,QS] = empcdf(Q,0);

  % Calculate p-values for comparisons adjusted to simultaneously control the FWER
  if isempty(ref)
    % Resampling version of Tukey's test
    % Note that Tukey's q-ratio here does not have the sqrt(2) factor. 
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
        c(i,7) = interp1(QS,1-cdf,c(i,6),'linear',0);
      end
      c(i,8) = c(i,5) - SED * interp1(cdf,QS,1-alpha,'linear');
      c(i,9) = c(i,5) + SED * interp1(cdf,QS,1-alpha,'linear');
    end
  else
    % Resampling version of Dunnett's test
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
        c(j,7) = interp1(QS,1-cdf,c(j,6),'linear',0);
      end
      c(j,8) = c(j,5) - SED * interp1(cdf,QS,1-alpha,'linear');
      c(j,9) = c(j,5) + SED * interp1(cdf,QS,1-alpha,'linear');
    end
    c(ref,:) = [];
  end

  % Calculate overall p-value
  q = max(c(:,6));
  p = sum(Q>=q)/nboot(1);

  % Prepare stats output structure
  % Include bootstrap-t confidence intervals for bootfun evaluated for each group
  % Within a small margin of error, these confidence intervals overlap at a FWER controlled p-value > 0.05
  stats = struct;
  stats.gnames = gnames;
  stats.ref = ref;
  stats.groups = zeros(k,5);
  stats.groups(:,1) = theta;
  stats.groups(:,2) = nk;
  stats.groups(:,3) = SE;
  stats.groups(:,4) = theta - sqrt((0.5*(w+1)).*Var/2) * interp1(cdf,QS,1-alpha,'linear');
  stats.groups(:,5) = theta + sqrt((0.5*(w+1)).*Var/2) * interp1(cdf,QS,1-alpha,'linear');
  stats.Var = Var;
  stats.stat = q;
  stats.nboot = nboot;
  stats.alpha = alpha;
  stats.clusters = clusters;
  stats.bootstat = Q;

  % Print output and plot graph with confidence intervals if no output arguments are requested
  cols = [1,2,5,6,7]; % columns in c that we want to print data for
  if (nargout == 0)
    fprintf (['\n',...
                    'Summary of bootstrap null hypothesis (H0) significance test(s)\n',...
                    '******************************************************************************\n']);
    if isempty(ref)
      fprintf ('H0: Groups of data are all sampled from the same population\n\n');
    else
      fprintf ('H0: Groups of data are all sampled from the same population as data in ref\n\n');
    end
    if (p < 0.001)
      fprintf (['Overall test p-value: <0.001 \n',...
                '------------------------------------------------------------------------------\n']);
    else
      fprintf (['Overall test p-value: %.3f \n',...
                '------------------------------------------------------------------------------\n'],p);
    end
    if size(c,1) > 1
      fprintf (['\n',...
                      'POST HOC TESTS (with simultaneous control of the FWER)\n',...
                      '------------------------------------------------------------------------------\n',...
                      '| Comparison |  Reference # |       Test # |  Difference | q-ratio |  p-value|\n',...
                      '|------------|--------------|--------------|-------------|---------|---------|\n']);
      if isempty(ref)
        for i = 1:n
          tmp = num2cell(c(i,cols));
          if (c(i,7) < 0.001)
            tmp(end) = [];
            fprintf ('| %10u | %12u | %12u | %+11.2e | %7.3f |  <0.001 |***\n',i,tmp{:});
          else
            fprintf ('| %10u | %12u | %12u | %+11.2e | %7.3f |   %.3f |',i,tmp{:});
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
          if (c(j,7) < 0.001)
            tmp(end) = [];
            fprintf ('| %10u | %12u | %12u | %+11.2e | %7.3f |  <0.001 |***\n',j,tmp{:});
          else
            fprintf ('| %10u | %12u | %12u | %+11.2e | %7.3f |   %.3f |',j,tmp{:});
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
      fprintf (['\n',...
                '------------------------------------------------------------------------------\n',...
                '|    GROUP # |                                                         GROUP |\n',...
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

