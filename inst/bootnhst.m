%  Function File: bootnhst
%
%  Bootstrap null hypothesis significance test (NHST)
%
%  p = bootnhst(data,group)
%  p = bootnhst(data,group,ref)
%  p = bootnhst(data,group,ref,nboot)
%  p = bootnhst(data,group,ref,nboot,bootfun)
%  p = bootnhst(data,group,ref,nboot,bootfun,paropt)
%  [p, c] = bootnhst(data,group,...)
%  [p, c, gnames] = bootnhst(data,group,...)
%  bootnhst(data,group,...);
%
%  This bootstrap function is for null hypothesis (H0) significance testing 
%  and can be used with univariate (vector) or multivatiate (matrix) data [1].
%  The assumptions of these bootstrap hypothesis tests are similar to the 
%  those of permutation tests.
% 
%  The specification of H0 (and the post-hoc tests if the number of groups 
%  is more than 2) depend on whether a reference group is specified with the 
%  ref input argument.
%
%  If ref is empty: 
%    H0: Groups of data are all sampled from the same population. 
%  If ref is specified:
%    H0: Groups of data are all sampled from the same population as data in ref.
%
%  The method used for bootstrap is balanced resampling (unless computations 
%  are accelerated by parallel processing). 
%
%  Note that matrix input for data will mean that calculation of bootfun 
%  will not be vectorized. When this is the case, the function will 
%  return a harmless warning saying so.
%
%  p = bootnhst(data,group) is a k-sample bootstrap permutation test where 
%  group is a vector or cell array the same size as y containing numbers or 
%  strings which denote group labels. If the number of groups (k) is 2, this 
%  is a 2-sample test. If k > 2, this is test for k groups (like in one-way 
%  ANOVA layout)
%   
%  p = bootnhst(data,group,nboot) sets the number of bootstrap resamples. 
%  The default is 5000.
%
%  p = bootnhst(data,group,nboot,bootfun) also sets the statistic calculated
%  from the bootstrap samples. This can be a function handle or string
%  corresponding to the function name. The default is @mean or 'mean'.
%
%  p = bootnhst(y,...,nboot,bootfun,paropt) specifies options that govern 
%  if and how to perform bootstrap iterations using multiple processors 
%  (if the Parallel Computing Toolbox or Octave Forge package is available).
%  This argument is a structure with the following recognised fields:
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
%  [p, c] = bootnhst(data,group,...) also returns a 6 column matrix that
%  summarises post-hoc test results. The family-wise error rate is 
%  simultaneously controlled since the null distribution for each test 
%  represents the maximum test statistic from each bootstrap sample. The 
%  test statistic is the difference between groups with respect to the 
%  values calculated using bootfun.
%    column 1: reference group number
%    column 2: test group number
%    column 3: value of bootfun evaluated using data from reference group
%    column 4: value of bootfun evaluated using data from test group
%    column 5: columns 4 minus column 3
%    column 6: p-value adjusted for multiple comparisons
%
%  [p, c, gnames] = bootnhst(data,group,...) also returns the group names
%  used in the group input argument. The index of gnames corresponds to the 
%  numbers used to identify groups in columns 1 and 2 of the output argument
%  c.
%
%  bootnhst(data,group,...); performs the calculations as per above but 
%  prints the results in a pretty table.
%
%   Bibliography:
%   [1] Efron and Tibshirani. Chapter 16 Hypothesis testing with the
%        bootstrap in An introduction to the bootstrap (CRC Press, 1994)
%
%  bootnhst v1.0.0.0 (01/01/2022)
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


function [p, c, gnames] = bootnhst (data, group, ref, nboot, bootfun, paropt)

  % Check and process bootnhst input arguments
  nvar = size(data,2);
  if (nargin < 2)
    error('bootnhst requires atleast two input arguments');
  end
  if ischar(group)
    group = cellstr(group);
  end
  if (size(group,1)>1) && (size(data,1) ~= size(group,1))
    error('data and group must be vectors the same size')
  end
  if iscell(group)
    if ~iscellstr(group)
      group = cell2mat(group);
    end
  end
  if (nargin < 3)
    ref = [];
  end
  if (nargin < 4) || isempty(nboot)
    nboot = 5000;
  end
  if any(size(nboot)>1)
    error('nboot must be scalar. bootnhst is not compatible with bootstrap iteration')
  end
  if (nargin < 5) || isempty(bootfun)
    bootfun = 'mean';
  end
  if (nargin < 6) || isempty(paropt)
    paropt = struct;
    paropt.UseParallel = false;
    % Initialise nproc if it doesn't exist
    if ~exist('nproc','builtin')
      nproc = 1;
    end
    paropt.nproc = nproc;
  end

  % Assign non-zero numbers to group labels
  [gnames,~,g] = unique(group);
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
  % H0: Data is exchangeable across all the groups labelled in group
  func = @(data) maxdiff(data,g,ref,bootfun,nvar);
  t = func(data);

  % Perform resampling and calculate bootstrap statistics
  T = bootstrp(nboot,func,data,'Options',paropt);

  % Calculate overall p-value
  p = sum(T>=t)/nboot;
  if (p == 0)
    warning('p-value too small to calculate. Try increasing nboot.')
  end

  % Calculate p-values for comparisons adjusted to control the FWER
  if isempty(ref)
    A = ones(k,1)*gk';
    B = tril(gk*ones(1,k),-1);
    M = [A(:) B(:)];
    ridx = (M(:,2)==0);
    M(ridx,:)=[];
    n = size(M,1);
    c = nan(n,6);
    c(:,1:2) = M;
    for i = 1:n
      c(i,3) = feval(bootfun,data(g==M(i,1),:));
      c(i,4) = feval(bootfun,data(g==M(i,2),:));
      c(i,5) = c(i,4)-c(i,3);
      c(i,6) = sum(T>=abs(c(i,5)))/nboot;
    end
  else
    c = nan(k,6);
    c(:,1) = ref;
    c(:,3) = feval(bootfun,data(g==ref,:));
    for j = 1:k
      c(j,2) = gk(j);
      c(j,4) = feval(bootfun,data(g==c(j,2),:));
      c(j,5) = c(j,4)-c(j,3);
      c(j,6) = sum(T>=abs(c(j,5)))/nboot;
    end
    c(ref,:) = [];
  end

  % Print output
  if (nargout == 0)
    disp (...
          sprintf (['\n',...
                    'Summary of bootstrap null hypothesis (H0) significance test(s)\n',...
                    '*******************************************************************************']));
    if isempty(ref)
      disp (...
            sprintf ('H0: Groups of data are all sampled from the same population.\n'));
    else
      disp (...
            sprintf ('H0: Groups of data are all sampled from the same population as data in ref.\n'));
    end
    if (p < 0.001)
      disp (sprintf (['Overall test p-value: <0.001 \n',...
                      '-------------------------------------------------------------------------------']));
    else
      disp (sprintf (['Overall test p-value: %.3f \n',...
                      '-------------------------------------------------------------------------------'],p));
    end
    if (k > 2)
      disp (...
            sprintf (['\n',...
                      'Post-hoc tests with FWER controlled by a single-step maxT procedure.\n',...
                      '-------------------------------------------------------------------------------\n',...
                      '|Reference # |     Test # |   Reference |       Test | Difference |     p-value\n',...
                      '|------------|------------|-------------|------------|------------|------------']));
      if isempty(ref)
        for i = 1:n
          tmp = num2cell(c(i,:));
          if (c(i,6) < 0.001)
            tmp(end) = [];
            disp (...
                  sprintf ('| %10u | %10u |    %+.2e|   %+.2e|   %+.2e|      <0.001',tmp{:}));
          else
            disp (...
                  sprintf ('| %10u | %10u |    %+.2e|   %+.2e|   %+.2e|       %.3f',tmp{:}));
          end
        end
      else
        for j = 1:k-1
          tmp = num2cell(c(j,:));
          if (c(j,6) < 0.001)
            tmp (end) = [];
            disp (...
                  sprintf ('| %10u | %10u |    %+.2e|   %+.2e|   %+.2e|      <0.001',tmp{:}));
          else
            disp (...
                  sprintf ('| %10u | %10u |    %+.2e|   %+.2e|   %+.2e|       %.3f',tmp{:}));
          end
        end
      end
      disp (sprintf (' ')) % this prints an empty line at the end
    end
  end

end
