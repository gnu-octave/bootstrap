%  Function File: bootknife
%
%  Bootknife resampling  
%
%  The function uses a variant of bootstrap called bootknife sampling, 
%  which involves creating leave-one-out jackknife samples of size n - 1 
%  and then drawing samples of size n with replacement from the jackknife 
%  samples [1]. The resampling of data values is balanced to reduce Monte 
%  Carlo error [2]. By default, the algorithm uses double bootstrap to
%  provide an accurate bias-corrected parameter estimate, standard error,
%  and 95% confidence intervals.
%
%  stats = bootknife(nboot,bootfun,d1,...,dN)
%  ... = bootknife(...,'alpha',alpha)
%  ... = bootknife(...,'Strata',strata)
%  [stats,bootstat] = bootknife(...)
%  [stats,bootstat] = bootknife(...)
%  [stats,bootstat,bootsam] = bootknife(...)
%
%  stats = bootknife(nboot,bootfun,d) draws nboot bootknife resamples from
%  rows of data and computes the statistic defined by bootfun to
%  approximate the sampling distribution, correct the statistic for bias, 
%  and estimate standard errors and confidence intervals.
%
%  nboot is a scalar, or vector of upto two positive integers indicating
%  the number of replicate samples for the first and second bootstraps.
%  By default, nboot is [2000,200], which implements double bootstrap. If 
%  nboot provided is scalar, the default number of resamples in the second 
%  bootstrap is 200. We can get away with a lower number of resamples in 
%  the second bootstrap (to reduce the computational expense of the double 
%  bootstrap) since linear interpolation is used to achieve near asymptotic 
%  calibration of confidence intervals [3].
%
%  bootfun is a function handle (e.g. specified with @), or a string 
%  indicating the function name. 
%  
%  The third input argument is data (column vector or a matrix), that is 
%  used to create inputs for bootfun. 
%  
%  The first output argument returns a row vector containing:
%    1) The bias-corrected statistic
%    2) The bootstrap standard error
%    3) Lower limit of the 95% bootstrap confidence interval
%    4) Upper limit of the 95% bootstrap confidence interval
%
%  stats = bootknife(nboot,bootfun,d1,...,dN) is as above except that 
%  the third and subsequent numeric input arguments are data vectors 
%  that are used to create inputs for bootfun. 
%
%  ... = bootknife(nboot,{bootfun,...},...,'strata',strata) specifies a
%  vector containing numeric identifiers of strata. The dimensions of
%  strata must be equal to that of the non-scalar input arguments to
%  bootfun. Bootstrap resampling is stratified so that every stratum is
%  represented in each bootstrap test statistic.
%
%  ... = bootknife(nboot,bootfun,...,'alpha',alpha) computes the
%  iterated bootstrap confidence interval of the statistic defined by the
%  function bootfun with lower and uppper interval ends calibrated to 
%  100 * alpha/2 and 100 * (1-alpha)/2 %, where alpha is a scalar
%  value between 0 and 1. The default value of alpha is 0.05 corresponding 
%  to intervals with a coverage of 95% confidence.
%
%  [stats,bootstat] = bootknife(...) also returns bootstat, a vector 
%  bootstrap statistics for the (first level of) bootknife resampling. 
%
%  [stats,bootstat,bootsam] = bootknife(...) also returns bootsam, a  
%  matrix of indices for bootknife sampling. Each column in bootsam
%  corresponds to one bootknife sample and contains the row 
%  indices of the values drawn from the nonscalar data argument 
%  to create that sample.
%
%  Bibliography:
%  [1] Hesterberg T.C. (2004) Unbiasing the Bootstrap—Bootknife Sampling 
%        vs. Smoothing; Proceedings of the Section on Statistics & the 
%        Environment. Alexandria, VA: American Statistical Association.
%  [2] Davison et al. (1986) Efficient Bootstrap Simulation.
%        Biometrika, 73: 555-66
%  [3] Hall, Lee and Young (2000) Importance of interpolation when
%        constructing double-bootstrap confidence intervals. Journal
%        of the Royal Statistical Society. Series B. 62(3): 479-491
%
%  bootknife v1.3.0.0 (16/05/2022)
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


function [stats, T1, idx] = bootknife (nboot, bootfun, varargin)
  
  % Error checking
  if nargin < 3
    error('not enough input argumemnts')
  end

  % Set defaults
  alpha = 0.05;
  strata = [];
  
  % Assign input arguments to function variables
  argin3 = varargin;
  narg = numel(argin3);
  if narg > 1
    while ischar(argin3{end-1})
      if strcmpi(argin3{end-1},'alpha')
        alpha = argin3{end};   
      elseif strcmpi(argin3{end-1},'strata')
        strata = argin3{end};   
      else
        error('unrecognised input argument to bootstrp')
      end
      argin3 = {argin3{1:end-2}};
      narg = numel(argin3);
      if narg < 3
        break
      end
    end
  end
  x = argin3;

  % Determine properies of data
  nvar = numel(x);
  sz = size(x{1});
  n = sz(1);

  % Initialize
  B = nboot(1);
  if (numel(nboot) > 1)
    C =  nboot(2);
  else
    C = 200;
  end
  T1 = zeros(1, B);
  U = zeros(1, B);
  M = zeros(1, B);
  V = zeros(1, B);
  X = cell(1, nvar);
  idx = zeros(n,B);
  c = ones(n,1) * B;

  % Perform balanced bootknife resampling
  % Octave or Matlab serial/vectorized computing
  %    Gleason, J.R. (1988) Algorithms for Balanced Bootstrap Simulations. 
  %    The American Statistician. Vol. 42, No. 4 pp. 263-266
  if ~isempty(strata)
    % Get strata IDs
    gid = unique(strata);  % strata ID
    K = numel(gid);        % number of strata
    % Create strata matrix
    g = zeros(n,K);
    for k = 1:K
      g(:,k) = (strata == gid(k));
      [~, ~, idx(g(:,k),:)] = bootknife (B, bootfun, x{:}(g(:,k),:));
    end
    % Perform data sampling
    for v = 1:nvar
      X{v} = x{v}(idx);
    end
    for b = 1:B
      T1(b) = feval(bootfun,X{:}(:,b));
    end
  else
    for b = 1:B
      % Choose which rows of the data to sample
      r = b - fix((b-1)/n) * n;
      for i = 1:n
        d = c;   
        d(r) = 0;
        if ~sum(d)
          d = c;
        end
        j = sum((rand(1) >= cumsum(d./sum(d)))) + 1;
        idx(i,b) = j;
        c(j) = c(j) - 1;
      end 
      % Perform data sampling
      for v = 1:nvar
        X{v} = x{v}(idx(:,b));
      end
      % Function evaluation on bootknife sample
      T1(b) = feval(bootfun,X{:}); 
    end
  end
  
  % Calculate the bootstrap standard error, bias and confidence intervals 
  % Bootstrap standard error estimation
  T0 = bootfun(x{:});
  if C > 0
    % Iterated bootstrap resampling for greater accuracy
    for b = 1:B
      [~,T2] = bootknife ([C,0], bootfun, x{:}(idx(:,b)));
      % Use quick interpolation to find the probability that T2 <= T0
      I = (T2<=T0);
      u = sum(I);
      U(b) =  interp1q([max([min(T2), max(T2(I))]);...
                        min([max(T2), min(T2(~I))])],...
                       [u; min(u+1,C)] / C,...
                       T0);
      if isnan(U(i))
        U(b) = u / C;
      end
      M(b) = mean(T2);
      V(b) = var(T2,1);
    end
    % Double bootstrap bias estimation
    % See Ouysee (2011) Economics Bulletin
    bias = mean(T1) - T0 - mean(M - T1);
    % Double bootstrap standard error
    se = sqrt(var(T1,1)^2 / mean(V));
    % Calibrate tail probabilities to half of alpha
    l = quantile (U, [alpha/2, 1-alpha/2]);
    % Calibrated percentile bootstrap confidence intervals
    ci = quantile (T1, l);
  else
    % Bootstrap bias estimation
    bias = mean(T1) - T0;
    % Bootstrap standard error
    se = std(T1,1);
    % Percentile bootstrap confidence intervals
    ci = quantile (T1, [alpha/2, 1-alpha/2]);
  end
  
  % Prepare output
  stats = [T0-bias; se; ci.'];
  
end
