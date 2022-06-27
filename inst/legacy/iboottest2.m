%  Function File: iboottest2
%
%  Two-sample bootstrap test
%
%   p = iboottest2(nboot,{bootfun,x,y})
%   p = iboottest2(nboot,{bootfun,x,y},Name,Value)
%   [p,ci] = iboottest2(...)
%   [p,ci,bootstat] = iboottest2(...)
%   [p,ci,bootstat,S] = iboottest2(...)
%
%  Two sample bootstrap test for unpaired data. The null hypothesis
%  is that the difference between the bootfun statistic calculated
%  for independent samples x and y is equal to zero. The test is
%  two-tailed.
%
%  Note that this function resamples the rows of data x and y.
%
%  See ibootci documentation for input argument definitions and
%  for Name-Value pairs. The following input arguments differ in
%  their format:
%
%  'Weights': The 'Weights' option must be provided as a cell array:
%  the first cell should contain weights for x and the second cell
%  should contain weights for y. An empty cell signifies that no
%  weights will be used in the bootstrap for that sample.
%
%  'Strata': The same format as for weights.
%
%  'Cluster': The same format as for weights.
%
%  'Block': The same format as for weights.
%
%  'Smooth': The same format as for weights.
%
%  iboottest2 v1.5.9.0 (01/10/2021)
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


function [p,ci,bootstat,S] = iboottest2(argin1,argin2,varargin)

  % Check and process iboottest2 input arguments
  nboot = argin1;
  bootfun = argin2{1};
  x = argin2{2};
  y = argin2{3};
  if numel(argin2)>3
    error('too many data variables provided')
  end
  iter = numel(nboot);
  if iter > 2
    error('Size of nboot exceeds maximum number of iterations supported by ibootci')
  end
  if iter==0
    B = 5000;
    C = 200;
    nboot = [B C];
  elseif iter==1
    B = nboot;
    C = 0;
    nboot = [B C];
  elseif iter==2
    B = nboot(1);
    C = nboot(2);
  end

  % Check number of output arguments requested
  if nargout > 4
    error('Too many output arguments requested')
  end
  
  % Retrieve some ibootci options
  options = varargin;
  type = 1+find(strcmpi('type',options));
  alpha = 1+find(strcmpi('alpha',options));
  weights = 1+find(strcmpi('Weights',options));
  strata = 1+find(sum([strcmpi('Strata',options);strcmpi('Stratum',options);strcmpi('Stratified',options)]));
  clusters = 1+find(sum([strcmpi('Cluster',options);strcmpi('Clusters',options)]));
  blocksize = 1+find(sum([strcmpi('Block',options);strcmpi('Blocks',options);strcmpi('Blocksize',options)]));
  bandwidth = 1+find(sum([strcmpi('smooth',options);strcmpi('smoothing',options);strcmpi('bandwidth',options)]));
  cellref = [];
  if ~isempty(type)
    try
      cellref = cat(2,cellref,[type-1,type]);
      type = options{type};
    catch
      type = 'per';
      cellref(end-1:end)=[];
    end
    if strcmpi(type,'bca')
      error('BCa correction is incompatible with iboottest2')
    end
  else
    type = 'per';
  end
  if ~isempty(alpha)
    try
      cellref = cat(2,cellref,[alpha-1,alpha]);
      alpha = options{alpha};
    catch
      alpha = 0.05;
      cellref(end-1:end)=[];
    end
  else
    alpha = 0.05;
  end
  if ~isa(alpha,'numeric') || numel(alpha)~=1
    error('The alpha value must be a numeric scalar value');
  end
  if (alpha <= 0) || (alpha >= 1)
    error('The alpha value must be a value between 0 and 1');
  end
  if ~isempty(weights)
    try
      cellref = cat(2,cellref,[weights-1,weights]);
      weights = options{weights};
    catch
      weights = {[],[]};
      cellref(end-1:end)=[];
    end
    if ~iscell(weights)
      error('The option argument for weights should be a cell array')
    else
      if numel(weights) ~= 2
        error('The option argument for weights should be an array of 2 cells')
      end
    end
  else
    weights = {[],[]};
  end
  if ~isempty(strata)
    try
      cellref = cat(2,cellref,[strata-1,strata]);
      strata = options{strata};
    catch
      strata = {[],[]};
      cellref(end-1:end)=[];
    end
    if ~iscell(strata)
      error('The option argument for strata should be a cell array')
    else
      if numel(strata) ~= 2
        error('The option argument for strata should be an array of 2 cells')
      end
    end
  else
    strata = {[],[]};
  end
  if ~isempty(clusters)
    try
      cellref = cat(2,cellref,[clusters-1,clusters]);
      clusters = options{clusters};
    catch
      clusters = {[],[]};
      cellref(end-1:end)=[];
    end
    if ~iscell(clusters)
      error('The option argument for clusters should be a cell array')
    else
      if numel(clusters) ~= 2
        error('The option argument for clusters should be an array of 2 cells')
      end
    end
  else
    clusters = {[],[]};
  end
  if ~isempty(blocksize)
    try
      cellref = cat(2,cellref,[blocksize-1,blocksize]);
      blocksize = options{blocksize};
    catch
      blocksize = {[],[]};
      cellref(end-1:end)=[];
    end
    if ~iscell(blocksize)
      error('The option argument for blocksize should be a cell array')
    else
      if numel(blocksize) ~= 2
        error('The option argument for blocksize should be an array of 2 cells')
      end
    end
  else
    blocksize = {[],[]};
  end
  if ~isempty(bandwidth)
    try
      cellref = cat(2,cellref,[bandwidth-1,bandwidth]);
      bandwidth = options{bandwidth};
    catch
      bandwidth = {[],[]};
      cellref(end-1:end)=[];
    end
    if ~iscell(bandwidth)
      error('The option argument for smooth should be a cell array')
    else
      if numel(bandwidth) ~= 2
        error('The option argument for smooth should be an array of 2 cells')
      end
    end
  else
    bandwidth = {[],[]};
  end
  options(cellref)=[];   % remove these evaluated options from the options array

  % Perform independent resampling from x and y
  state = warning;
  warning off;
  [ciX,bootstatX,SX] = ibootci(nboot,{bootfun,x},'type',type,'Weights',weights{1},'Strata',strata{1},'Cluster',clusters{1},'Block',blocksize{1},'Smooth',bandwidth{1},options{:});
  [ciY,bootstatY,SY] = ibootci(nboot,{bootfun,y},'type',type,'Weights',weights{2},'Strata',strata{2},'Cluster',clusters{2},'Block',blocksize{2},'Smooth',bandwidth{2},options{:});

  if C>0
    if ~isempty(weights{1})
      [ciX,bootstatX{1}] = ibootci(B,{bootfun,x},'type',type,'alpha',SX.cal,'Weights',weights{1},'Strata',strata{1},'Cluster',clusters{1},'Smooth',bandwidth{1},options{:});
    end
    if ~isempty(weights{2})
      [ciY,bootstatY{1}] = ibootci(B,{bootfun,y},'type',type,'alpha',SY.cal,'Weights',weights{2},'Strata',strata{2},'Cluster',clusters{2},'Smooth',bandwidth{2},options{:});
    end
  end
  warning(state);

  % Calculate differences between bootfun evaluated for bootstrap
  % sample sets derived from x and y
  if iscell(bootstatX)
    bootstat = cell(2,1);
    bootstat{1} = bootstatX{1} - bootstatY{1};
    if size(bootstatX{2},1) > 1
      bootstat{2} = bootstatX{2} - bootstatY{2};
    else
      bootstat{2} = sqrt(bootstatX{2}.^2 + bootstatY{2}.^2);  % propagation of error for stderr
    end
  else
    bootstat = bootstatX - bootstatY;
  end

  % Create template settings structure and calculate the sample test statistic
  S = SX;
  T0 = SX.stat - SY.stat;
  S.stat = T0;                       % assign correct sample test statistic in S
  S.alpha = alpha;                   % reset alpha in S
  S.n = SX.n + SY.n;                 % calculate total sample size
  S.SE = sqrt(SX.SE.^2 + SY.SE.^2);  % propagation of error (required if stderr option used)

  % Calculate confidence interval using ibootci
  [ci,bootstat,S,calcurve] = ibootci(bootstat, S);

  % Update output structure
  S.ICC = [SX.ICC, SY.ICC];
  S.DEFF = [SX.DEFF, SY.DEFF];
  S.xcorr = {SX.xcorr, SY.xcorr};
  S.strata = {SX.strata,SY.strata};
  S.clusters = {SX.clusters,SY.clusters};
  S.blocksize = [SX.blocksize,SY.blocksize];
  S.bandwidth = [SX.bandwidth,SY.bandwidth];
  S.weights = {SX.weights,SY.weights};

  % Calculate p-value using ibootp
  p = ibootp(0,bootstat,S,calcurve);

end
