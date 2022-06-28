%     Function file:
%
%     [h0, h1] = kdehist(x, kernel, h0, h1)
%
%     Plotting overlaid histograms and kernel density estimates.
%
%     Bibliography:
%     Scott (1979) Biometrika. 66(3):605-610
%     Scott (1992) Multivariate Density Estimation. John Wiley & Sons, Inc.
%     Terrell and Scott (1985) JASA. 80(389):209-214
%
%     kdehist v0.2 (last updated: 21/08/2018)
%     Author: Andrew Charles Penn
%     https://www.researchgate.net/profile/Andrew_Penn/

function [h0, h1] = kdehist(x, kernel, h0, h1)

  % Helper function file for plotboot

  if nargin < 2
    kernel = 'epanechnikov';
  end

  x=x(:);

  % Calculate total number of data points
  n=numel(x);

  % Calculate number of unique data points
  % Use this for calculation of bandwidth and bin widths
  nu=numel(unique(x));

  % Scott's rules
  % Also see Scott (1979) Biometrika. 66(3):605-610
  % Calculate standard deviation of the data in reference set (x1)
  % In cases of skewness or kurtosis, the use of IQR is more robust and limits oversmoothing
  IQR=prctile(x,75)-prctile(x,25);
  sigma=min(std(x),IQR/1.349);

  % Settings for KDE and histogram
  switch lower(kernel)
    case {'epan','epanechnikov'}
      % Bin width bounds; nearly optimal for a variety of densities
      if nargin < 3
        h0=3.55*sigma*nu^(-1/3);  % For histogram
      end
      if nargin < 4
        h1=2.44*sigma*nu^(-1/5);  % For kernel density estimate using Epanechnikov kernel
      end
      % Epanechnikov (parabolic) kernel
      K = @(u) max(0,(3/4)*(1-u.^2));
    case {'norm','normal','gaussian'}
      % Optimal bin width and bandwidth for Gaussian data
      if nargin < 3
        h0=3.49*sigma*nu^(-1/3);  % For histogram
      end
      if nargin < 4
        h1=1.059*sigma*nu^(-1/5); % For kernel density estimate using Gaussian kernel
      end
      % Gaussian kernel
      K = @(u) 1/sqrt(2*pi)*exp(-0.5*u.^2);
  end

  % Calculate vector of bin centers
  lo=min(x)-h1;
  hi=max(x)+h1;
  nbins=round((hi-lo)/h0);
  binEdges=linspace(lo,hi,nbins);
  centers=diff(binEdges)/2+binEdges(1:nbins-1);

  % Calculate and plot histograms
  [y] = hist(x,centers);
  y=y/(h0*n);
  stairs([binEdges(1) binEdges(1:end)],[0 y 0],'color',[0.2 0.2 0.6],'lineWidth',1.0);

  hold on
  % Modify graph properties
  set(gca,'Box','off');
  hold off

  % Calculate and plot Kernel Density estimates
  N=1000;
  u=linspace(lo,hi,N)';
  Z=zeros(N,n);
  for i=1:n
    Z(:,i)=K((u-x(i))/h1);
  end
  z=sum(Z,2);
  z=z/(h1*n); % Kh(u) = 1/hn K(u/h)
  hold on; plot(u,z,'-','color',[0.2 0.2 0.6],'linewidth',3.0); hold off;
