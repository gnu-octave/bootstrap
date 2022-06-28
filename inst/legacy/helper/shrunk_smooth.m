function x_shrunk = shrunk_smooth (x,bandwidth,xbar,xvar,noise)

  % Helper function file required for ibootci

  % Bootstrap smoothing with shrinkage correction to maintain
  % the sample variance
  %
  % Fryer (1976) Some errors associated with the non-parametric
  %   estimation of density functions. J. Inst. Maths Applics.
  %   18, 371-380
  %
  % Jones (1991) On correcting for variance inflation in kernel
  %   density estimation. Comput Stat Data An. 11, 3-15
  %
  % Wang (1995) OPTIMIZING THE SMOOTHED BOOTSTRAP
  %   Ann. Inst. Statist. Math. Vol. 47, No. 1, 65-80
  %
  % Tymoteusz Wolodzko 2019-11-07
  % https://www.rdocumentation.org/packages/kernelboot/versions/0.1.6/topics/kernelboot

  % Correction by shrinking x residuals
  %x_shrunk = (x - xbar) .* sqrt(1 - bandwidth^2 / xvar) +...
  %           xbar + noise;

  % Correction by shrinking both x residuals and the bandwidth
  x_shrunk = sqrt(xvar) * (x - xbar + noise) ./...
             sqrt(bandwidth^2 + xvar) + xbar;

end
