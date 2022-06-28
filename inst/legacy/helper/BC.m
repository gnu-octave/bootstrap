function [m1, m2, S] = BC (B, T1, T0, alpha, S)

  % Private function file required for ibootci

  % Note that alpha input argument is nominal coverage

  % Create distribution functions
  stdnormcdf = @(x) 0.5*(1+erf(x/sqrt(2)));
  stdnorminv = @(p) sqrt(2)*erfinv(2*p-1);

  % Calculate bias correction z0
  z0 = stdnorminv(sum(T1<T0)/B);
  if isinf(z0) || isnan(z0)
    error('Unable to calculate bias correction z0, please use percentile intervals instead')
  end

  % Calculate Bias-corrected percentiles
  z1 = stdnorminv(0.5*(1+alpha));
  m1 = stdnormcdf(2*z0+z1);
  z2 = stdnorminv(0.5*(1-alpha));
  m2 = stdnormcdf(2*z0+z2);
  S.z0 = z0;
  S.a = 0;

end
