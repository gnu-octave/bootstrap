function [m1, m2, S] = BCa (nboot, func, data, bootstat, stat, alpha, S, paropt, opt)

  % Private function file required for ibootci

  % Note that alpha input argument is nominal coverage

  % Create distribution functions
  stdnormcdf = @(x) 0.5 * (1 + erf (x / sqrt(2)));
  stdnorminv = @(p) sqrt (2) * erfinv (2 * p-1);

  % Get number of primary sampling units
  n = S.n(end);

  % Calculate bias correction z0
  z0 = stdnorminv (sum (bootstat < stat)/ nboot);
  if isinf(z0) || isnan(z0)
    error('Unable to calculate bias correction z0, please use percentile intervals instead')
  end

  % Use the Jackknife to calculate acceleration constant
  if nargin > 8
    [SE, T, U] = jack (data, func, paropt, opt);
    if any(diff(opt.weights)) && isempty(opt.clusters)
      w = opt.weights/sum(opt.weights);
    else
      w = 1/n .* ones(n,1);
    end
  else
    [SE, T, U] = jack (data, func, paropt, opt);
    w = 1/n .* ones(n,1);
  end
  a = ((1 / 6) * ((sum(w .* U.^3)) / (sum(w .* U.^2))^(3/2))) / sqrt(n);

  % Calculate BCa percentiles
  z1 = stdnorminv(0.5 * (1 + alpha));
  m1 = stdnormcdf(z0 + ((z0 + z1)/(1 - a * (z0 + z1))));
  z2 = stdnorminv(0.5 * (1 - alpha));
  m2 = stdnormcdf(z0 + ((z0 + z2)/(1 - a * (z0 + z2))));
  S.z0 = z0;
  S.a = a;

end
