function bootout = iboot (x, T0, nboot, bootfun, strata, isoctave)

    % Helper function file for bootknife required for accelerating double bootstrap by parallel processing

    % Initialize output
    bootout = cell(1);
    bootout{1} = zeros(2,1);
    
    % Perform inner level of resampling
    [junk, T2] = bootknife (x, nboot, bootfun, [], strata, 0, isoctave);
    
    % Use quick interpolation to find the probability that T2 <= T0
    I = (T2 <= T0);
    u = sum (I);
    t2 = [max([min(T2), max(T2(I))]),...
          min([max(T2), min(T2(~I))])];
    if (u < nboot) && ((t2(2) - t2(1)) > 0)
      % Linear interpolation to calculate U, which is required to calibrate alpha and improving confidence interval coverage 
      % U is the probability that T2 (for this iteration of the inner layer of resamplng) <= T0
      bootout{1}(1) = ((t2(2) - T0) * u / nboot + (T0 - t2(1)) * min ((u + 1) / nboot, 1)) /...
                       (t2(2) - t2(1));
    else
      bootout{1}(1) = u / nboot;
    end
    
    % Calculate M required for double bootstrap estimate of bias 
    % M is the mean of the bootstrap distribution (for this iteration of the inner layer of resamplng) 
    bootout{1}(2) = mean (T2);    


end