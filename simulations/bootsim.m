% BOOTSTRAP CONFIDENCE INTERVAL SIMULATION SCRIPT
% EVALUATION OF INTERVAL COVERAGE, TAIL PROBABILITIES, LENGTH AND SHAPE
clear

% Save warning states and turn off all warnings
state = warning('query','all');
warning('off','all');

% Set significance level
alpha = .1;         % 95% confidence interval
%alpha = [.025,.975]; % 95% confidence interval

%--------------------------------------------------------------
% Uncomment one of the following example simulation conditions

% 1) Mean of a normal distribution (no skewness or excess kurtosis)
%bootfun = @mean; nvar = 1; rnd = @(n) random ('norm', 0, 1, [n, 1]); theta = 0;

% 2) Mean of a folded normal distribution (high skewness, low excess kurtosis)
%bootfun = @mean; nvar = 1; rnd = @(n) abs (random ('norm', 0, 1, [n, 1])); theta = sqrt (2/pi);

% 3) Mean of an exponential distribution (no skewness, high excess kurtosis)
%bootfun = @mean; nvar = 1; rnd = @(n) random ('exp', 1, [n, 1]); theta = 1;   

% 4) Mean of a log-normal distribution (high skewness, high kurtosis)
%bootfun = @mean; nvar = 1; rnd = @(n) random ('logn', 0, 1, [n, 1]); theta = exp (0.5); 

% 5) Variance of a normal distribution
bootfun = @(x) var(x,1); nvar = 1; rnd = @(n) random ('norm', 0, 1, [n, 1]); theta = 1;

% 6) Correlation coefficient of bivariate normal distribution (rho = 0)
%bootfun =  @cor; nvar = 2; rnd = @(n) random ('norm', 0, 1, [n, 1]); theta = 0;

% 7) Linear regression (through the origin) of bivariate normal distribution (rho = 0)
%bootfun =  @regress; nvar = 2; rnd = @(n) random ('norm', 0, 1, [n, 1]); theta = 0;
%--------------------------------------------------------------

% Define sample size
n = 26;

% Define number of simulations
sim = 1000;

% Reset outcome counters
accept = 0;
reject = 0;
above = 0;
below = 0;

% Bootstrap resampling
nboot = 2000;
type = 'stud';

% Print settings
fprintf('----- BOOTSTRAP CONFIDENCE INTERVAL SIMULATION -----\n')
fprintf('Simulation size: %u\n',sim);
fprintf('Sample size: %u\n',n);
fprintf('Algorithm: bootci\n');
fprintf('nboot: %u\n',nboot);
fprintf('Type: %s\n',type);
fprintf('Statistic: %s\n',char(bootfun));
fprintf('Alpha: %.3f\n',alpha);

% Initialize simulation variables
length = nan (sim,1);
shape  = nan (sim,1);
coverage  = nan (sim,1);

% Parallel processing
ncpus = 4;
paropt = struct ('UseParallel', true, 'nproc', ncpus);
  
% Start simulation
for i=1:sim

  % Create random sample
  x = rnd (n);
  if (nvar > 1)
    y = rnd (n);
  end

  % Calculate statistic from the sample
  if (nvar > 1)
    stat = feval (bootfun, x, y);
  else
    stat = feval (bootfun, x);
  end
  
  % Bootstrap confidence interval
  ci = bootci (nboot, {bootfun,x}, 'alpha', alpha, 'type', type, 'Options', paropt);
  %if (nvar > 1)
  %  S = bootknife ({x, y}, nboot, bootfun, alpha, [], ncpus); ci = [S.CI_lower; S.CI_upper];
  %else
  %  S = bootknife (x, nboot, bootfun, alpha, [], ncpus); ci = [S.CI_lower; S.CI_upper];
  %end
  

  % Coverage counter
  if theta >= ci(1) && theta <= ci(2)
   accept = accept + 1;
  else
   reject = reject + 1;
  end
  if theta < ci(1)
    below = below + 1;
  elseif theta > ci(2)
    above = above + 1;
  end
  coverage(i) = accept/(accept+reject);
  if i>1
    fprintf(repmat('\b', 1, 14))
  end
  fprintf('%s: % 5s%s',...
          sprintf('%06d',i),...
          sprintf('%.1f',round(1000*coverage(i))/10),'%');

  % Interval lengths
  length(i) = ci(2) - ci(1);

  % Interval shape
  shape(i) = (ci(2) - stat) / (stat - ci(1));

end

% Print results
fprintf(repmat('\b', 1, 14))
fprintf(['%s: %.1f%s\n',...
         '%s: %.2f (%.2f-%.2f)\n',...
         '%s: %.2f (%.2f-%.2f)\n',...
         '%s: %.4f\n',...
         '%s: %.4f\n',...
         '%s: %.4f\n'],...
          'Coverage',100 * coverage(end),'%',...
          'Length',median(length),quantile(length,0.25),quantile(length,0.75),...
          'Shape',median(shape),quantile(shape,0.25),quantile(shape,0.75),...
          'Lower tail rejection',below/sim,...
          'Upper tail rejection',above/sim);

% Restore initial warning states
warning(state);
