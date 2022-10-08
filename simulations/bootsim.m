% BOOTSTRAP CONFIDENCE INTERVAL SIMULATION SCRIPT
% EVALUATION OF INTERVAL COVERAGE, TAIL PROBABILITIES, LENGTH AND SHAPE
clear

% Save warning states and turn off all warnings
state = warning('query','all');
warning('off','all');

% Set significance level
alpha = 0.05;

% Function of the data
func = @mean;

% Population parameter
population_param = 0;

% Define sample size
n = 7;

% Define number of simulations
sim = 1000;

% Reset outcome counters
accept = 0;
reject = 0;
above = 0;
below = 0;

% Bootstrap resampling
nboot = 2000;
type = 'cal';

% Print settings
fprintf('----- BOOTSTRAP CONFIDENCE INTERVAL SIMULATION -----\n')
fprintf('Simulation size: %u\n',sim);
fprintf('Sample size: %u\n',n);
fprintf('Algorithm: bootci\n');
fprintf('nboot: %u\n',nboot);
fprintf('Type: %s\n',type);
fprintf('Statistic: %s\n',char(func));
fprintf('Alpha: %.3f\n',alpha);

% Initialize simulation variables
length = nan(sim,1);
shape  = nan(sim,1);
coverage  = nan(sim,1);

% Parallel processing
paropt = struct ('UseParallel', true, 'nproc', 4);
  
% Start simulation
for i=1:sim

  % Create random sample
  x = randn(n,1);

  % Bootstrap confidence interval
  ci = bootci(nboot,{func,x},'alpha',alpha,'type',type,'Options',paropt);
  stat = func(x);

  % Coverage counter
  if population_param >= ci(1) && population_param <= ci(2)
   accept = accept + 1;
  else
   reject = reject + 1;
  end
  if population_param < ci(1)
    below = below + 1;
  elseif population_param > ci(2)
    above = above + 1;
  end
  coverage(i) = accept/(accept+reject);
  if i>1
    fprintf(repmat('\b', 1, 13))
  end
  fprintf('%s: %.1f%s',...
          sprintf('%05s',num2str(i)),...
          100*coverage(i),'%');

  % Interval lengths
  length(i) = ci(2) - ci(1);

  % Interval shape
  shape(i) = (ci(2) - stat) / (stat - ci(1));

end

% Print results
fprintf(repmat('\b', 1, 13))
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