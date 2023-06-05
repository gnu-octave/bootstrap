% BAYESIAN CREDIBLE INTERVAL SIMULATION SCRIPT
clear

% bootbayes settings
nboot = 10000;
prob = 0.95;

% Simulation settings
n = 7;
mu = 0;
sigma = 1;
nsamp = 1e+04;                  % Total sample size desired from simulation
nsim = 1e+07;                   % Size of simulation block
eps = 5e-02 * sigma * sqrt (n); % Stringency/precision
mu_good = [];

% Create random sample
y = randn (n, 1) * sigma + mu; 
mu_within_range_of_y = (min(y) <= mu) & (mu <= max(y));
range = max(y) - min(y);
s = std (y, 1); % Standard deviation (sample as a population)

% Compute credible interval
stats = bootbayes (y, [], [], nboot, prob);
CI = [stats.CI_lower stats.CI_upper];
mu_within_range_of_CI = (CI(1) <= mu) & (mu <= CI(2));  % Frequentist inference

% Print settings
LogicalStr = {'false', 'true'};
fprintf ('----- BAYESIAN BOOTSTRAP CREDIBLE INTERVAL (CI) SIMULATION -----\n')
fprintf ('Sample size: %u\n', n);
fprintf ('nboot: %u\n', nboot);
fprintf ('Prior: %.2f\n', 1.0);
if (numel(prob) > 1)
  fprintf (' Credible interval (CI) type: Percentile interval\n');;
  fprintf (' Credible interval: %.3g%% (%.1f%%, %.1f%%)\n', mass, 100 * prob);
else
  fprintf ('Credible interval (CI) type: Shortest probability interval\n');
  fprintf ('Credible interval: %.3g%%\n', 100 * prob);
end
fprintf ('Population mean within range of sample: %s\n', ...
         LogicalStr{mu_within_range_of_y + 1});
fprintf ('Population mean within range of CI (Frequentist): %s\n',...
         LogicalStr{mu_within_range_of_CI + 1});
fprintf ('Simulation size: %u\n', nsamp);
fprintf ('Simulation progress:       ');


% Perform simulation (in blocks of size nsim)
count = 0;
while (numel (mu_good) < nsamp)

  % Update counter
  count = count + 1;
  fprintf (repmat ('\b', 1, 7));
  fprintf ('% 6s%%', sprintf ('%.1f', round (1000 * numel (mu_good) / nsamp) / 10));

  % Sample some mean values from the (flat) prior within 
  % the range of the observed data (non-parametric)
  mu_sim = rand (1, nsim) * range + min(y);

  % Simulate data for each of these mean values
  Y = randn (n, nsim) * s;
  Y = bsxfun (@plus, Y, mu_sim);

  % Find simulated data that matches our "observed" data
  Y = sort(Y);
  y = sort(y);
  i = all (abs (bsxfun (@minus, Y, y)) <= eps);
  mu_good = cat (2, mu_good, mu_sim(i));

end

% Calculate statistics for Bayesian inference
mu_good_within_CI = (CI(1) <= mu_good(1:nsamp)) & (mu_good(1:nsamp) <= CI(2));

% Print results
% Probability that population mean lies within the CI (Bayesian)'
fprintf (repmat ('\b', 1, 7));
fprintf ('% 6s%%', sprintf ('%.1f', 100));
fprintf (['\nGiven our prior belief that the population mean could be equally\n',...
          'likely within the range of the data, our updated belief (posterior)\n',...
          'is summarized as a credible interval within which the population\n',...
          'mean effectively exists with %.2f%% probability.\n'],...
          sum (mu_good_within_CI) / nsamp * 100);

