% BAYESIAN BOOTSTRAP CREDIBLE INTERVAL SIMULATION SCRIPT
clear

% bootbayes settings
nboot = 10000;
prob = 0.95;

% Simulation settings
n = 4;
mu = 0;
sigma = 1;
nsamp = 1e+05;                  % Total sample size desired from simulation
nsim = 1e+07;                   % Size of simulation block
eps = 5e-02 * sigma * sqrt (n); % Stringency/precision
mu_good = [];

% Create random sample
y = randn (n, 1) * sigma + mu; 
mu_within_range_of_y = (min (y) <= mu) & (mu <= max (y));
range = max (y) - min (y);
s = std (y, 1);        % Standard deviation (2nd moment, sample as a population)
if (exist ('pearsrnd'))
  skewflag = true;
  g = skewness (y, 1); % Skewness (3rd moment, sample as population)
else
  skewflag = false;
end

% Compute credible interval
stats = bootbayes (y, [], [], nboot, prob, 1);
CI = [stats.CI_lower stats.CI_upper];
mu_within_range_of_CI = (CI(1) <= mu) & (mu <= CI(2));  % Frequentist inference

% Print settings
LogicalStr = {'false', 'true'};
fprintf ('----- BAYESIAN BOOTSTRAP CREDIBLE INTERVAL (CI) SIMULATION -----\n')
fprintf ('Sample size: %u\n', n);
fprintf ('nboot: %u\n', nboot);
fprintf ('Prior: %.2f\n', 1.0);
if (numel(prob) > 1)
  fprintf ('Credible interval (CI) type: Percentile interval\n');
  mass =  100 * abs (prob(2) - prob(1));
  fprintf ('Credible interval: %.3g%% (%.1f%%, %.1f%%)\n', ...
           mass, 100 * prob);
else
  fprintf ('Credible interval (CI) type: Shortest probability interval\n');
  mass = 100 * prob;
  fprintf ('Credible interval: %.3g%%\n', mass);
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
  mu_sim = rand (1, nsim) * range + min (y);

  % Simulate data for each of these mean values with the same scale as the sample
  if (skewflag)
    Y = pearsrnd (0, s, g, 3, n, nsim); % 1st, 2nd and 3rd moment
  else
    Y = randn (n, nsim) * s;            % 1st and 2nd moment
  end
  Y = bsxfun (@plus, Y, mu_sim);

  % Find simulated data that matches our "observed" data
  Y = sort (Y);
  y = sort (y);
  i = all (abs (bsxfun (@minus, Y, y)) <= eps);
  mu_good = cat (2, mu_good, mu_sim(i));

end

% Calculate statistics for Bayesian inference
mu_good_within_CI = (CI(1) <= mu_good(1:nsamp)) & (mu_good(1:nsamp) <= CI(2));

% Print results
% Probability that population mean lies within the CI (Bayesian)'
fprintf (repmat ('\b', 1, 7));
fprintf ('% 6s%%', sprintf ('%.1f', 100));
fprintf (['\n Our prior belief is that the location of the population\n',...
          ' mean is equally likely anywhere within the range of the data.\n',...
          ' Our updated belief (i.e. posterior) after non-parametric\n',...
          ' Bayesian bootstrap is summarized as a credible interval within\n',...
          ' which the population mean exists with effectively %.2f%%\n',...
          ' probability (nominally %.1f%%).\n'],...
          sum (mu_good_within_CI) / nsamp * 100, mass);

