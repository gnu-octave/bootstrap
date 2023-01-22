% Package error checking

info = ver; 
ISOCTAVE = any (ismember ({info.Name}, 'Octave'));

% Turn off some warnings to keep tests tidy
warning ('off', 'bootknife:parallel')
%if ISOCTAVE
%  warning ('off', 'Octave:divide-by-zero')
%  warning ('off', 'Octave:nearly-singular-matrix')
%  warning ('off', 'Octave:broadcast')
%else
%  warning ('off', 'MATLAB:rankDeficientMatrix')
%  warning ('off', 'MATLAB:divideByZero')
%end

try 
  % boot
  boot (3, 20);
  boot (3, 20, false, 1);
  boot (3, 20, true, 1);
  boot (3, 20, [], 1);
  boot (3, 20, true, 1, [30,30,0]);
  
  % bootknife 
  % bootknife:test:1
  y = randn (20,1); 
  strata = [1;1;1;1;1;1;1;1;1;1;2;2;2;2;2;3;3;3;3;3];
  stats = bootknife (y, 2000, @mean);
  stats = bootknife (y, 2000, 'mean');
  stats = bootknife (y, 2000, {@var,1});
  stats = bootknife (y, 2000, {'var',1});
  stats = bootknife (y, 2000, @mean, [], strata);
  stats = bootknife (y, 2000, {'var',1}, [], strata);
  stats = bootknife (y, 2000, {@var,1}, [], strata, 2);
  stats = bootknife (y, 2000, @mean, .1, strata, 2);
  stats = bootknife (y, 2000, @mean, [.05,.95], strata, 2);
  stats = bootknife (y, [2000,200], @mean, .1, strata, 2);
  stats = bootknife (y, [2000,200], @mean, [.05,.95], strata, 2);
  stats = bootknife (y(1:5), 2000, @mean, .1);
  stats = bootknife (y(1:5), 2000, @mean, [.05,.95]);
  stats = bootknife (y(1:5), [2000,200], @mean, .1);
  stats = bootknife (y(1:5), [2000,200], @mean, [.05,.95]);
  % bootknife:test:2
  Y = randn (20); 
  strata = [1;1;1;1;1;1;1;1;1;1;2;2;2;2;2;3;3;3;3;3];
  stats = bootknife (Y, 2000, @mean);
  stats = bootknife (Y, 2000, 'mean');
  stats = bootknife (Y, 2000, {@var, 1});
  stats = bootknife (Y, 2000, {'var',1});
  stats = bootknife (Y, 2000, @mean, [], strata);
  stats = bootknife (Y, 2000, {'var',1}, [], strata);
  stats = bootknife (Y, 2000, {@var,1}, [], strata, 2);
  stats = bootknife (Y, 2000, @mean, .1, strata, 2);
  stats = bootknife (Y, 2000, @mean, [.05,.95], strata, 2);
  stats = bootknife (Y, [2000,200], @mean, .1, strata, 2);
  stats = bootknife (Y, [2000,200], @mean, [.05,.95], strata, 2);
  stats = bootknife (Y(1:5,:), 2000, @mean, .1);
  stats = bootknife (Y(1:5,:), 2000, @mean, [.05,.95]);
  stats = bootknife (Y(1:5,:), [2000,200], @mean, .1);
  stats = bootknife (Y(1:5,:), [2000,200], @mean, [.05,.95]);
  stats = bootknife (Y, 2000, @(Y) mean(Y(:),1)); % Cluster/block resampling
  % Y(1,end) = NaN; % Unequal clustersize
  %stats = bootknife (Y, 2000, @(Y) mean(Y(:),1,'omitnan'));
  % bootknife:test:3
  y = randn (20,1); x = randn (20,1); X = [ones(20,1), x];
  stats = bootknife ({x,y}, 2000, @cor);
  stats = bootknife ({x,y}, 2000, @cor, [], strata);
  stats = bootknife ({y,x}, 2000, @(y,x) pinv(x)*y); % Could also use @regress
  stats = bootknife ({y,X}, 2000, @(y,X) pinv(X)*y);
  stats = bootknife ({y,X}, 2000, @(y,X) pinv(X)*y, [], strata);
  stats = bootknife ({y,X}, 2000, @(y,X) pinv(X)*y, [], strata, 2);
  stats = bootknife ({y,X}, 2000, @(y,X) pinv(X)*y, [.05,.95], strata);
  
  % bootci
  % bootci:test:1
  y = randn (20, 1); 
  bootci (2000, 'mean', y);
  bootci (2000, @mean, y);
  bootci (2000, @mean, y, 'alpha', 0.1);
  bootci (2000, {'mean', y}, 'alpha', 0.1);
  bootci (2000, {@mean, y}, 'alpha', 0.1);
  bootci (2000, {@mean, y}, 'alpha', 0.1, 'seed', 1);
  bootci (2000, {@mean, y}, 'alpha', 0.1, 'type', 'norm');
  bootci (2000, {@mean, y}, 'alpha', 0.1, 'type', 'per');
  bootci (2000, {@mean, y}, 'alpha', 0.1, 'type', 'basic');
  bootci (2000, {@mean, y}, 'alpha', 0.1, 'type', 'bca');
  bootci (2000, {@mean, y}, 'alpha', 0.1, 'type', 'stud');
  bootci (2000, {@mean, y}, 'alpha', 0.1, 'type', 'stud', 'nbootstd', 100);
  bootci (2000, {@mean, y}, 'alpha', 0.1, 'type', 'cal');
  bootci (2000, {@mean, y}, 'alpha', 0.1, 'type', 'cal', 'nbootcal', 200);
  % bootci:test:2
  Y = randn (20); 
  bootci (2000, 'mean', Y);
  bootci (2000, @mean, Y);
  bootci (2000, @mean, Y, 'alpha', 0.1);
  bootci (2000, {'mean', Y}, 'alpha', 0.1);
  bootci (2000, {@mean, Y}, 'alpha', 0.1);
  bootci (2000, {@mean, Y}, 'alpha', 0.1, 'seed', 1);
  bootci (2000, {@mean, Y}, 'alpha', 0.1, 'type', 'norm');
  bootci (2000, {@mean, Y}, 'alpha', 0.1, 'type', 'per');
  bootci (2000, {@mean, Y}, 'alpha', 0.1, 'type', 'basic');
  bootci (2000, {@mean, Y}, 'alpha', 0.1, 'type', 'bca');
  bootci (2000, {@mean, Y}, 'alpha', 0.1, 'type', 'stud');
  bootci (2000, {@mean, Y}, 'alpha', 0.1, 'type', 'cal');
  % bootci:test:3
  y = randn (20,1); x = randn (20,1); X = [ones(20,1),x];
  bootci (2000, @cor, x, y);
  bootci (2000, @(y,X) pinv(X)*y, y, X);
  bootci (2000, @(y,X) pinv(X)*y, y, X, 'alpha', 0.1);
  bootci (2000, {@(y,X) pinv(X)*y, y, X}, 'alpha', 0.1);
  bootci (2000, {@(y,X) pinv(X)*y, y, X}, 'alpha', 0.1, 'type', 'norm');
  bootci (2000, {@(y,X) pinv(X)*y, y, X}, 'alpha', 0.1, 'type', 'per');
  bootci (2000, {@(y,X) pinv(X)*y, y, X}, 'alpha', 0.1, 'type', 'basic');
  bootci (2000, {@(y,X) pinv(X)*y, y, X}, 'alpha', 0.1, 'type', 'bca');
  bootci (2000, {@(y,X) pinv(X)*y, y, X}, 'alpha', 0.1, 'type', 'stud');
  bootci (2000, {@(y,X) pinv(X)*y, y, X}, 'alpha', 0.1, 'type', 'cal');
  
  % bootstrp
  y = randn (20,1);
  bootstat = bootstrp (50, @mean, y);
  
  % bootnhst 
  % bootnhst:test:1
  y = [111.39 110.21  89.21  76.64  95.35  90.97  62.78;
       112.93  60.36  92.29  59.54  98.93  97.03  79.65;
        85.24 109.63  64.93  75.69  95.28  57.41  75.83;
       111.96 103.40  75.49  76.69  77.95  93.32  78.70];
  g = [1 2 3 4 5 6 7;
       1 2 3 4 5 6 7;
       1 2 3 4 5 6 7;
       1 2 3 4 5 6 7];
  p = bootnhst (y(:),g(:),'ref',1,'nboot',[1000,0],'DisplayOpt',false);
  p = bootnhst (y(:),g(:),'nboot',[1000,0],'DisplayOpt',false);
  % bootnhst:test:2
  y = [54       43
       23       34 
       45       65
       54       77
       45       46
      NaN       65];
  g = {'male' 'female'
       'male' 'female'
       'male' 'female'
       'male' 'female'
       'male' 'female'
       'male' 'female'};
  p = bootnhst (y(:),g(:),'ref','male','nboot',[1000,0],'DisplayOpt',false);
  p = bootnhst (y(:),g(:),'nboot',[1000,0],'DisplayOpt',false);
  % bootnhst:test:3
  y = [54  87  45
       23  98  39
       45  64  51
       54  77  49
       45  89  50
       47 NaN  55];
  g = [ 1   2   3
        1   2   3
        1   2   3
        1   2   3
        1   2   3
        1   2   3];
  p = bootnhst (y(:),g(:),'nboot',[1000,0],'DisplayOpt',false);
  p = bootnhst (y(:),g(:),'bootfun',@(y)std(y,1),'DisplayOpt',false);
  p = bootnhst (y(:),g(:),'bootfun',{@std,1},'DisplayOpt',false);
  % bootnhst:test:4
  Y = randn (20, 2); g = [zeros(10, 1); ones(10, 1)];
  func = @(M) cor (M(:,1), M(:,2));
  p = bootnhst (Y, g, 'bootfun', func, 'DisplayOpt', false);
  % bootnhst:test:5
  y = randn (20, 1); x = randn (20, 1); X = [ones(20, 1), x];
  g = [zeros(10, 1); ones(10, 1)];
  func = @(M) subsref (M(:,2:end) \ M(:,1), struct ('type', '()', 'subs', {{2}}));
  p = bootnhst ([y, X], g, 'bootfun', func, 'DisplayOpt', false);

catch exception

  % Turn warnings back on 
  warning ('on', 'bootknife:parallel')
  %if ISOCTAVE
  %  warning ('on', 'Octave:divide-by-zero')
  %  warning ('on', 'Octave:nearly-singular-matrix')
  %else
  %  warning ('on', 'MATLAB:rankDeficientMatrix')
  %  warning ('on', 'MATLAB:divideByZero')
  %end
  %rethrow (exception)


end

% Turn warnings back on 
warning ('on', 'bootknife:parallel');
%if ISOCTAVE
%  warning ('on', 'Octave:divide-by-zero')
%  warning ('on', 'Octave:nearly-singular-matrix')
%  warning ('on', 'Octave:broadcast')
%else
%  warning ('on', 'MATLAB:rankDeficientMatrix')
%  warning ('on', 'MATLAB:divideByZero')
%end
