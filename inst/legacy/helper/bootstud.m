function  p = bootstud(m,bootstat,S)

  % Helper function file required for ibootp

  % Use bootstrap-t method with variance stabilization for small samples
  % Polansky (2000) Can J Stat. 28(3):501-516
  se = nanfun(@std,bootstat{1});
  if size(bootstat{2},1) > 1
    SE1 = nanfun(@std,bootstat{2});
  else
    SE1 = bootstat{2};
  end
  a = S.n(1)^(-3/2) * se;  % additive correction factor

  % Calculate Studentized statistics
  ridx = isnan(bootstat{1}); bootstat{1}(ridx)=[]; SE1(ridx)=[];
  T = (bootstat{1} - S.stat)./(SE1 + a);
  t = (S.stat - m)/se;

  % Calculate p value from empirical distribution of the Studentized bootstrap statistics
  [cdf,T] = empcdf(T,0);
  p = 1-interp1(T,cdf,t,'linear');

end

