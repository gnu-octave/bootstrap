function r = autocorr (x, maxlag)

  % Private function file required for ibootci

  % Efficient calculation of autocorrelation by fast fourier transform
  [m,n] = size(x);
  if nargin > 1
    if ~isa(maxlag,'numeric') || maxlag < 1 || maxlag ~= abs(maxlag)
      error('maxlag must be an integer >= 1')
    end
  else
    maxlag = m;
  end
  if n > 1
    error('x must be a vector')
  end
  X = fft(x,2^nextpow2(2*m-1));
  r = ifft(abs(X).^2);
  r(min(m,maxlag)+1:end) = [];
  r = r./r(1);

end
