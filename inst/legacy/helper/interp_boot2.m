function  U = interp_boot2 (T2, T0, C)

  % Helper function file required for ibootci

    U = sum(T2<=T0);
    if U < 1
      U = 0;
    elseif U == C
      U = C;
    else
      % Quick linear interpolation to approximate asymptotic calibration
      t2 = zeros(1,2);
      I = (T2<=T0);
      if any(I)
        t2(1) = max(T2(I));
      else
        t2(1) = min(T2);
      end
      I = (T2>T0);
      if any(I)
        t2(2) = min(T2(I));
      else
        t2(2) = max(T2);
      end
      if (t2(2)-t2(1) == 0)
        U = t2(1);
      else
        U = ((t2(2)-T0)*U + (T0-t2(1))*(U+1)) /...
                (t2(2) - t2(1));
      end
    end
    U = U/C;

end
