% -- Function File: L = trt_vs_ctrl (L_EMM)
% -- Function File: L = trt_vs_ctrl (L_EMM, REF)
% -- Function File: [L, PAIRS] = trt_vs_ctrl (...)
%
%     Function file to create a hypothesis matrix for computing comparisons
%     of treatment versus control from the regression coefficients of a linear
%     model given the hypothesis matrix used to compute the estimated marginal
%     means. By default, the control is set to group 1 corresponding to the
%     first row in L_EMM. This setting can be changed using the REF input
%     argument. To specify an alternative control The function can also return
%     PAIRS, which is a 2-column matrix defining which groups (i.e. rows of
%     L_EMM) are being compared in each column of L.
%
%     The syntax in this function code is known to be compatible with
%     recent versions of Octave (v3.2.4 on Debian 6 Linux 2.6.32) and
%     Matlab (v6.5.0 and v7.4.0 on Windows XP).
%
%  pairwise (version 2023.07.04)
%  Author: Andrew Charles Penn
%  https://www.researchgate.net/profile/Andrew_Penn/
%
%  Copyright 2019 Andrew Charles Penn
%  This program is free software: you can redistribute it and/or modify
%  it under the terms of the GNU General Public License as published by
%  the Free Software Foundation, either version 3 of the License, or
%  (at your option) any later version.
%
%  This program is distributed in the hope that it will be useful,
%  but WITHOUT ANY WARRANTY; without even the implied warranty of
%  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%  GNU General Public License for more details.
%
%  You should have received a copy of the GNU General Public License
%  along with this program.  If not, see <http://www.gnu.org/licenses/>.

function [L, pairs] = trt_vs_ctrl (L_EMM, REF)

  if (nargin < 2)
    REF = 1;
  end

  % Get number of group members from the hypothesis matrix used 
  % to generate estimated marginal means
  Ng = size (unique (L_EMM','rows'), 1);
  if (REF > Ng)
    error ('trt_vs_ctrl: REF exceeds number of groups (i.e. rows in L_EMM)');
  end

  % Create pairs matrix for pairwise comparisons
  gid = (1 : Ng)';  % Create numeric group ID
  pairs = zeros (Ng - 1, 2);
  pairs(:, 1) = REF;
  pairs(:, 2) = gid(gid ~= REF);

  % Calculate hypothesis matrix for pairwise comparisons from the
  % estimated marginal means
  Np = size (pairs, 1);
  L_PWC = zeros (Np, Ng);
  for j = 1:Np
    L_PWC(j, pairs(j,:)) = [1,-1];
  end
  
  % Create hypothesis matrix to generate pairwise comparisons directly
  % from the regression coefficients. Note that the contrasts used to
  % fit the original model must sum to zero
  L = (L_PWC * L_EMM')';

end