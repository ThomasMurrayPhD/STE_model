function [pstruct] = obs_suttonK1_comb_obs_namep(pvec)
% [pstruct] = m1_comb_obs_namep(pvec)
%
% Creates a struct with a field for each parameter from parameter vector.
% (Designed to be compatible with the HGF Toolbox as part of TAPAS).
%
% INPUT
%   pvec       vector        1xP vector containing model params
%
%   OPTIONAL:
%
% OUTPUT    
%   pstruct      struct      Struct with fields for mdoel params
%
% _________________________________________________________________________
% Author: Alex Hess
%
% Copyright (C) 2023 Translational Neuromodeling Unit
%                    Institute for Biomedical Engineering
%                    University of Zurich & ETH Zurich
%
% This file is released under the terms of the GNU General Public Licence
% (GPL), version 3. You can redistribute it and/or modify it under the
% terms of the GNU General Public License as published by the Free Software
% Foundation, either version 3 of the License, or (at your option) any
% later version.
%
% This file is distributed in the hope that it will be useful, but WITHOUT
% ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
% FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for
% more details.
% 
% You should have received a copy of the GNU General Public License along
% with this program. If not, see <https://www.gnu.org/licenses/>.
% _________________________________________________________________________

pstruct = struct;

pstruct.ze = pvec(1);
pstruct.beta0 = pvec(2);
pstruct.beta1 = pvec(3);
pstruct.beta2 = pvec(4);
pstruct.sa = pvec(length(pvec));

return;
