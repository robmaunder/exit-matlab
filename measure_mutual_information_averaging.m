% Measure the mutual information of some LLRs using the averaging method.
% This method assumes that the LLRs are self-consistent and
% well-conditioned.
% Copyright (C) 2008  Robert G. Maunder

% This program is free software: you can redistribute it and/or modify it 
% under the terms of the GNU General Public License as published by the
% Free Software Foundation, either version 3 of the License, or (at your 
% option) any later version.

% This program is distributed in the hope that it will be useful, but 
% WITHOUT ANY WARRANTY; without even the implied warranty of 
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General 
% Public License for more details.

% The GNU General Public License can be seen at http://www.gnu.org/licenses/.



% llrs is a 1xK vector of LLRs
% mutual_information is a scalar in the range 0 to 1
function mutual_information = measure_mutual_information_averaging(llrs)
    P0 = exp(llrs)./(1+exp(llrs));
    P1 = 1-P0;
    entropies = -P0.*log2(P0)-P1.*log2(P1);
    mutual_information = 1-sum(entropies(~isnan(entropies)))/length(entropies);
end