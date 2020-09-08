% Jacobian logarithm
% If A = log(a) and B = log(b), then log(a+b) = max(A,B) + log(1+exp(-abs(A-B)))
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



function C = jac(A,B)

    C = max(A,B) + log(1+exp(-abs(A-B)));
%    C = max(A,B);
    


end
	