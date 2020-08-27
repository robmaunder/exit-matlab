% Generate Gaussian distributed a priori LLRs
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



function llrs = generate_llrs(bits, mutual_information)
    if(mutual_information < 0 || mutual_information >= 1)
        error('mutual_information must be in the range [0,1)');
    end
	sigma = (-1.0/0.3073*log(1.0-mutual_information^(1.0/1.1064))/log(2.0))^(1.0/(2.0*0.8935));
	llrs = randn(1,length(bits))*sigma - (bits-0.5)*sigma^2;
end