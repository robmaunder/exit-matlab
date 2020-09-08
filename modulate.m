% QPSK modulator using natural mapping
% Copyright (C) 2010  Robert G. Maunder

% This program is free software: you can redistribute it and/or modify it 
% under the terms of the GNU General Public License as published by the
% Free Software Foundation, either version 3 of the License, or (at your 
% option) any later version.

% This program is distributed in the hope that it will be useful, but 
% WITHOUT ANY WARRANTY; without even the implied warranty of 
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General 
% Public License for more details.

% The GNU General Public License can be seen at http://www.gnu.org/licenses/.

% bits is a k*n vector of bits
% tx is a vector of complex symbols
function tx = modulate(bits)

    % Specify the constellation points and the bit mapping here
    constellation_points = [+1+1i; -1+1i; -1-1i; +1-1i]/sqrt(2);
    bit_labels = [0,0; 0,1; 1,1; 1,0];
    
    % Determine the number of bits per symbol and the number of constellation points here
    k = size(bit_labels,2);
    M = 2^k;
    N = length(bits)/k;

    
    % Check that all the vectors and matrices have the correct dimensions
    if ~isequal(size(constellation_points),[M,1]) || ~isequal(size(bit_labels),[M,k])
        error('wrong dimensions');
    end

    symbols = bin2dec(num2str(reshape(bits,[k,N])'))'+1;
    tx = constellation_points(symbols);
    
    
end
    
