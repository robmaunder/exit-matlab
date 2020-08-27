% Encoder for a half-rate systematic recursive convolutional code
% having 1 memory element, a generator polynomial of [1,0] and a feedback
% polynomial of [1,1].
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



function [encoded1_bits, encoded2_bits] = convolutional_encoder(uncoded_bits)

    % Systematic bits
    encoded1_bits = uncoded_bits;
    
    % Parity bits
    encoded2_bits=zeros(1,length(uncoded_bits));
    encoded2_bits(1) = uncoded_bits(1);
    for i = 2:length(uncoded_bits)
        encoded2_bits(i) = mod(encoded2_bits(i-1)+uncoded_bits(i),2);
    end

end