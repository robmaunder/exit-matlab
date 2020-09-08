% Soft QPSK demodulator using natural mapping
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

% apriori_llrs is a 1xk vector of a priori LLRs
% rx is a complex symbol
% channel is a complex channel coefficient
% N0 is the noise power spectral density
% extrinsic_llrs is a 1xk vector of extrinsic LLRs
function extrinsic_llrs = soft_demodulate(rx, channel, N0, apriori_llrs)

    % Specify the constellation points and the bit mapping here
    constellation_points = [+1+1i; -1+1i; -1-1i; +1-1i]/sqrt(2);
    bit_labels = [0,0; 0,1; 1,0; 1,1];
    
    % Determine the number of bits per symbol and the number of constellation points here
    k = size(bit_labels,2);
    M = 2^k;
    N = length(rx);

    % Check that all the vectors and matrices have the correct dimensions
    if ~isequal(size(constellation_points),[M,1]) || ~isequal(size(bit_labels),[M,k])
        error('wrong dimensions');
    end
    
    if length(channel) ~= length(rx) && length(channel) ~= 1
        error('wrong dimensions');
    end
    
    
    
aposteriori_symbol_LLRs = zeros(M,N);

% Put the influence of the received signals into the symbol LLRs
for perm_index = 1:M
    aposteriori_symbol_LLRs(perm_index,:) = -abs(rx-channel*constellation_points(perm_index)).^2./N0;
end

if exist('apriori_llrs','var')
    % Put the influence of the apriori LLRs into the symbol LLRs
    for bit_index = 1:k       
%        aposteriori_symbol_LLRs(:,bit_permutations(bit_index,:) == 0) = aposteriori_symbol_LLRs(:,bit_permutations(bit_index,:) == 0) + repmat(apriori_llrs(bit_index:bits_per_symbol:end),[1,permutations/2]);
        
        for perm_index = 1:M
            if bit_labels(perm_index,bit_index) == 0
                aposteriori_symbol_LLRs(perm_index,:) = aposteriori_symbol_LLRs(perm_index,:) + apriori_llrs(bit_index:k:end);
            end
        end   
    end
end

% Extract the aposteriori LLRs from the symbol LLRs
aposteriori_llrs = zeros(1,N*k);
for bit_index = 1:k
    p0 = -inf(1,N);
    p1 = -inf(1,N);

    for perm_index = 1:M
        if bit_labels(perm_index,bit_index) == 0
            p0 = jac(p0, aposteriori_symbol_LLRs(perm_index,:));
        else
            p1 = jac(p1, aposteriori_symbol_LLRs(perm_index,:));
        end
    end   
    
    aposteriori_llrs(bit_index:k:end) = p0-p1;
end
        
if exist('apriori_llrs','var')
    % Remove the apriori from the aposteriori to get the extrinsic
    extrinsic_llrs = aposteriori_llrs - apriori_llrs;
else
    extrinsic_llrs = aposteriori_llrs;
end


    
 
end
