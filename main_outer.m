% EXIT function for a convolutional code used as an outer code
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



% Number of bits to encode
bit_count=100000;

% Number of a priori mutual informations to consider
IA_count=11;

% Generate some random bits
uncoded_bits  = round(rand(1,bit_count));

% Encode using a half-rate systematic recursive convolutional code having a single memory element
[encoded1_bits, encoded2_bits] = convolutional_encoder(uncoded_bits);

% A priori mutual informations to consider
IA = 0.999*(0:1/(IA_count-1):1);

% Initialise results
IE_hist=zeros(1,IA_count);
IE_av=zeros(1,IA_count);
BER=zeros(1,IA_count);
area=0.0;

% Consider each a priori mutual information
for IA_index = 1:IA_count

    % Generate the a priori LLRs having the a priori mutual information considered
    apriori_encoded1_llrs = generate_llrs(encoded1_bits, IA(IA_index));
    apriori_encoded2_llrs = generate_llrs(encoded2_bits, IA(IA_index));

    % No a priori information for the uncoded bits when operating as an outer code
    apriori_uncoded_llrs = zeros(1,length(uncoded_bits));

    % Do the BCJR
    [aposteriori_uncoded_llrs, aposteriori_encoded1_llrs, aposteriori_encoded2_llrs] = bcjr_decoder(apriori_uncoded_llrs, apriori_encoded1_llrs, apriori_encoded2_llrs);

    % Calculate the new information
    extrinsic_encoded1_llrs = aposteriori_encoded1_llrs-apriori_encoded1_llrs;
    extrinsic_encoded2_llrs = aposteriori_encoded2_llrs-apriori_encoded2_llrs;

    % Measure the mutual information of the extrinsic LLRs
    IE_hist(IA_index) = (measure_mutual_information_histogram(extrinsic_encoded1_llrs, encoded1_bits) + measure_mutual_information_histogram(extrinsic_encoded2_llrs, encoded2_bits))/2;
    IE_av(IA_index) = (measure_mutual_information_averaging(extrinsic_encoded1_llrs) + measure_mutual_information_averaging(extrinsic_encoded2_llrs))/2;
    
    % Calculate the BER
    decoded_bits = aposteriori_uncoded_llrs < 0;
    BER(IA_index) = sum(uncoded_bits ~= decoded_bits)/length(uncoded_bits);

    % Update the area beneath the EXIT function
    if(IA_index > 1)
       area = area + (IE_av(IA_index)+IE_av(IA_index-1))*(IA(IA_index)-IA(IA_index-1))/2;
    end
    
end

% Plot BER
figure
semilogy(IA,BER);
xlim([0 1]);
ylim([min(100/bit_count,0.1) 1]);
xlabel('Quality of input LLRs (a priori mutual information I_A)');
ylabel('BER');

% Plot inverted EXIT function
figure
xlim([0 1]);
ylim([0 1]);
xlabel('Quality of output LLRs (extrinsic mutual information I_E)');
ylabel('Quality of input LLRs (a priori mutual information I_A)');
title('Inverted EXIT function');
hold on
plot(IE_hist,IA,'r');
plot(IE_av,IA,'b');
legend({'True quality','Claimed quality'},'Location','northwest');

% Display the area beneath the inverted EXIT function
annotation('textbox','String',{['Area = ', num2str(1-area)]},'LineStyle','none','Position',[0.7 0.1 0.2 0.1]);

