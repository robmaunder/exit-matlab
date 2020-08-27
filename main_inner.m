% EXIT function for a convolutional code used as an inner code
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

% Channel SNR in dB
SNR = -6;

% Noise variance
N0 = 1/10^(SNR/10);

% Generate some random bits
uncoded_bits  = round(rand(1,bit_count));

% Encode using a half-rate systematic recursive convolutional code having a single memory element
[encoded1_bits, encoded2_bits] = convolutional_encoder(uncoded_bits);

% BPSK modulator
tx1 = -2*(encoded1_bits-0.5);
tx2 = -2*(encoded2_bits-0.5);

% Send the two BPSK signals one at a time over an AWGN channel
rx1 = tx1 + sqrt(N0/2)*(randn(1,length(tx1))+1i*randn(1,length(tx1)));
rx2 = tx2 + sqrt(N0/2)*(randn(1,length(tx2))+1i*randn(1,length(tx2)));

% BPSK demodulator
apriori_encoded1_llrs = (abs(rx1+1).^2-abs(rx1-1).^2)/N0;
apriori_encoded2_llrs = (abs(rx2+1).^2-abs(rx2-1).^2)/N0;

% Plot the LLR histograms
display_llr_histograms([apriori_encoded1_llrs,apriori_encoded2_llrs],[encoded1_bits,encoded2_bits]);

% A priori mutual informations to consider
IA = 0.999*(0:1/(IA_count-1):1);

% Initialise results
IE_av=zeros(1,IA_count);
IE_hist=zeros(1,IA_count);
area=0.0;

% Consider each a priori mutual information
for IA_index = 1:IA_count

    % Generate the a priori LLRs having the a priori mutual information considered
    apriori_uncoded_llrs = generate_llrs(uncoded_bits, IA(IA_index));
  
    % Do the BCJR
    [aposteriori_uncoded_llrs, aposteriori_encoded1_llrs, aposteriori_encoded2_llrs] = bcjr_decoder(apriori_uncoded_llrs, apriori_encoded1_llrs, apriori_encoded2_llrs);

    % Calculate the new information
    extrinsic_uncoded_llrs = aposteriori_uncoded_llrs-apriori_uncoded_llrs;

    % Measure the mutual information of the extrinsic LLRs
    IE_hist(IA_index) = measure_mutual_information_histogram(extrinsic_uncoded_llrs, uncoded_bits);
    IE_av(IA_index) = measure_mutual_information_averaging(extrinsic_uncoded_llrs);

    % Update the area beneath the EXIT function
    if(IA_index > 1)
       area = area + (IE_av(IA_index)+IE_av(IA_index-1))*(IA(IA_index)-IA(IA_index-1))/2;
    end
end


% Plot EXIT function
figure
xlim([0 1]);
ylim([0 1]);
xlabel('Quality of input LLRs (a priori mutual information I_A)');
ylabel('Quality of output LLRs (extrinsic mutual information I_E)');
title(['EXIT function for SNR = ', num2str(SNR), ' dB']);
hold on
plot(IA,IE_hist,'r');
plot(IA,IE_av,'b');
legend({'True quality','Claimed quality'},'Location','northwest');

% Display the area beneath the EXIT function
annotation('textbox','String',{['Area = ', num2str(area)]},'LineStyle','none','Position',[0.7 0.1 0.2 0.1]);

