% EXIT function for a soft-input soft-output demodulator
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
SNR = 0;

% Noise variance
N0 = 1/10^(SNR/10);

% A priori mutual informations to consider
IA = 0.999*(0:1/(IA_count-1):1);

% Initialise results
IE_av=zeros(1,IA_count);
IE_hist=zeros(1,IA_count);
area=0.0;

% Consider each a priori mutual information
for IA_index = 1:IA_count

    % Generate some random bits
    bits  = round(rand(1,bit_count));

    % Encode using a half-rate systematic recursive convolutional code having a single memory element
    tx = modulate(bits);

    % Rayleigh fading 
    h = sqrt(1/2)*(randn(size(tx))+1i*randn(size(tx)));

    % Noise
    n = sqrt(N0/2)*(randn(size(tx))+1i*randn(size(tx)));

    % Uncorrelated narrowband Rayleigh fading channel
    rx = h.*tx + n;


    % Generate the a priori LLRs having the a priori mutual information considered
    apriori_llrs = generate_llrs(bits, IA(IA_index));
  
    % Do the BCJR
    extrinsic_llrs = soft_demodulate(rx, h, N0, apriori_llrs);

    % Measure the mutual information of the extrinsic LLRs
    IE_hist(IA_index) = measure_mutual_information_histogram(extrinsic_llrs, bits);
    IE_av(IA_index) = measure_mutual_information_averaging(extrinsic_llrs);

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

