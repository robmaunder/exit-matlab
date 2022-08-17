% Plots the Continuous-input Continuous-output Memoryless Channel (CCMC)
% and Discrete-input Continuous-output Memoryless Channel (DCMC) capacity
% of AWGN and uncorrelated Rayleigh fading channels for BPSK, QPSK, 8PSK
% and 16QAM.
% Copyright (C) 2011  Robert G. Maunder

% This program is free software: you can redistribute it and/or modify it
% under the terms of the GNU General Public License as published by the
% Free Software Foundation, either version 3 of the License, or (at your
% option) any later version.

% This program is distributed in the hope that it will be useful, but
% WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General
% Public License for more details.

% The GNU General Public License can be seen at http://www.gnu.org/licenses/.

clear all;
close all;

% Control the accuracy and duration of the simulation
symbol_count = 10000;

% Set range of channel SNRs
snr = -10:0.1:30; % dB

% Select channel
channel = 'AWGN';
% channel = 'Rayleigh';

% Setup modulation schemes
modulation_name{1} = 'BPSK';
modulation{1} = [+1, -1];

modulation_name{2} = 'QPSK';
modulation{2} = [+1, +1i, -1, -1i];

modulation_name{3} = '8PSK';
modulation{3} = [+1, sqrt(1/2)*(+1+1i), +1i, sqrt(1/2)*(-1+1i), -1, sqrt(1/2)*(-1-1i), -1i, sqrt(1/2)*(+1-1i)];

modulation_name{4} = '16QAM';
modulation{4} = sqrt(1/10)*[-3+3*1i, -1+3*1i, +1+3*1i, +3+3*1i, -3+1*1i, -1+1*1i, +1+1*1i, +3+1*1i, -3-1*1i, -1-1*1i, +1-1*1i, +3-1*1i, -3-3*1i, -1-3*1i, +1-3*1i, +3-3*1i];

% If you add more modulation schemes here, make sure their average transmit power is normalised to unity



% Calculate the CCMC capacity
CCMC_capacity = nan(size(snr));
for snr_index = 1:length(snr)
    Gamma = 10^(snr(snr_index)/10);
    if strcmp(channel, 'AWGN')
        CCMC_capacity(snr_index)=log2(1+Gamma);
    elseif strcmp(channel, 'Rayleigh')
        % Refer to the following paper for CCMC capacity of Rayleigh fading
        % channel
        % https://ieeexplore.ieee.org/stamp/stamp.jsp?tp=&arnumber=130999
    else
        error('Unsupported channel');
    end
end

% Plot vs SNR
my_legend = {'CCMC'};
figure(1);
plot(snr,CCMC_capacity);
xlabel('SNR [dB]');
ylabel('Capacity [bit/s/Hz]')
legend(my_legend);
title(channel);
hold on

% Plot vs Eb/N0
figure(2);
plot(snr-10*log10(CCMC_capacity),CCMC_capacity);
xlabel('E_b/N_0 [dB]');
ylabel('Capacity [bit/s/Hz]')
legend(my_legend);
title(channel);
hold on

% DCMC simulation
for modulation_index = 1:length(modulation)
    DCMC_capacity = nan(size(snr));
    for snr_index = 1:length(snr)

        % Generate some random symbols
        symbols = ceil(length(modulation{modulation_index})*rand(1,symbol_count));

        % Generate the transmitted signal
        x = modulation{modulation_index}(symbols);

        % Generate the channel gains
        if strcmp(channel, 'AWGN')
            h = ones(1,symbol_count);
        elseif strcmp(channel, 'Rayleigh')
            % Uncorrelated narrowband Rayleigh fading channel
            h = sqrt(1/2)*(randn(1,symbol_count)+1i*randn(1,symbol_count));
        else
            error('Unsupported channel');
        end

        % Generate some noise
        N0 = 1/(10^(snr(snr_index)/10));
        n = sqrt(N0/2)*(randn(1,symbol_count)+1i*randn(1,symbol_count));

        % Generate the received signal
        y = x.*h+n;

        % Calculate the symbol probabilities
        probabilities = max(exp(-(abs(ones(length(modulation{modulation_index}),1)*y - modulation{modulation_index}.'*h).^2)/N0),realmin);

        % Normalise the symbol probabilities
        probabilities = probabilities ./ (ones(length(modulation{modulation_index}),1)*sum(probabilities));

        % Calculate the DCMC capacity
        DCMC_capacity(snr_index) = log2(length(modulation{modulation_index}))+mean(sum(probabilities.*log2(probabilities)));
    end

    % Plot vs SNR
    my_legend{end+1} = modulation_name{modulation_index};
    figure(1)
    plot(snr,DCMC_capacity);
    legend(my_legend);

    % Plot vs Eb/N0
    figure(2)
    plot(snr-10*log10(DCMC_capacity),DCMC_capacity);
    legend(my_legend);

end