% Plots the Continuous-input Continuous-output Memoryless Channel (CCMC)
% and Discrete-input Continuous-output Memoryless Channel (DCMC) capacity
% of AWGN and uncorrelated Rayleigh fading channels for BPSK, QPSK, 8PSK
% and 16QAM. Copyright (C) 2011  Robert G. Maunder

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

% +---------------------------------------------+
% | Choose the SNR, modulation and channel here |
% +---------------------------------------------+

% Control the accuracy and duration of the simulation
symbol_count = 10000;

% Set range of channel SNRs
snr = -10:0.1:30; % dB

% Select channel
channel = 'AWGN';
% channel = 'Rayleigh';


% Seelct modulation scheme
% BPSK
%modulation = [+1, -1];

% QPSK
modulation = [+1, +1i, -1, -1i];

% 8PSK
% modulation = [+1, sqrt(1/2)*(+1+1i), +1i, sqrt(1/2)*(-1+1i), -1, sqrt(1/2)*(-1-1i), -1i, sqrt(1/2)*(+1-1i)];

% 16QAM
% modulation = sqrt(1/10)*[-3+3*1i, -1+3*1i, +1+3*1i, +3+3*1i, -3+1*1i, -1+1*1i, +1+1*1i, +3+1*1i, -3-1*1i, -1-1*1i, +1-1*1i, +3-1*1i, -3-3*1i, -1-3*1i, +1-3*1i, +3-3*1i];

% If you add more modulation schemes here, make sure their average transmit power is normalised to unity



% +------------------------+
% | Simulation starts here |
% +------------------------+

DCMC_capacity = nan(size(snr));
CCMC_capacity = nan(size(snr));

for index = 1:length(snr)

    % Generate some random symbols
    symbols = ceil(length(modulation)*rand(1,symbol_count));

    % Generate the transmitted signal
    x = modulation(symbols);

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
    N0 = 1/(10^(snr(index)/10));
    n = sqrt(N0/2)*(randn(1,symbol_count)+1i*randn(1,symbol_count));

    % Generate the received signal
    y = x.*h+n;

    % Calculate the symbol probabilities
    probabilities = max(exp(-(abs(ones(length(modulation),1)*y - modulation.'*h).^2)/N0),realmin);

    % Normalise the symbol probabilities
    probabilities = probabilities ./ (ones(length(modulation),1)*sum(probabilities));

    % Calculate the DCMC capacity
    DCMC_capacity(index) = log2(length(modulation))+mean(sum(probabilities.*log2(probabilities)));

    % Calculate the CCMC capacity
    if strcmp(channel, 'AWGN')
        CCMC_capacity(index)=log2(1+(10^(snr(index)/10)));
    elseif strcmp(channel, 'Rayleigh')
        % Symbolic Math Toolbox required to calculate CCMC capacity for
        % Rayleigh fading channel
        try
            CCMC_capacity(index) = log2(exp(1)) * exp (-1/(10^(snr(index)/10))) * (-vpa(eulergamma) + log (10^(snr(index)/10)) + 1/(10^(snr(index)/10)) ) ; %eq 8 from https://ieeexplore.ieee.org/stamp/stamp.jsp?tp=&arnumber=130999
        catch
    
        end
    else
        error('Unsupported channel');
    end
end

figure
plot(snr,DCMC_capacity);
hold on
plot(snr,CCMC_capacity);
xlabel('SNR [dB]');
ylabel('Capacity [bit/s/Hz]')

figure
plot(snr-10*log10(DCMC_capacity),DCMC_capacity);
hold on
plot(snr-10*log10(CCMC_capacity),CCMC_capacity);
xlabel('E_b/N_0 [dB]');
ylabel('Capacity [bit/s/Hz]')

