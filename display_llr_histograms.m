% Display the histograms of LLRs. This can be used to check that the LLRs
% are self-consistent and well-conditioned.
% Copyright (C) 2009  Robert G. Maunder

% This program is free software: you can redistribute it and/or modify it
% under the terms of the GNU General Public License as published by the
% Free Software Foundation, either version 3 of the License, or (at your
% option) any later version.

% This program is distributed in the hope that it will be useful, but
% WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General
% Public License for more details.

% The GNU General Public License can be seen at http://www.gnu.org/licenses/.



% llrs is a 1xK vector of LLRs
% bits is a 1xK vector of the correct bit values
% bin_width is an optional input, which can be set to the difference
%   between consecutive LLR values in the case of fixed-point LLRs. If
%   bit_width is omitted, then an appropriate value is calculated
%   automatically.
function display_llr_histograms(llrs, bits, bin_width)


if(length(llrs) ~= length(bits))
    error('Must have same number of llrs and bits!');
end


bit_1_count = sum(bits);
bit_0_count = length(bits) - bit_1_count;
if(bit_0_count == 0 || bit_1_count == 0)
    error('All bits have the same value');
end

llr_0_noninfinite_count = 0;
llr_1_noninfinite_count = 0;
llr_0_max = -Inf;
llr_0_min = Inf;
llr_1_max = -Inf;
llr_1_min = Inf;
for bit_index = 1:length(bits)
    if(llrs(bit_index) ~= -Inf && llrs(bit_index) ~= Inf)
        if(bits(bit_index) == 0)
            llr_0_noninfinite_count = llr_0_noninfinite_count+1;
            
            if(llrs(bit_index) > llr_0_max)
                llr_0_max = llrs(bit_index);
            end
            if(llrs(bit_index) < llr_0_min)
                llr_0_min = llrs(bit_index);
            end
        else
            llr_1_noninfinite_count = llr_1_noninfinite_count+1;
            
            if(llrs(bit_index) > llr_1_max)
                llr_1_max = llrs(bit_index);
            end
            if(llrs(bit_index) < llr_1_min)
                llr_1_min = llrs(bit_index);
            end
        end
    end
end

if(llr_0_noninfinite_count > 0 && llr_1_noninfinite_count > 0)
    if ~exist('bin_width','var')
        llr_0_mean = 0.0;
        llr_1_mean = 0.0;
        for bit_index = 1:length(bits)
            if(llrs(bit_index) ~= -Inf && llrs(bit_index) ~= Inf)
                if(bits(bit_index) == 0)
                    llr_0_mean = llr_0_mean+llrs(bit_index);
                else
                    llr_1_mean = llr_1_mean+llrs(bit_index);
                end
            end
            
        end
        llr_0_mean = llr_0_mean/llr_0_noninfinite_count;
        llr_1_mean = llr_1_mean/llr_1_noninfinite_count;
        
        llr_0_variance = 0.0;
        llr_1_variance = 0.0;
        for bit_index = 1:length(bits)
            if(llrs(bit_index) ~= -Inf && llrs(bit_index) ~= Inf)
                
                if(bits(bit_index) == 0)
                    
                    llr_0_variance = llr_0_variance + (llrs(bit_index) - llr_0_mean)^2;
                else
                    
                    llr_1_variance = llr_1_variance + (llrs(bit_index) - llr_1_mean)^2;
                end
            end
        end
        llr_0_variance = llr_0_variance/llr_0_noninfinite_count;
        llr_1_variance = llr_1_variance/llr_1_noninfinite_count;
        
        bin_width = 0.5*(3.49*sqrt(llr_0_variance)*(llr_0_noninfinite_count^(-1.0/3.0)) + 3.49*sqrt(llr_1_variance)*(llr_1_noninfinite_count^(-1.0/3.0)));
    end
    if(bin_width > 0.0)
        
        bin_offset = floor(min(llr_0_min, llr_1_min)/bin_width)-1;
        temp = max(llr_0_max, llr_1_max)/bin_width-bin_offset+1;
        bin_count = ceil(temp);
        if(bin_count == temp)
            bin_count = bin_count+1;
        end
        
    else
        
        bin_offset = -1;
        bin_count = 3;
    end
    lots_of_bins = true;
    
else
    lots_of_bins = false;
    bin_count = 4;
end

histogram = zeros(2,bin_count);

for bit_index = 1:length(bits)
    if(llrs(bit_index) == -Inf)
        histogram(bits(bit_index)+1,1) = histogram(bits(bit_index)+1,1)+1;
    elseif(llrs(bit_index) == Inf)
        histogram(bits(bit_index)+1,bin_count) = histogram(bits(bit_index)+1,bin_count)+1;
    else
        if(lots_of_bins == true)
            if(bin_width > 0.0)
                histogram(bits(bit_index)+1,floor(llrs(bit_index)/bin_width)-bin_offset+1) = histogram(bits(bit_index)+1,floor(llrs(bit_index)/bin_width)-bin_offset+1)+1;
            else
                histogram(bits(bit_index)+1,2) = histogram(bits(bit_index)+1,2)+1;
            end
        else
            histogram(bits(bit_index)+1,bits(bit_index)+2) = histogram(bits(bit_index)+1,bits(bit_index)+2)+1;
        end
    end
end

results = zeros(4, bin_count);


for bin_index=1:bin_count
%    if(histogram(1,bin_index) > 0 || histogram(2,bin_index) > 0)
        if(bin_index == 1)
            results(1,bin_index) = -inf;
%            fprintf('          -inf ');
        elseif(bin_index == bin_count)
            results(1,bin_index) = inf;
%            fprintf('           inf ');
        else
            if(lots_of_bins == true)
                if(bin_width > 0.0)
                    results(1,bin_index) = (bin_index+bin_offset-1)*bin_width+bin_width/2.0;
%                    fprintf('%14.6f ', (bin_index+bin_offset-1)*bin_width+bin_width/2.0);
                else
                    results(1,bin_index) = 0.0;
%                    fprintf('%14.6f ', 0.0);
                end
            else
                if(bin_index == 2)
                    results(1,bin_index) = -1;
%                    fprintf('           neg ');
                else
                    results(1,bin_index) = -2;
%                    fprintf('           pos ');
                end
            end
        end
        p0 = histogram(1,bin_index)/bit_0_count;
        p1 = histogram(2,bin_index)/bit_1_count;
        
        results(2,bin_index) = p0;
%        fprintf('%14.6f ', p0);
        results(3,bin_index) = p1;
%        fprintf('%14.6f ', p1);
        
        
        if(p0 == 0.0)
            results(4,bin_index) = -inf;
%            fprintf('          -inf \n');
        elseif(p1 == 0.0)
            results(4,bin_index) = inf;
%            fprintf('           inf \n');
        else
            results(4,bin_index) = log(p0/p1);
%            fprintf('%14.6f \n', log(p0/p1));
        end
%    end
end


figure;
plot(results(1,:),results(2,:),results(1,:),results(3,:))
xlabel('The values that the LLRs have');
ylabel('Histogram');
legend({'LLRs of 0-valued bits','LLRs of 1-valued bits'},'Location','northwest')
hold on

figure;
plot(results(1,:),results(4,:),'-');
axis equal
axis manual
hold on
plot([-1000,1000],[-1000,1000],'--');
xlabel('The values that the LLRs have');
ylabel('The values that the LLRs should have');

end
