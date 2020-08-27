% BCJR algorithm for a half-rate systematic recursive convolutional code
% having 1 memory element, a generator polynomial of [1,0] and a feedback
% polynomial of [1,1]. For more information, see Section 1.3.2.2 of Rob's
% thesis (http://eprints.ecs.soton.ac.uk/14980) or the BCJR paper 
% (http://ieeexplore.ieee.org/xpls/abs_all.jsp?arnumber=1055186).
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



function [aposteriori_uncoded_llrs, aposteriori_encoded1_llrs, aposteriori_encoded2_llrs] = bcjr_decoder(apriori_uncoded_llrs, apriori_encoded1_llrs, apriori_encoded2_llrs)

    if(length(apriori_uncoded_llrs) ~= length(apriori_encoded1_llrs) || length(apriori_encoded1_llrs) ~= length(apriori_encoded2_llrs))
        error('LLR sequences must have the same length');
    end


    % All calculations are performed in the logarithmic domain in order to
    % avoid numerical issues. These occur in the normal domain, because some of 
    % the confidences can get smaller than the smallest number the computer can
    % store. See Section 1.3.2.4 of Rob's thesis for more information on this.
    %
    % A multiplication of two confidences is achieved using the addition of the
    % corresponding log-confidences. If A = log(a) and B = log(b), then
    % log(a*b) = A+B (Equation 1.17 in Rob's thesis).
    %
    % An addition of two confidences is achieved using the Jacobian logarithm
    % of the corresponding log-confidences. The Jacobian logarithm is defined
    % in the jac.m file. If A = log(a) and B = log(b), then 
    % log(a+b) = max(A,B) + log(1+exp(-abs(A-B))) (Equation 1.19 in Rob's
    % thesis).

    % Matrix to describe the trellis
    % Each row describes one transition in the trellis
    % Each state is allocated an index 1,2,3,... Note that this list starts 
    % from 1 rather than 0.
    %              FromState,   ToState,     UncodedBit,  Encoded1Bit, Encoded2Bit
    transitions = [1,           1,           0,           0,           0;
                   1,           2,           1,           1,           1;
                   2,           1,           1,           1,           0;
                   2,           2,           0,           0,           1];

    % Find the largest state index in the transitions matrix           
    % In this example, we have two states since the code has one memory element
    state_count = max(max(transitions(:,1)),max(transitions(:,2)));

    % Calculate the a priori transition log-confidences by adding the
    % log-confidences associated with each corresponding bit value. This is
    % similar to Equation 1.12 in Rob's thesis or Equation 9 in the BCJR paper.
    gammas=zeros(size(transitions,1),length(apriori_uncoded_llrs));
    for bit_index = 1:length(apriori_uncoded_llrs)
       for transition_index = 1:size(transitions,1)
          if transitions(transition_index, 3)==0
              gammas(transition_index, bit_index) = gammas(transition_index, bit_index) + apriori_uncoded_llrs(bit_index)/2; 
              % Dividing the LLR by 2 gives a value that is log-proportional to
              % the actual log-confidence it instills. Log-proportional is fine
              % for us though, because we're after the log-ratio of
              % confidences. Don't worry too much about this confusing
              % feature.
          else
              gammas(transition_index, bit_index) = gammas(transition_index, bit_index) - apriori_uncoded_llrs(bit_index)/2;
          end

          if transitions(transition_index, 4)==0
              gammas(transition_index, bit_index) = gammas(transition_index, bit_index) + apriori_encoded1_llrs(bit_index)/2;
          else
              gammas(transition_index, bit_index) = gammas(transition_index, bit_index) - apriori_encoded1_llrs(bit_index)/2;
          end

          if transitions(transition_index, 5)==0
              gammas(transition_index, bit_index) = gammas(transition_index, bit_index) + apriori_encoded2_llrs(bit_index)/2;
          else
              gammas(transition_index, bit_index) = gammas(transition_index, bit_index) - apriori_encoded2_llrs(bit_index)/2;
          end
       end
    end

    % Forward recursion to calculate state log-confidences. This is similar to
    % Equation 1.13 in Rob's thesis or Equations 5 and 6 in the BCJR paper.
    alphas=zeros(state_count,length(apriori_uncoded_llrs));
    alphas=alphas-inf;
    alphas(1,1)=0; % We know that this is the first state
    for state_index = 2:state_count
        alphas(state_index,1)=-inf; % We know that this is *not* the first state (a log-confidence of minus infinity is equivalent to a confidence of 0)
    end
    for bit_index = 2:length(apriori_uncoded_llrs)
       for transition_index = 1:size(transitions,1)
           alphas(transitions(transition_index,2),bit_index) = jac(alphas(transitions(transition_index,2),bit_index),alphas(transitions(transition_index,1),bit_index-1) + gammas(transition_index, bit_index-1));   
       end
    end

    % Backwards recursion to calculate state log-confidences. This is similar
    % to Equation 1.14 in Rob's thesis or Equations 7 and 8 in the BCJR paper.
    betas=zeros(state_count,length(apriori_uncoded_llrs));
    betas=betas-inf;
    for state_index = 1:state_count
        betas(state_index,length(apriori_uncoded_llrs))=0; % The final state could be any one of these
    end
    for bit_index = length(apriori_uncoded_llrs)-1:-1:1
       for transition_index = 1:size(transitions,1)
           betas(transitions(transition_index,1),bit_index) = jac(betas(transitions(transition_index,1),bit_index),betas(transitions(transition_index,2),bit_index+1) + gammas(transition_index, bit_index+1));   
       end
    end

    % Calculate a posteriori transition log-confidences. This is similar to
    % Equation 1.15 in Rob's thesis or Equation 4 in the BCJR paper.
    deltas=zeros(size(transitions,1),length(apriori_uncoded_llrs));
    for bit_index = 1:length(apriori_uncoded_llrs)
       for transition_index = 1:size(transitions,1)
           deltas(transition_index, bit_index) = alphas(transitions(transition_index,1),bit_index) + gammas(transition_index, bit_index) + betas(transitions(transition_index,2),bit_index);
       end
    end

    % Calculate the a posteriori LLRs. This is similar to Equation 1.16 in
    % Rob's thesis.
    aposteriori_uncoded_llrs = zeros(1,length(apriori_uncoded_llrs));
    for bit_index = 1:length(apriori_uncoded_llrs)    
       prob0=-inf;
       prob1=-inf;
       for transition_index = 1:size(transitions,1)
           if transitions(transition_index,3)==0
               prob0 = jac(prob0, deltas(transition_index,bit_index));
           else
               prob1 = jac(prob1, deltas(transition_index,bit_index));
           end      
       end
       aposteriori_uncoded_llrs(bit_index) = prob0-prob1;
    end

    aposteriori_encoded1_llrs = zeros(1,length(apriori_uncoded_llrs));
    for bit_index = 1:length(apriori_uncoded_llrs)    
       prob0=-inf;
       prob1=-inf;
       for transition_index = 1:size(transitions,1)
           if transitions(transition_index,4)==0
               prob0 = jac(prob0, deltas(transition_index,bit_index));
           else
               prob1 = jac(prob1, deltas(transition_index,bit_index));
           end      
       end   
       aposteriori_encoded1_llrs(bit_index) = prob0-prob1;
    end

    aposteriori_encoded2_llrs = zeros(1,length(apriori_uncoded_llrs));
    for bit_index = 1:length(apriori_uncoded_llrs)    
       prob0=-inf;
       prob1=-inf;
       for transition_index = 1:size(transitions,1)
           if transitions(transition_index,5)==0
               prob0 = jac(prob0, deltas(transition_index,bit_index));
           else
               prob1 = jac(prob1, deltas(transition_index,bit_index));
           end      
       end   
       aposteriori_encoded2_llrs(bit_index) = prob0-prob1;
    end

end