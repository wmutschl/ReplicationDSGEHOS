% Original author: Martin M. Andreasen
% Modified by Willi Mutschler, September 26, 2016. Email: willi@mutschler.eu

% This function computes the optimal weigthing matrix by a Bartlett kernel
% with max lag q_lag
function Wopt = getOptimalWeighting(moments,data,autoLagsIdx,HOSThirdOrder,HOSFourthOrder,qLag)
% We compute the h-function for all observations
T = size(data,1);

hFunc = get_hFunc(data,autoLagsIdx,moments,HOSThirdOrder,HOSFourthOrder);

% The required correlation matrices
numMom = size(hFunc,2); %Number of moments
GAMA_array = zeros(numMom,numMom,qLag);
GAMA0 = CorrMatrix(hFunc,T,numMom,0);
if qLag > 0
    for i=1:qLag
        GAMA_array(:,:,i) = CorrMatrix(hFunc,T,numMom,i);
    end
end

% The estimate of S
S = GAMA0;
if qLag > 0
    for i=1:qLag
        S = S + (1-i/(qLag+1))*(GAMA_array(:,:,i) + GAMA_array(:,:,i)');
    end
end
Wopt = S\eye(size(S,1));

end

% The correlation matrix
function GAMAcorr = CorrMatrix(hFunc,T,numMom,v)
GAMAcorr = zeros(numMom,numMom);
for t=1+v:T
    GAMAcorr = GAMAcorr + hFunc(t-v,:)'*hFunc(t,:);
end
GAMAcorr = GAMAcorr/T;    
end
