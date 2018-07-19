% Original author: Martin M. Andreasen
% Modified by Willi Mutschler, September 26, 2016. Email: willi@mutschler.eu
function hFunc = get_hFunc(data,autoLagsIdx,modelMoments,HOSThirdOrder,HOSFourthOrder)
% Thus function reports the h-function for GMM estimation, i.e. 
% h = m(data) - m(theta). 

% Number of observations (T) and number of observables (ny)
[T,ny] = size(data);

% Number of moments
numLags      = size(autoLagsIdx,2);
numMom       = ny + ny*(ny+1)/2 + ny*numLags;

% Preallocation of memory
hFunc = zeros(T,numMom);

% For Ey
idxMom = 0;
hFunc(:,idxMom+1:idxMom+ny) = data;
idxMom = idxMom + ny;

% For Eyy
for i=1:ny
    for j=i:ny
        idxMom = idxMom + 1;
        hFunc(:,idxMom) = data(:,i).*data(:,j);        
    end
end

% For autoEyy
maxLag   = max(autoLagsIdx);
autoEyy  = zeros(T,ny,maxLag);
for i=1:maxLag
    for j=1:ny
        autoEyy(1:i,j,i)   = mean(data(1+i:T,j).*data(1:T-i,j));
        autoEyy(1+i:T,j,i) = data(1+i:T,j).*data(1:T-i,j);
    end
end

numLags      = size(autoLagsIdx,2);
for i=1:numLags
    for j=1:ny
        idxMom = idxMom + 1;
        hFunc(:,idxMom) = autoEyy(:,j,autoLagsIdx(1,i));
    end
end

%% For HOS

% Own product moments    
if HOSThirdOrder
    EYYY  = zeros(T,ny);
    for i=1:ny
        dati = data(:,i);
        EYYY(:,i) = dati.*dati.*dati;
    end
    hFunc = [hFunc EYYY];
end
        
if HOSFourthOrder
    EYYYY = zeros(T,ny);                        
    for i=1:ny
        dati = data(:,i);            
        EYYYY(:,i) = dati.*dati.*dati.*dati;                
    end
    hFunc = [hFunc EYYYY];
end

hFunc = hFunc - repmat(modelMoments',T,1);

end

