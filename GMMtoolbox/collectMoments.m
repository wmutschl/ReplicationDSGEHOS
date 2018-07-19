% Original author: Martin M. Andreasen
% Modified by Willi Mutschler, September 26, 2016. Email: willi@mutschler.eu

function [moments,nameMoments] = collectMoments(Ey,Eyy,Eyyy,Eyyyy,autoEyy,autoLagsIdx,HOSThirdOrder,HOSFourthOrder)

% Selecting the moments for GMM estimation 
ny           = size(Ey,1);
numLags      = size(autoLagsIdx,2);
numMom       = ny + ny*(ny+1)/2 + ny*numLags;
moments      = zeros(numMom,1);
idxMom       = 0;
moments(idxMom+1:ny,1)= Ey;
idxMom       = idxMom + ny;
for i=1:ny
    for j=i:ny
        idxMom = idxMom + 1;
        moments(idxMom,1) = Eyy(i,j);        
    end
end

for i=1:numLags
    moments(idxMom+1:idxMom+ny,1) = diag(autoEyy(:,:,autoLagsIdx(1,i)));
    idxMom = idxMom + ny;
end


% Own higher-order product moments
if HOSThirdOrder
    [~,D3y,~] = DiagionalizationMatrix(ny,0,1,0);
    moments = [moments;D3y'*Eyyy];        
end
if HOSFourthOrder
    [~,~,D4y] = DiagionalizationMatrix(ny,0,0,1);
    moments = [moments;D4y'*Eyyyy];          
end

if nargout > 1
    % We construct a cell with general labels of the moments
    for i=1:ny
        nameMoments{i} = ['E[y_t(',num2str(i),')]'];
    end
    idxMom = ny;
    for i=1:ny
        for j=i:ny
            idxMom = idxMom + 1;
            %moments(idxMom,1) = Eyy(i,j);
            nameMoments{idxMom} = ['E[y_t(',num2str(i),')y_t(',num2str(j),')]'];
        end
    end
    for k=1:numLags
        for i=1:ny
            idxMom = idxMom + 1;
            %moments(idxMom+1:idxMom+ny,1) = diag(autoEyy(:,:,autoLagsIdx(1,i)));
            nameMoments{idxMom} = ['E[y_t-',num2str(autoLagsIdx(1,k)),'(',num2str(i),')^2]'];
        end
    end    
    if HOSThirdOrder
        for i=1:ny
            idxMom = idxMom + 1;            
            nameMoments{idxMom} = ['E[y_t(',num2str(i),')y_t(',num2str(i),')y_t(',num2str(i),')]'];            
        end
    end
    if HOSFourthOrder
        for i=1:ny
            idxMom = idxMom + 1;            
            nameMoments{idxMom} = ['E[y_t(',num2str(i),')y_t(',num2str(i),')y_t(',num2str(i),')y_t(',num2str(i),')]'];             
        end
    end    
end


end

