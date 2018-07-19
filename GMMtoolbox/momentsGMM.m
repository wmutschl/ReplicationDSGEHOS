% Original author: Martin M. Andreasen
% Modified by Willi Mutschler, September 26, 2016. Email: willi@mutschler.eu

% This function collects the following moments for GMM estimation
%  - E[y]
%  - E[y(i)*y(j)]       for i=1:ny and j=i,ny
%  - E[y(i)_t*y(i)_t-k] for i=1:ny and k=1,2,...autoLoag 
%  - E[y(i)*y(i)*y(i)]       for i=1:ny
%  - E[y(i)*y(i)*y(i)*y(i)]       for i=1:ny

function modelMoments = momentsGMM(moments,autoLagsIdx,HOSThirdOrder,HOSFourthOrder,MomComp)

if MomComp == 0
    maxLag  = max(autoLagsIdx);
    ny      = size(moments.E_y,1);
    Ey      = model.g0 + moments.E_y;
    VarY    = moments.Var_y;
    Eyy     = VarY + Ey*Ey';
    autoEyy = zeros(ny,ny,maxLag);
    for i=1:maxLag
        autoEyy(:,:,i) = moments.Cov_y(:,:,i) + Ey*Ey';
    end
    Eyyy = [];
    Eyyyy = [];
elseif MomComp == 1
    maxLag  = max(autoLagsIdx);
    ny      = size(moments.Ey,1);
    Ey      = moments.Ey;    
    Eyy     = reshape(moments.C2y0,ny,ny)+ Ey*Ey';
    autoEyy = zeros(ny,ny,maxLag);
    for i=1:maxLag
        autoEyy(:,:,i) = reshape(moments.C2yt(:,i),ny,ny)+ Ey*Ey';        
    end    
    % Note that for mean zero cumulants C2y0, C3y0 and C4y0 we have the following relationships to nondemeaned productmoments Eyy, Eyyy, Eyyyy:
    %     C2y0 = Eyy - kron(Ey,Ey);
    %     C3y0 = Eyyy  - (eye(ny^3) + commutation(ny,ny^2) + commutation(ny^2,ny))*kron(Ey,Eyy) + (3-1)*kron(Ey,kron(Ey,Ey));
    %     C4y0 = Eyyyy - (eye(ny^4) + commutation(ny,ny^3) + commutation(ny^2,ny^2) +commutation(ny^3,ny))*kron(Ey,Eyyy)...
    %                  + (eye(ny^4) + kron(eye(ny^2),commutation(ny,ny)) + kron(commutation(ny,ny),eye(ny^2)) + kron(commutation(ny,ny),commutation(ny,ny)))*kron(Ey,kron(Eyy,Ey))...
    %                  + (eye(ny^4) + commutation(ny^2,ny^2))*kron(Ey,kron(Ey,Eyy))...
    %                  - (4-1)*kron(Ey,kron(Ey,kron(Ey,Ey)))...
    %                  - (eye(ny^4)+PermutMat(ny)'+PermutMat(ny))*kron(C2y0,C2y0);
    % Note that for mean zero cumulants C2y0, C3y0 and C4y0 we have the following relationships to demeaned productmoments Eyy, Eyyy, Eyyyy:
    %     C2y0 = Eyy;
    %     C3y0 = Eyyy;
    %     C4y0 = Eyyyy - (eye(ny^4)+PermutMat(ny)'+PermutMat(ny))*kron(C2y0,C2y0);
    
    if HOSThirdOrder
        Eyyy = moments.C3y0 + (eye(ny^3) + commutation(ny,ny^2) + commutation(ny^2,ny))*kron(Ey(:),Eyy(:)) - (3-1)*kron(Ey,kron(Ey,Ey));        
    else
        Eyyy = [];        
    end
    if HOSFourthOrder
        Eyyy = moments.C3y0 + (eye(ny^3) + commutation(ny,ny^2) + commutation(ny^2,ny))*kron(Ey(:),Eyy(:)) - (3-1)*kron(Ey,kron(Ey,Ey));
        Eyyyy = moments.C4y0 + (eye(ny^4) + commutation(ny,ny^3) + commutation(ny^2,ny^2) +commutation(ny^3,ny))*kron(Ey,Eyyy(:))...
                                 - (eye(ny^4) + kron(eye(ny^2),commutation(ny,ny)) + kron(commutation(ny,ny),eye(ny^2)) + kron(commutation(ny,ny),commutation(ny,ny)))*kron(Ey,kron(Eyy(:),Ey))...
                                 - (eye(ny^4) + commutation(ny^2,ny^2))*kron(Ey,kron(Ey,Eyy(:)))...
                                 + (4-1)*kron(Ey,kron(Ey,kron(Ey,Ey)))...
                                 + (eye(ny^4)+PermutMat(ny)'+PermutMat(ny))*kron(moments.C2y0,moments.C2y0);        
    else
        Eyyyy = [];        
    end        
end
% We collect the moments we need for estimation
modelMoments = collectMoments(Ey,Eyy,Eyyy,Eyyyy,autoEyy,autoLagsIdx,HOSThirdOrder,HOSFourthOrder);
