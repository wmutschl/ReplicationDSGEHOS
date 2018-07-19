% Original author: Martin M. Andreasen
% Modified by Willi Mutschler, September 26, 2016. Email: willi@mutschler.eu

function [dataMoments,nameMoments] = momentsGMMData(data,autoLagsIdx,HOSThirdOrder,HOSFourthOrder)
% This function computes the following empirical moments from data
%  - E[y]
%  - E[y(i)*y(j)]       for i=1:ny and j=i,ny
%  - E[y(i)_t*y(i)_t-k] for i=1:ny and k=1,2,...autoLoag 
%  - E[y(i)*y(i)*y(i)]       for i=1:ny
%  - E[y(i)*y(i)*y(i)*y(i)]  for i=1:ny

% Number of observations (T) and number of observables (ny)
[T,ny] = size(data);

% The means
Ey       = sum(data,1)'/T;
Eyy      = data'*data/T;
maxLag   = max(autoLagsIdx);
autoEyy  = zeros(ny,ny,maxLag);
for i=1:maxLag
    autoEyy(:,:,i) = data(1+i:T,:)'*data(1:T-i,:)/(T-i);
end

if HOSThirdOrder
    Eyyy  = zeros(ny,ny^2,1);
    for i = 1:ny
        dati = data(:,i);            
        Eyyy(1:ny, ((i-1)*ny+1):(i*ny))= 1/T*(data'*diag(dati)*data);
    end
    Eyyy = vec(Eyyy);    
else
    Eyyy = [];
end

if HOSFourthOrder
    Eyyyy = zeros(ny,ny^3);                
    id = 1;
    for i = 1:ny
        for j = 1:ny        
            dati = data(:,i);
            datj = data(:,j);                
            Eyyyy(1:ny,((id-1)*ny+1):(id*ny)) = 1/T*(data'*diag(dati)*diag(datj)*data); 
            id = id+1;
        end
    end
    Eyyyy = vec(Eyyyy);    
else
    Eyyyy = [];
end

% We collect the moments we need for estimation
[dataMoments,nameMoments] = collectMoments(Ey,Eyy,Eyyy,Eyyyy,autoEyy,autoLagsIdx,HOSThirdOrder,HOSFourthOrder);

%% Alternative computations
% Note that for mean zero cumulants C2y0, C3y0 and C4y0 we have the following relationships to demeaned productmoments Eyy, Eyyy, Eyyyy:
%     C2y0 = Eyy;
%     C3y0 = Eyyy;
%     C4y0 = Eyyyy - (eye(ny^4)+PermutMat(ny)'+PermutMat(ny))*kron(C2y0,C2y0);
% Note that for mean zero cumulants C2y0, C3y0 and C4y0 we have the following relationships to nondemeaned productmoments Eyy, Eyyy, Eyyyy:
%     C2y0 = Eyy - kron(Ey,Ey);
%     C3y0 = Eyyy  - (eye(ny^3) + commutation(ny,ny^2) + commutation(ny^2,ny))*kron(Ey,Eyy) + (3-1)*kron(Ey,kron(Ey,Ey));
%     C4y0 = Eyyyy - (eye(ny^4) + commutation(ny,ny^3) + commutation(ny^2,ny^2) +commutation(ny^3,ny))*kron(Ey,Eyyy)...
%                  + (eye(ny^4) + kron(eye(ny^2),commutation(ny,ny)) + kron(commutation(ny,ny),eye(ny^2)) + kron(commutation(ny,ny),commutation(ny,ny)))*kron(Ey,kron(Eyy,Ey))...
%                  + (eye(ny^4) + commutation(ny^2,ny^2))*kron(Ey,kron(Ey,Eyy))...
%                  - (4-1)*kron(Ey,kron(Ey,kron(Ey,Ey)))...
%                  - (eye(ny^4)+PermutMat(ny)'+PermutMat(ny))*kron(C2y0,C2y0);

% Alternative 1 for mean-zero second, third and fourth cumulants
%     dataHOS = bsxfun(@minus,data,mean(data));
%     Cum2y = zeros(ny*(ny+1)/(1*2),1); id2=1;    
%     Cum3y = zeros(ny*(ny+1)*(ny+2)/(1*2*3),1); id3=1;    
%     Cum4y = zeros(ny*(ny+1)*(ny+2)*(ny+3)/(1*2*3*4),1); id4=1;            
%     for i=1:ny
%         ws = dataHOS(:,i);
%         for j = i:ny
%             xs = dataHOS(:,j);
%             for k = j:ny      
%                 ys = dataHOS(:,k);
%                 for l = k:ny
%                     zs = dataHOS(:,l);
%                     Cum4y(id4) = cum4x(data(:,i),data(:,j),data(:,k),data(:,l));                                            
%                     R_wy = ws' * ys /T;                    
%                     M_wz = ws' * zs /T;
%                     R_zy = zs' * ys /T;                                        
%                     R_wx = ws' * xs /T;
%                     R_zx = zs' * xs /T;
%                     M_yx = ys' * xs /T;
%                     Cum4y(id4) = (transpose(ws .* ys .* zs) * xs   - R_zy * R_wx *T - R_wy * R_zx*T - M_wz'* M_yx *T )/T;
%                     id4 = id4+1;
%                 end
%                 Cum3y(id3) = cum3x(data(:,i),data(:,j),data(:,k));
%                 Cum3y(id3) = transpose(ws .* ys) * xs /T;                
%                 id3 = id3+1;
%             end
%             Cum2y(id2) = cum2x(dataHOS(:,i),dataHOS(:,j));
%             Cum2y(id2) = ws' * xs /T;
%             id2 = id2+1;
%         end
%     end
%
% Alternative 2 for mean-zero second, third and fourth cumulants
%     dataHOS = bsxfun(@minus,data,mean(data));
%     Cum2y = zeros(ny^2,1);    
%     Cum2y = Eyy(:) - vec(Ey*Ey');     
%     Cum3y = zeros(ny^3,1);    
%     Cum4y = zeros(ny^4,1); Pny = PermutMat(ny); auxCum4= speye(ny^4)+Pny'+Pny;
%     aux = auxCum4*kron(Cum2y,Cum2y)/T;
%     for t=1:T
%         dat_kron_dat = kron(dataHOS(t,:)',dataHOS(t,:)');
%         Cum2y = Cum2y + dat_kron_dat/T;
%         Cum3y = Cum3y + kron(dataHOS(t,:)',dat_kron_dat)/T;        
%         Cum4y = Cum4y + kron(dat_kron_dat,dat_kron_dat)/T - aux ;
%     end    

% % Test for lags
%     t1=1; t2=1;t3=1;
%     dataHOS = bsxfun(@minus,data,mean(data));
%     Cum2y = zeros(ny*(ny+1)/(1*2),1); id2=1;    
%     Cum3y = zeros(ny*(ny+1)*(ny+2)/(1*2*3),1); id3=1;    
%     Cum4y = zeros(ny*(ny+1)*(ny+2)*(ny+3)/(1*2*3*4),1); id4=1;    
%     for i=1:ny
%         ws = dataHOS(1:(T-t1-t2-t3),i);
%         for j = i:ny
%             xs = dataHOS((1+t1):(T-t2-t3),j);
%             for k = j:ny      
%                 ys = dataHOS((1+t2):(T-t1-t3),k);
%                 for l = k:ny
%                     zs = dataHOS((1+t3):(T-t1-t2),l);
%                     Cum4y(id4) = cum4x(data(:,i),data(:,j),data(:,k),data(:,l));
%                     R_wy = ws' * ys /T;                    
%                     M_wz = ws' * zs /T;
%                     R_zy = zs' * ys /T;                                        
%                     R_wx = ws' * xs /T;
%                     R_zx = zs' * xs /T;
%                     M_yx = ys' * xs /T;
%                     Cum4y(id4) = (transpose(ws .* ys .* zs) * xs   - R_zy * R_wx *T - R_wy * R_zx*T - M_wz'* M_yx *T )/T;
%                     id4 = id4+1;
%                 end
%                 Cum3y(id3) = cum3x(data(:,i),data(:,j),data(:,k));
%                 Cum3y(id3) = transpose(ws .* ys) * xs /T;                
%                 id3 = id3+1;
%             end
%             Cum2y(id2) = cum2x(dataHOS(:,i),dataHOS(:,j));
%             Cum2y(id2) = ws' * xs /T;
%             id2 = id2+1;
%         end
%     end
%     DP = duplication(ny); TP = triplication(ny); QP = quadruplication(ny);
%     Pny = PermutMat(ny); auxCum4= speye(ny^4)+Pny'+Pny;    
%     % mean zero product moments
%     PM2y = DP*Cum2y;
%     PM3y = TP*Cum3y;
%     PM4y = QP*Cum4y + auxCum4*kron(DP*Cum2y,DP*Cum2y);
end

