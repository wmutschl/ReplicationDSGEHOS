% By Willi Mutschler, September 26, 2016. Email: willi@mutschler.eu
function HOS = HOSrun(setup,HOSlowerorder,C2z0_old,C3z0_old,C4z0_old,auxOut)
global M_

% Get results from lower order approximation if possible
if isempty(HOSlowerorder) == 0
    if setup.orderApp == 2
        E_XF1 = HOSlowerorder.E_XF1;
        E_XF2min = HOSlowerorder.E_XF2min;
        E_XF3min = HOSlowerorder.E_XF3min;
        E_XF4min = HOSlowerorder.E_XF4min;        
    elseif setup.orderApp == 3
        E_XF1 = HOSlowerorder.E_XF1;
        E_XF2min = HOSlowerorder.E_XF2min; 
        E_XF3min = HOSlowerorder.E_XF3min; 
        E_XF4min = HOSlowerorder.E_XF4min; 
        E_XF5min = HOSlowerorder.E_XF5min; 
        E_XF6min = HOSlowerorder.E_XF6min;
        E_XF7min = HOSlowerorder.E_XF7min;
        E_XF8min = HOSlowerorder.E_XF8min;
        E_XS1 = HOSlowerorder.E_XS1;
        E_XS2min = HOSlowerorder.E_XS2min;
        E_XS3min = HOSlowerorder.E_XS3min;
        E_XS4min = HOSlowerorder.E_XS4min;
        E_XF1_XS1 = HOSlowerorder.E_XF1_XS1; 
        E_XF2_XS1min = HOSlowerorder.E_XF2_XS1min;
        E_XF1_XS2min = HOSlowerorder.E_XF1_XS2min ; 
        E_XF3_XS1min = HOSlowerorder.E_XF3_XS1min;
        E_XF2_XS2min = HOSlowerorder.E_XF2_XS2min; 
        E_XF1_XS3min = HOSlowerorder.E_XF1_XS3min;
        E_XF4_XS1min = HOSlowerorder.E_XF4_XS1min;
        E_XF3_XS2min = HOSlowerorder.E_XF3_XS2min;
        E_XF2_XS3min = HOSlowerorder.E_XF2_XS3min;
        E_XF5_XS1min = HOSlowerorder.E_XF5_XS1min;
        E_XF4_XS2min = HOSlowerorder.E_XF4_XS2min;
        E_XF6_XS1min = HOSlowerorder.E_XF6_XS1min;
    end
end
%% Numerical computation of cumulants and moments

[A,B,C,D,c,d,ybar,Fxi,Fz,nx,ny,nu,nximin,id_xf,id_xs,id_xf_xf] = GetStateSpaceRep(setup);
if setup.orderApp == 1
    nameXiFile2 = sprintf('Xi_%s_prodmom2_orderApp%i_nu%i',setup.distrib,1,nu);
    nameXiFile3 = sprintf('Xi_%s_prodmom3_orderApp%i_nu%i',setup.distrib,1,nu);
    nameXiFile4 = sprintf('Xi_%s_prodmom4_orderApp%i_nu%i',setup.distrib,1,nu);
else    
    nameXiFile2 = sprintf('Xi_%s_prodmom2_orderApp%i_nx%i_nu%i',setup.distrib,setup.orderApp,nx,nu);
    nameXiFile3 = sprintf('Xi_%s_prodmom3_orderApp%i_nx%i_nu%i',setup.distrib,setup.orderApp,nx,nu);
    nameXiFile4 = sprintf('Xi_%s_prodmom4_orderApp%i_nx%i_nu%i',setup.distrib,setup.orderApp,nx,nu);    
end

% Auxiliary parameters for distribution
[~,DPuinv] = duplication(nu);
switch setup.distrib
    case 'Gaussian'
        argdistrib=DPuinv*M_.Sigma_e(:);
    case 'Student_t'
        for i = 1:M_.param_nbr
            if strcmp(deblank(M_.param_names(i,:)),'df_studt')
                df_studt = M_.params(i);                
            end
        end
        argdistrib =[df_studt;(df_studt-2)/df_studt*DPuinv*M_.Sigma_e(:)];
end

%% First-order moments, ie expectation of variables
IminA = speye(size(A,1))-A;
Ez   = IminA\c;
Ey   = ybar + C*Ez+d; % recall y = yss + C*z + d
    

%% Compute Zero-Lag Cumulants of innovations, states and controls
nz = size(A,1);

BFxi = B*Fxi;
DFxi = D*Fxi;
if setup.HOSThirdOrder || setup.HOSFourthOrder
    if strcmp(setup.sylv_method,'disclyap') == 0
        AkronA = kron(A,A);
    end
    CkronC = kron(C,C);
    BFxikronBFxi= kron(BFxi,BFxi);    
    DFxikronDFxi= kron(DFxi,DFxi);
end

% Auxiliary parameters to evaluate script file for GAMMA2xi
if setup.orderApp == 1
    arg2 = argdistrib;
elseif setup.orderApp == 2
    arg2 = [argdistrib;E_XF1;E_XF2min];    
elseif setup.orderApp == 3    
    arg2 = [argdistrib; E_XF1; E_XF2min; E_XF3min; E_XF4min; E_XS1;E_XS2min; E_XF1_XS1; E_XF2_XS1min]; % These are the product moments we need to compute GAMMA2XI        
end

% Second-order unique product moments of innovation vector
if exist(nameXiFile2,'file') == 0
    fprintf('Sorry, we need to first create 2nd order HOS script files for %d-order approximation.\n This could take a while, but note that this needs to be done only once for any model with the same number of shocks and states.\n',setup.orderApp);        
    createHOSScriptFiles(nu,nx,setup.orderApp,setup.distrib,2);
end
evalFile2 = str2func(nameXiFile2);
GAMMA2ximin = evalFile2(arg2);
try 
    load(sprintf('DP_%i',nximin),'DP'); 
catch
    [DP,DPinv] = duplication(nximin);
    save(sprintf('DP_%i',nximin),'DP','DPinv');
    fprintf('You may copy DP_%i.mat into the HOSScriptFiles/Other folder\n',nximin);
end
GAMMA2XI = reshape(DP*GAMMA2ximin,nximin,nximin); %Note GAMMA2 is a matrix, not vector
CC = BFxi*GAMMA2XI*transpose(BFxi);
% Calculate zero-lag cumulant of states and controls for second-order
if setup.sparseflag == 0
    CC = full(CC);
end
if sum(isnan(C2z0_old)) > 0; C2z0_old = []; end
if isempty(C2z0_old) == 0; C2z0_old = reshape(C2z0_old,nz,nz); end
if strcmp(setup.sylv_method,'disclyap')
    C2z0 = disclyap(A,CC,C2z0_old,setup.sylv_tol,setup.sylv_max_iter,setup.sylv_sparseflag,1,2);      
else
    C2z0 = solve_gensylv(A,transpose(A),CC,C2z0_old,setup.sylv_tol,setup.sylv_max_iter,setup.sylv_sparseflag,0,setup.sylv_method);      
end
C2y0 = C*C2z0*transpose(C) + DFxi*GAMMA2XI*transpose(DFxi);


if setup.HOSThirdOrder || setup.HOSFourthOrder
    % Third-order unique product moments of innovation vector
    % Auxiliary parameters to evaluate script file for GAMMA3xi
    if setup.orderApp == 1
        arg3 = argdistrib;
    elseif setup.orderApp == 2
        arg3 = [argdistrib;E_XF1;E_XF2min;E_XF3min];        
    elseif setup.orderApp == 3
        arg3 = [argdistrib;E_XF1;E_XF2min; E_XF3min; E_XF4min; E_XF5min; E_XF6min;...
                E_XS1; E_XS2min; E_XS3min;...
                E_XF1_XS1; 
                E_XF2_XS1min; E_XF1_XS2min; 
                E_XF3_XS1min; E_XF2_XS2min; 
                E_XF4_XS1min]; % These are the product moments we need to compute GAMMA3XI        
    end
    % Calculate zero-lag cumulant of states and controls for third-order
    if exist(nameXiFile3,'file') == 0
        fprintf('Sorry, we need to first create 3rd order HOS script files for %d-order approximation.\n This could take a while, but note that this needs to be done only once for any model with the same number of shocks and states.\n',setup.orderApp);        
        createHOSScriptFiles(nu,nx,setup.orderApp,setup.distrib,3);
    end
    evalFile3 = str2func(nameXiFile3);
    GAMMA3ximin = evalFile3(arg3);
    try 
        load(sprintf('TP_%i',nximin),'TP'); 
    catch
        [TP,TPinv] = triplication(nximin,1);
        save(sprintf('TP_%i',nximin),'TP','TPinv'); 
        fprintf('You may copy TP_%i.mat into the HOSScriptFiles/Other folder\n',nximin);
    end
    GAMMA3XI = reshape(TP*GAMMA3ximin,nximin^2,nximin); %Note GAMMA3 is a matrix, not vector
    CC = BFxikronBFxi*GAMMA3XI*transpose(BFxi);
    if setup.sparseflag == 0
        CC = full(CC);
    end
    if sum(isnan(C3z0_old)) > 0; C3z0_old = []; end
    if isempty(C3z0_old) == 0; C3z0_old = reshape(C3z0_old,nz^2,nz); end
    if strcmp(setup.sylv_method,'disclyap')
        C3z0 = disclyap(A,CC,C3z0_old,setup.sylv_tol,setup.sylv_max_iter,setup.sylv_sparseflag,0,3);        
    else
        C3z0 = solve_gensylv(AkronA,transpose(A),CC,C3z0_old,setup.sylv_tol,setup.sylv_max_iter,setup.sylv_sparseflag,0,setup.sylv_method);          
    end
    C3y0 = CkronC*C3z0*transpose(C) + DFxikronDFxi*GAMMA3XI*transpose(DFxi);
else
    C3z0=[];
    C3y0=[];
    GAMMA3XI=[];
end

if setup.HOSFourthOrder
    % Fourth-order unique product moments of innovation vector
    % Auxiliary parameters to evaluate script file for GAMMA4xi
    if setup.orderApp == 1
        arg4 = argdistrib;
    elseif setup.orderApp == 2
        arg4 = [argdistrib; E_XF1;E_XF2min;E_XF3min;E_XF4min];        
    elseif setup.orderApp == 3
        arg4 = [argdistrib;
                E_XF1; E_XF2min; E_XF3min; E_XF4min; E_XF5min; E_XF6min; E_XF7min; E_XF8min;
                E_XS1; E_XS2min; E_XS3min; E_XS4min;
                E_XF1_XS1; 
                E_XF2_XS1min; E_XF1_XS2min;
                E_XF3_XS1min; E_XF2_XS2min; E_XF1_XS3min;
                E_XF4_XS1min; E_XF3_XS2min; E_XF2_XS3min;
                E_XF5_XS1min; E_XF4_XS2min;
                E_XF6_XS1min]; % These are the product moments we need to compute GAMMA4XI                        
    end
    % Calculate zero-lag cumulant for fourth-order
    if exist(nameXiFile4,'file') == 0
        fprintf('Sorry, we need to first create 4th order HOS script files for %d-order approximation.\n This could take a while, but note that this needs to be done only once for any model with the same number of shocks and states.\n',setup.orderApp);        
        createHOSScriptFiles(nu,nx,setup.orderApp,setup.distrib,4);
    end
    evalFile4 = str2func(nameXiFile4);
    try 
        load(sprintf('QP_%i',nximin),'QP','QPinv'); 
    catch
        [QP,QPinv] = quadruplication(nximin,1);
        save(sprintf('QP_%i',nximin),'QP','QPinv'); 
        fprintf('You may copy QP_%i.mat into the HOSScriptFiles/Other folder\n',nximin);
    end
    try
        load(sprintf('auxCUM4xi_%i',nximin),'auxCUM4xi');
    catch
        PermMat = PermutMat(nximin); % Permutation matrix for fourth-order cumulant
        auxCUM4xi = QPinv*(speye(nximin^4)+transpose(PermMat)+PermMat)*kron(DP,DP);
        save(sprintf('auxCUM4xi_%i',nximin),'auxCUM4xi');
        fprintf('You may copy auxCUM4xi_%i.mat into the HOSScriptFiles/Other folder\n',nximin);
    end
    GAMMA4ximin = evalFile4(arg4) - auxCUM4xi*kron(GAMMA2ximin,GAMMA2ximin);                
    GAMMA4XI = reshape(QP*GAMMA4ximin,nximin^2,nximin^2);%Note GAMMA4 is a matrix, not vector    
    CC = BFxikronBFxi*GAMMA4XI*transpose(BFxikronBFxi);
    if setup.sparseflag == 0
        CC = full(CC);
    end
    if sum(isnan(C4z0_old)) > 0; C4z0_old = []; end
    if isempty(C4z0_old) == 0; C4z0_old = reshape(C4z0_old,nz^2,nz^2); end    
    if strcmp(setup.sylv_method,'disclyap')
        C4z0 = disclyap(A,CC,C4z0_old,setup.sylv_tol,setup.sylv_max_iter,setup.sylv_sparseflag,1,4);        
    else
        C4z0 = solve_gensylv(AkronA,transpose(AkronA),CC,C4z0_old,setup.sylv_tol,setup.sylv_max_iter,setup.sylv_sparseflag,1,setup.sylv_method);        
    end
    C4y0 = CkronC*C4z0*transpose(CkronC) + DFxikronDFxi*GAMMA4XI*transpose(DFxikronDFxi);
else
    C4z0=[];
    C4y0=[];
    GAMMA4XI=[];
end
    
%% Get autocumulants for lag order T
max2lags = max(setup.autoLagsIdx);
C2yt = zeros(ny^2,max2lags); 
At1 = A;
for t1 = 1:max2lags          
    C2yt(:,t1) = vec(C*At1*C2z0*transpose(C));
    At1 = At1*A;    
end


%% Store everything into structure
HOS.Ez = Ez;
HOS.Ey = Ey;

HOS.C2z0 = C2z0(:);
HOS.C2y0 = C2y0(:);
HOS.C2yt = C2yt;
HOS.GAMMA2XI = GAMMA2XI(:);

HOS.C3z0 = C3z0(:);
HOS.C3y0 = C3y0(:);
HOS.GAMMA3XI = GAMMA3XI(:);

HOS.C4z0 = C4z0(:); 
HOS.C4y0 = C4y0(:);
HOS.GAMMA4XI = GAMMA4XI(:);

% Additional results
if auxOut
    if setup.orderApp == 1
        E_XF1 = Ez;
        E_XF2 = C2z0(:); %since E(xf) = 0
        HOS.E_XF1 = E_XF1;
        try 
            load(sprintf('DP_%i',nx),'DPinv'); 
        catch
            [DP,DPinv] = duplication(nx); 
            save(sprintf('DP_%i',nx),'DP','DPinv');
            fprintf('You may copy DP_%i.mat into the HOSScriptFiles/Other folder\n',nx);
        end
        HOS.E_XF2min = DPinv*E_XF2;
        if setup.HOSThirdOrder || setup.HOSFourthOrder
            E_XF3 = C3z0(:); %since E(xf) = 0
            try 
                load(sprintf('TP_%i',nx),'TPinv'); 
            catch
                [TP,TPinv] = triplication(nx,1);
                save(sprintf('TP_%i',nx),'TP','TPinv');
                fprintf('You may copy TP_%i.mat into the HOSScriptFiles/Other folder\n',nx);
            end
            HOS.E_XF3min = TPinv*E_XF3;
        else
            HOS.E_XF3min = [];
        end
        if setup.HOSFourthOrder
            PermMat_z = PermutMat(nz); % Permutation matrix for fourth-order cumulant    
            auxCUM4z = (speye(nz^4)+transpose(PermMat_z)+PermMat_z)*kron(C2z0(:),C2z0(:));    
            E_XF4 = C4z0(:)+auxCUM4z;
            try 
                load(sprintf('QP_%i',nx),'QPinv'); 
            catch
                [QP,QPinv] = quadruplication(nx,1); 
                save(sprintf('QP_%i',nx),'QP','QPinv'); 
                fprintf('You may copy QP_%i.mat into the HOSScriptFiles/Other folder\n',nx);
            end
            HOS.E_XF4min = QPinv*E_XF4;
        else
            HOS.E_XF4min=[];        
        end
    elseif setup.orderApp == 2
        if setup.sparseflag
            Inx = speye(nx);
        else
            Inx = eye(nx);
        end
        % From First-order Approximation
        HOS.E_XF1 = E_XF1;
        HOS.E_XF2min = E_XF2min;    
        HOS.E_XF3min = E_XF3min;
        HOS.E_XF4min = E_XF4min;
        try 
            load(sprintf('DP_%i',nx),'DP','DPinv'); 
        catch
            [DP,DPinv] = duplication(nx); 
            save(sprintf('DP_%i',nx),'DP','DPinv'); 
            fprintf('You may copy DP_%i.mat into the HOSScriptFiles/Other folder\n',nx);
        end
        try 
            load(sprintf('TP_%i',nx),'TP','TPinv'); 
        catch
            [TP,TPinv] = triplication(nx,1); 
            save(sprintf('TP_%i',nx),'TP','TPinv');
            fprintf('You may copy TP_%i.mat into the HOSScriptFiles/Other folder\n',nx);
        end
        try 
            load(sprintf('QP_%i',nx),'QP','QPinv'); 
        catch
            [QP,QPinv] = quadruplication(nx,1);
            save(sprintf('QP_%i',nx),'QP','QPinv');
            fprintf('You may copy QP_%i.mat into the HOSScriptFiles/Other folder\n',nx);
        end
        E_XF2 = DP*E_XF2min;
        E_XF3 = TP*E_XF3min;
        E_XF4 = QP*E_XF4min;
        % From expectation    
        E_XS1 = Ez(id_xs);
        % From C2z0
        if setup.dimreduce
            Fz_Fz = kron(Fz,Fz);
            C2z0 = Fz_Fz*C2z0(:);
            nz = size(Fz,1);
        end
        E_XS2     = GetKronMom2(GetSubMoms2(C2z0(:),id_xs,id_xs,nz),E_XS1,E_XS1);
        E_XF1_XS1 = GetKronMom2(GetSubMoms2(C2z0(:),id_xf,id_xs,nz),E_XF1,E_XS1);
        E_XF2_XS1 = GetKronMom2(GetSubMoms2(C2z0(:),id_xf_xf,id_xs,nz),E_XF2,E_XS1);
        HOS.E_XS1 = E_XS1;
        HOS.E_XS2min = DPinv*E_XS2;
        HOS.E_XF1_XS1 = E_XF1_XS1;
        HOS.E_XF2_XS1min = kron(DPinv,Inx)*E_XF2_XS1;
        if setup.HOSThirdOrder || setup.HOSFourthOrder
            try 
                load(sprintf('QP5_%i',nx),'QP5inv'); 
            catch
                [QP5,QP5inv] = Q5_plication(nx,1); 
                save(sprintf('QP5_%i',nx),'QP5','QP5inv'); 
                fprintf('You may copy QP5_%i.mat into the HOSScriptFiles/Other folder\n',nx);
            end
            try 
                load(sprintf('QP6_%i',nx),'QP6inv'); 
            catch
                [QP6,QP6inv] = Q6_plication(nx,1); 
                save(sprintf('QP6_%i',nx),'QP6','QP6inv'); 
                fprintf('You may copy QP6_%i.mat into the HOSScriptFiles/Other folder\n',nx);
            end
            % From C3z0 
            if setup.dimreduce
                Fz_Fz_Fz = kron(Fz,Fz_Fz);
                C3z0 = Fz_Fz_Fz*C3z0(:);
            end
            E_XF5     = GetKronMom3(GetSubMoms3(C3z0(:),id_xf_xf,id_xf_xf,id_xf,nz),E_XF2,E_XF2,E_XF1,E_XF4,E_XF3,E_XF3);
            HOS.E_XF5min = QP5inv*E_XF5;
            E_XF6     = GetKronMom3(GetSubMoms3(C3z0(:),id_xf_xf,id_xf_xf,id_xf_xf,nz),E_XF2,E_XF2,E_XF2,E_XF4,E_XF4,E_XF4);
            HOS.E_XF6min = QP6inv*E_XF6;
            E_XS3     = GetKronMom3(GetSubMoms3(C3z0(:),id_xs,id_xs,id_xs,nz),E_XS1,E_XS1,E_XS1,E_XS2,E_XS2,E_XS2);
            HOS.E_XS3min = TPinv*E_XS3;
            E_XF1_XS2 = GetKronMom3(GetSubMoms3(C3z0(:),id_xf,id_xs,id_xs,nz),E_XF1,E_XS1,E_XS1,E_XF1_XS1,E_XF1_XS1,E_XS2);
            HOS.E_XF1_XS2min = kron(Inx,DPinv)*E_XF1_XS2;
            E_XF2_XS2 = GetKronMom3(GetSubMoms3(C3z0(:),id_xf_xf,id_xs,id_xs,nz),E_XF2,E_XS1,E_XS1,E_XF2_XS1,E_XF2_XS1,E_XS2);
            HOS.E_XF2_XS2min = kron(DPinv,DPinv)*E_XF2_XS2;
            E_XF3_XS1 = GetKronMom3(GetSubMoms3(C3z0(:),id_xf_xf,id_xf,id_xs,nz),E_XF2,E_XF1,E_XS1,E_XF3,E_XF2_XS1,E_XF1_XS1);
            HOS.E_XF3_XS1min = kron(TPinv,Inx)*E_XF3_XS1;
            E_XF4_XS1 = GetKronMom3(GetSubMoms3(C3z0(:),id_xf_xf,id_xf_xf,id_xs,nz),E_XF2,E_XF2,E_XS1,E_XF4,E_XF2_XS1,E_XF2_XS1);
            HOS.E_XF4_XS1min = kron(QPinv,Inx)*E_XF4_XS1;
        else
            HOS.E_XF5min=[]; HOS.E_XF6min=[]; HOS.E_XS3min=[]; HOS.E_XF1_XS2min =[]; HOS.E_XF2_XS2min = []; HOS.E_XF3_XS1min = []; HOS.E_XF4_XS1min=[];
        end
        if setup.HOSFourthOrder
            % From C4z0
            if setup.dimreduce
                Fz_Fz_Fz_Fz = kron(Fz,Fz_Fz_Fz);
                C4z0 = Fz_Fz_Fz_Fz*C4z0(:);
            end
            PermMat_z = PermutMat(nz); % Permutation matrix for fourth-order cumulant
            auxCUM4z = (speye(nz^4)+transpose(PermMat_z)+PermMat_z)*kron(C2z0(:),C2z0(:));    
            Ezzzz = C4z0(:)+auxCUM4z;
            E_XS4     = GetKronMom4(GetSubMoms4(Ezzzz,id_xs,id_xs,id_xs,id_xs,nz),E_XS1,E_XS1,E_XS1,E_XS1,...
                                                E_XS2,E_XS2,E_XS2,E_XS2,E_XS2,E_XS2,...
                                                E_XS3,E_XS3,E_XS3,E_XS3);
            E_XF2_XS3 = GetKronMom4(GetSubMoms4(Ezzzz,id_xf_xf,id_xs,id_xs,id_xs,nz),E_XF2,E_XS1,E_XS1,E_XS1,...
                                                E_XF2_XS1,E_XF2_XS1,E_XF2_XS1,E_XS2,E_XS2,E_XS2,...
                                                E_XF2_XS2,E_XF2_XS2,E_XF2_XS2,E_XS3);
            E_XF3_XS2 = GetKronMom4(GetSubMoms4(Ezzzz,id_xf_xf,id_xf,id_xs,id_xs,nz),E_XF2,E_XF1,E_XS1,E_XS1,...
                                                E_XF3,E_XF2_XS1,E_XF2_XS1,E_XF1_XS1,E_XF1_XS1,E_XS2,...
                                                E_XF3_XS1,E_XF3_XS1,E_XF2_XS2,E_XF1_XS2);
            E_XF4_XS2 = GetKronMom4(GetSubMoms4(Ezzzz,id_xf_xf,id_xf_xf,id_xs,id_xs,nz),E_XF2,E_XF2,E_XS1,E_XS1,...
                                                E_XF4,E_XF2_XS1,E_XF2_XS1,E_XF2_XS1,E_XF2_XS1,E_XS2,...
                                                E_XF4_XS1,E_XF4_XS1,E_XF2_XS2,E_XF2_XS2);    
            E_XF1_XS3 = GetKronMom4(GetSubMoms4(Ezzzz,id_xf,id_xs,id_xs,id_xs,nz),E_XF1,E_XS1,E_XS1,E_XS1,...
                                                E_XF1_XS1,E_XF1_XS1,E_XF1_XS1,E_XS2,E_XS2,E_XS2,...
                                                E_XF1_XS2,E_XF1_XS2,E_XF1_XS2,E_XS3);    
            E_XF5_XS1 = GetKronMom4(GetSubMoms4(Ezzzz,id_xf_xf,id_xf_xf,id_xf,id_xs,nz),E_XF2,E_XF2,E_XF1,E_XS1,...
                                                E_XF4,E_XF3,E_XF2_XS1,E_XF3,E_XF2_XS1,E_XF1_XS1,...
                                                E_XF5,E_XF4_XS1,E_XF3_XS1,E_XF3_XS1); 
            E_XF6_XS1 = GetKronMom4(GetSubMoms4(Ezzzz,id_xf_xf,id_xf_xf,id_xf_xf,id_xs,nz),E_XF2,E_XF2,E_XF2,E_XS1,...
                                                E_XF4,E_XF4,E_XF2_XS1,E_XF4,E_XF2_XS1,E_XF2_XS1,...
                                                E_XF6,E_XF4_XS1,E_XF4_XS1,E_XF4_XS1); 
            E_XF7     = GetKronMom4(GetSubMoms4(Ezzzz,id_xf_xf,id_xf_xf,id_xf_xf,id_xf,nz),E_XF2,E_XF2,E_XF2,E_XF1,...
                                                E_XF4,E_XF4,E_XF3,E_XF4,E_XF3,E_XF3,...
                                                E_XF6,E_XF5,E_XF5,E_XF5); 
            E_XF8     = GetKronMom4(GetSubMoms4(Ezzzz,id_xf_xf,id_xf_xf,id_xf_xf,id_xf_xf,nz),E_XF2,E_XF2,E_XF2,E_XF2,...
                                                E_XF4,E_XF4,E_XF4,E_XF4,E_XF4,E_XF4,...
                                                E_XF6,E_XF6,E_XF6,E_XF6);
            try 
                load(sprintf('QP7_%i',nx),'QP7inv'); 
            catch
                [QP7,QP7inv] = Q7_plication(nx,1); 
                save(sprintf('QP7_%i',nx),'QP7','QP7inv');
                fprintf('You may copy QP7_%i.mat into the HOSScriptFiles/Other folder\n',nx);
            end
            try 
                load(sprintf('QP8_%i',nx),'QP8inv'); 
            catch
                [QP8,QP8inv] = Q8_plication(nx,1); 
                save(sprintf('QP8_%i',nx),'QP8','QP8inv');
                fprintf('You may copy QP8_%i.mat into the HOSScriptFiles/Other folder\n',nx);
            end
            HOS.E_XF7min = QP7inv*E_XF7;
            HOS.E_XF8min = QP8inv*E_XF8;
            HOS.E_XS4min = QPinv*E_XS4;
            HOS.E_XF1_XS3min = kron(Inx,TPinv)*E_XF1_XS3;
            HOS.E_XF2_XS1min = kron(DPinv,Inx)*E_XF2_XS1;
            HOS.E_XF2_XS3min  = kron(DPinv,TPinv)*E_XF2_XS3;
            HOS.E_XF3_XS2min = kron(TPinv,DPinv)*E_XF3_XS2;
            HOS.E_XF4_XS2min = kron(QPinv,DPinv)*E_XF4_XS2;
            HOS.E_XF5_XS1min = kron(QP5inv,Inx)*E_XF5_XS1;
            HOS.E_XF6_XS1min = kron(QP6inv,Inx)*E_XF6_XS1;
        else
            % From C4z0
            HOS.E_XS4min     = []; HOS.E_XF2_XS3min = []; HOS.E_XF3_XS2min = []; HOS.E_XF4_XS2min = [];    
            HOS.E_XF1_XS3min = []; HOS.E_XF5_XS1min = []; HOS.E_XF6_XS1min = []; HOS.E_XF7min     = []; HOS.E_XF8min     = [];
        end    
    end
end

end %MAIN function end
