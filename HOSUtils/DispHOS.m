% By Willi Mutschler, September 26, 2016. Email: willi@mutschler.eu
function DispHOS(opt)
global M_ oo_ options_


%% Section 1: Set default options
% Check if default options are overwritten by opt
opt = updateOptions(opt);
% Read options from Dynare
opt.orderApp = options_.order;
opt.indx_obs = oo_.dr.order_var; %note that we display all variables
opt.indx_states = transpose(M_.nstatic+(1:M_.nspred));
opt.nu = M_.exo_nbr; % number of shocks and measurement errors
opt.nx = M_.nspred; %number of states
opt.ny = M_.endo_nbr; %number of all endogenous variables

 
%% Section 3: Some checks and precomputations
CheckExistMoments(opt.orderApp,opt.HOSThirdOrder,opt.HOSFourthOrder,opt.distrib,M_)


%% Print theoretical statistics
% For efficient starting values when computing unconditional second moments
global C2z0_old_2 C3z0_old_1 C4z0_old_1 C2z0_old_1 C3z0_old_2 C4z0_old_2 C2z0_old_3 C3z0_old_3 C4z0_old_3
if ~exist('C2z0_old_1','var'); C2z0_old_1 = []; end
if ~exist('C3z0_old_1','var'); C3z0_old_1 = []; end
if ~exist('C4z0_old_1','var'); C4z0_old_1 = []; end
if ~exist('C2z0_old_2','var'); C2z0_old_2 = []; end
if ~exist('C3z0_old_2','var'); C3z0_old_2 = []; end
if ~exist('C4z0_old_2','var'); C4z0_old_2 = []; end
if ~exist('C2z0_old_3','var'); C2z0_old_3 = []; end
if ~exist('C3z0_old_3','var'); C3z0_old_3 = []; end
if ~exist('C4z0_old_3','var'); C4z0_old_3 = []; end

timerHOS = tic;
if opt.orderApp == 1
    HOS = HOSrun(opt,[],C2z0_old_1,C3z0_old_1,C4z0_old_1,0);
    C2z0_old_1 = HOS.C2z0;
    if opt.HOSThirdOrder; C3z0_old_1 = HOS.C3z0; end
    if opt.HOSFourthOrder; C4z0_old_1 = HOS.C4z0; end
elseif opt.orderApp == 2
    oldThirdOrder = opt.HOSThirdOrder; oldFourthOrder = opt.HOSFourthOrder;
    opt.orderApp = 1;
    opt.HOSThirdOrder = 1; opt.HOSFourthOrder = 1;
    HOS1st = HOSrun(opt,[],C2z0_old_1,C3z0_old_1,C4z0_old_1,1);        
    C2z0_old_1 = HOS1st.C2z0; C3z0_old_1 = HOS1st.C3z0; C4z0_old_1 = HOS1st.C4z0;
    opt.orderApp = 2;
    opt.HOSThirdOrder = oldThirdOrder; opt.HOSFourthOrder = oldFourthOrder;
    HOS2nd = HOSrun(opt,HOS1st,C2z0_old_2,C3z0_old_2,C4z0_old_2,0);
    C2z0_old_2 = HOS2nd.C2z0; 
    if opt.HOSThirdOrder; C3z0_old_2 = HOS2nd.C3z0; end
    if opt.HOSFourthOrder; C4z0_old_2 = HOS2nd.C4z0; end
    HOS = HOS2nd;
elseif opt.orderApp == 3
    oldThirdOrder = opt.HOSThirdOrder; oldFourthOrder = opt.HOSFourthOrder;        
    opt.HOSThirdOrder = 1; opt.HOSFourthOrder = 1;         

    opt.orderApp = 1;
    HOS1st = HOSrun(opt,[],C2z0_old_1,C3z0_old_1,C4z0_old_1,1);        
    C2z0_old_1 = HOS1st.C2z0; 
    if opt.HOSThirdOrder; C3z0_old_1 = HOS1st.C3z0; end
    if opt.HOSFourthOrder;C4z0_old_1 = HOS1st.C4z0; end

    opt.orderApp = 2;
    HOS2nd = HOSrun(opt,HOS1st,C2z0_old_2,C3z0_old_2,C4z0_old_2,1);        
    C2z0_old_2 = HOS2nd.C2z0; 
    if opt.HOSThirdOrder; C3z0_old_2 = HOS2nd.C3z0; end
    if opt.HOSFourthOrder;C4z0_old_2 = HOS2nd.C4z0; end

    opt.orderApp = 3;
    opt.HOSThirdOrder = oldThirdOrder; opt.HOSFourthOrder = oldFourthOrder;
    HOS3rd = HOSrun(opt,HOS2nd,C2z0_old_3,C3z0_old_3,C4z0_old_3,0);        
    C2z0_old_3 = HOS3rd.C2z0;        
    if opt.HOSThirdOrder; C3z0_old_3 = HOS3rd.C3z0; end
    if opt.HOSFourthOrder; C4z0_old_3 = HOS3rd.C4z0; end
    HOS = HOS3rd; 
end

% Theoretical variance, skewness and excess kurtosis of innovations and variables
nximin = sqrt(size(HOS.GAMMA2XI,1));
ny = length(opt.indx_obs);
nu = opt.nu;
[D2Gam,D3Gam,D4Gam] = DiagionalizationMatrix(nximin,1,1,1);
[D2y,D3y,D4y] = DiagionalizationMatrix(ny,1,1,1);

VARxi = D2Gam'*HOS.GAMMA2XI;
VARy = D2y'*HOS.C2y0;

if opt.HOSThirdOrder
    SKEWxi = (D3Gam'*HOS.GAMMA3XI)./VARxi.^1.5;
    SKEWy = (D3y'*HOS.C3y0)./VARy.^1.5;
else
    SKEWxi = nan(nximin,1);
    SKEWy  = nan(ny,1);
end
if opt.HOSFourthOrder
    KURTxi = (D4Gam'*HOS.GAMMA4XI)./(VARxi.*VARxi); %excess kurtosis
    KURTy = (D4y'*HOS.C4y0)./(VARy.*VARy); % excess kurtosis
else
    KURTxi = nan(nximin,1);
    KURTy = nan(ny,1);
end
elapsedTimeHOS = toc(timerHOS);
fprintf('Elapsed time for theoretical HOS computations is %.4f seconds.\n',elapsedTimeHOS)

if opt.simdata==0    
    oldnoprint = options_.noprint;
    options_.noprint = 0;
    labels = char(M_.exo_names);    
    zzz = [ zeros(nu,1) full(VARxi(1:nu)) full(SKEWxi(1:nu)) full(KURTxi(1:nu)) ];        
    title=sprintf('THEORETICAL STATISTICS OF INNOVATIONS FOR %s DISTRIBUTION',upper(opt.distrib));
    headers=char('VARIABLE','MEAN','VARIANCE','SKEWNESS','EXC.KURTOSIS');
    dyntable(title,headers,labels,zzz,size(labels,2)+2,16,6);
    labels = deblank(M_.endo_names(oo_.dr.order_var,:));
    zzz = [full(oo_.dr.ys(oo_.dr.order_var)) full(HOS.Ey) full(VARy) full(SKEWy) full(KURTy)];        
    title=sprintf('THEORETICAL STATISTICS OF VARIABLES FOR %s DISTRIBUTION',upper(opt.distrib));
    headers=char('VARIABLE','STEADY-STATE','MEAN','VARIANCE','SKEWNESS','EXC.KURTOSIS');
    dyntable(title,headers,labels,zzz,size(labels,2)+2,16,6);
    options_.noprint = oldnoprint;
elseif opt.simdata == 1
    oldnoprint = options_.noprint;
    options_.noprint = 0;
    %% Simulate data  
    R = opt.MCruns; numSim = opt.numSim; burnin = opt.burnin;
    numSim = 2*ceil(numSim/2);% make sure numSim is even
    TotSim = 2*ceil((numSim+burnin)/2);% make sure TotSim is even
    fprintf('Simulate data and compute empirical variance, skewness and excess kurtosis for innovations and observables\n'); 
    reverseStr = '';
    ny = length(opt.indx_obs);
    meanU = nan(nu,R); varU =  nan(nu,R); skewU =  nan(nu,R); kurtU =  nan(nu,R);
    meanY = nan(ny,R); varY =  nan(ny,R); skewY =  nan(ny,R); kurtY =  nan(ny,R);
    timerData = tic;
    for r = 1:R
        if r == 1 || rem(r,10) == 0
            msg = sprintf('   MC Sample %d of %d with %d observations each',r,R,numSim); fprintf([reverseStr, msg]); reverseStr = repmat(sprintf('\b'), 1, length(msg));
        elseif r==R
            msg = sprintf('   MC Sample %d of %d with %d observations each\n',r,R,numSim); fprintf([reverseStr, msg]); reverseStr = repmat(sprintf('\b'), 1, length(msg));                
        end        
        if strcmp(opt.distrib,'Gaussian')
            exo_shocks = transpose(mvnrnd(zeros(M_.exo_nbr,1),M_.Sigma_e,TotSim));    
        elseif strcmp(opt.distrib,'Student_t')
            for i = 1:M_.param_nbr
                if strcmp(deblank(M_.param_names(i,:)),'df_studt')
                    df_studt = M_.params(i);                
                end
            end
            V = gamrnd(df_studt/2,2/df_studt,TotSim,1);        
            U = mvnrnd(zeros(M_.exo_nbr,1),(df_studt-2)/df_studt*M_.Sigma_e,TotSim);
            exo_shocks = transpose(repmat(V.^(-1/2),1,M_.exo_nbr).*U);        
        end;
        if opt.antithetic            
            exo_shocks(:,TotSim/2+1:TotSim) = -exo_shocks(:,1:TotSim/2); % Antithetic shock generation - match skewness
            for i=1:M_.exo_nbr
                exo_shocks(i,:) = (exo_shocks(i,:) - mean(exo_shocks(i,:)'))/std(exo_shocks(i,:)')*sqrt(M_.Sigma_e(i,i)); % Quadratic resampling - match variance
            end            
        end
        Y_sim = transpose(simult_(oo_.dr.ys,oo_.dr,exo_shocks',opt.orderApp));
        % Get rid of burnin phase
        Y_sim(1:burnin,:) = [];
        Y_sim = Y_sim(:,oo_.dr.order_var);
        % Compute statistics        
        meanY(:,r) = mean(Y_sim); meanU(:,r) = mean(exo_shocks');
        varY(:,r) = var(Y_sim);   varU(:,r) = var(exo_shocks');
        if opt.HOSThirdOrder            
            skewY(:,r) = skewness(Y_sim);
            skewU(:,r) = skewness(exo_shocks');
        end
        if opt.HOSFourthOrder            
            kurtY(:,r) = kurtosis(Y_sim) - 3;
            kurtU(:,r) = kurtosis(exo_shocks') - 3;
        end
    end
    
    % Empirical variance, skewness and excess kurtosis of innovations
    meanu = mean(meanU,2);  sdmeanu = std(meanU')';  
    varu = mean(varU,2);    sdvaru = std(varU')';    
    skewu = mean(skewU,2);  sdskewu = std(skewU')';
    kurtu = mean(kurtU,2);  sdkurtu = std(kurtU')';
    % Empirical mean,variance, skewness and excess kurtosis of observables
    meany = mean(meanY,2); sdmeany = std(meanY')';
    vary = mean(varY,2);   sdvary = std(varY')';
    skewy = mean(skewY,2); sdskewy = std(skewY')';
    kurty = mean(kurtY,2); sdkurty = std(kurtY')';
    elapsedTimeData = toc(timerData);    

    % Compare to theoretical values
    labels = char(M_.exo_names);
    zzz = [ zeros(nu,1) meanu sdmeanu full(VARxi(1:nu)) varu sdvaru full(SKEWxi(1:nu)) skewu sdskewu full(KURTxi(1:nu)) kurtu sdkurtu];        
    title=sprintf('EMPIRICAL (E) AND THEORETICAL (T) STATISTICS OF INNOVATIONS FOR %s DISTRIBUTION',upper(opt.distrib));
    headers=char('VARIABLE','MEAN(T)','MEAN(E)','std(MEAN(E))','VARIANCE(T)','VARIANCE(E)','std(VARIANCE(E))','SKEWNESS(T)','SKEWNESS(E)','std(SKEWNESS(E))','EXC.KURTOSIS(T)','EXC.KURTOSIS(E)','std(EXC.KURTOSIS(E))');
    dyntable(title,headers,labels,zzz,size(labels,2)+2,16,6);
    labels = deblank(M_.endo_names(oo_.dr.order_var,:));
    zzz = [ full(HOS.Ey) meany sdmeany full(VARy) vary sdvary full(SKEWy) skewy sdskewy full(KURTy) kurty sdkurty];    
    title=sprintf('EMPIRICAL (E) AND THEORETICAL (T) STATISTICS OF VARIABLES FOR %s DISTRIBUTION',upper(opt.distrib));
    headers=char('VARIABLE','MEAN(T)','MEAN(E)','std(MEAN(E))','VARIANCE(T)','VARIANCE(E)','std(VARIANCE(E))','SKEWNESS(T)','SKEWNESS(E)','std(SKEWNESS(E))','EXC.KURTOSIS(T)','EXC.KURTOSIS(E)','std(EXC.KURTOSIS(E))');
    dyntable(title,headers,labels,zzz,size(labels,2)+2,16,6);
    fprintf('Elapsed time for data simulation is %.4f seconds.\n',elapsedTimeData)
    options_.noprint = oldnoprint;
end

end