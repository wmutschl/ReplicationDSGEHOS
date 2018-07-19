% By Willi Mutschler, September 26, 2016. Email: willi@mutschler.eu
function [results,params0Step2] = RunGMM(opt)
global M_ oo_ options_ estim_params_
% This script estimates DSGE models solved up to third order by GMM
% using the following moments
%  - E[y]
%  - E[y(i)*y(j)]            for i=1:ny and j=i,ny
%  - E[y(i)_t*y(i)_t-k]      for i=1:ny and k=1,2,...autoLagsIdx
%  - E[y(i)*y(i)*y(i)]       for i=1:ny
%  - E[y(i)*y(i)*y(i)*y(i)]  for i=1:ny
% shocks can be either multiariate Gaussian or Student-t distributed
%

%% Section 1: Set default options
% Check if default options are overwritten by opt
opt = updateOptions(opt);
opt.orderApp              = options_.order;
opt.indx_obs              = options_.varobs_id;
opt.indx_states = transpose(M_.nstatic+(1:M_.nspred));
opt.nu = M_.exo_nbr; % number of shocks and measurement errors
opt.nx = M_.nspred; %number of states
opt.ny = M_.endo_nbr; %number of all endogenous variables
opt.indx_params0 = estim_params_.param_vals(:,1);
opt.calibrateParams = M_.params(setdiff(1:M_.param_nbr,estim_params_.param_vals(:,1)));
opt.params0 = estim_params_.param_vals(:,2);
opt.lowerBounds = estim_params_.param_vals(:,3);
opt.upperBounds = estim_params_.param_vals(:,4);
opt.tr =  estim_params_.param_vals(:,10);
if isempty(options_.qz_criterium)
    options_.qz_criterium = 1-1e-6;
end
 
%% Section 3: Some checks and precomputations
if opt.HOSThirdOrder || opt.HOSFourthOrder
    opt.MomComp = 1; % Andreasen's method does not compute third- and fourth order moments
end
CheckExistMoments(opt.orderApp,opt.HOSThirdOrder,opt.HOSFourthOrder,opt.distrib,M_);

%% Section 4: We simulate the model using M_.params or load data
if isempty(opt.seed_nr) == 0
    randn('seed',opt.seed_nr);
end
opt.numSim = 2*ceil(opt.numSim/2);% make sure numSim is even
TotSim = 2*ceil((opt.numSim+opt.burnin)/2);% make sure TotSim is even
if strcmp(opt.datadistrib,'Gaussian');    
    exo_shocks = transpose(mvnrnd(zeros(M_.exo_nbr,1),M_.Sigma_e,TotSim));
elseif strcmp(opt.datadistrib,'Student_t')
    for i = 1:M_.param_nbr
        try
            df_studt = opt.df_studt;
        catch
            if strcmp(deblank(M_.param_names(i,:)),'df_studt')
                df_studt = M_.params(i);                
            end
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
opt.data = Y_sim(opt.burnin+1:end-1,opt.indx_obs);
% Compute empirical moments
[opt.dataMoments,opt.nameMoments] = momentsGMMData(opt.data,opt.autoLagsIdx,opt.HOSThirdOrder,opt.HOSFourthOrder);

disp(['Parameters to estimate = ', num2str(length(opt.params0)) ,'. Moments for estimation = ', num2str(size(opt.dataMoments,1))]);
if size(opt.dataMoments,1) < length(opt.params0)
    error('We must have at least as many moments as parameters for GMM')    
end


%% Section 5: Constructing the struct setupStep1
setupStep1 = opt;

%% AUXILIARY TEST IF EMPIRICAL MOMENTS FIT TO THEORETICAL MOMENTS
if opt.AUXTEST
    setupStep1.params0 = M_.params(opt.indx_params0);
    setupStep1.Wstep1  = getOptimalWeighting(opt.dataMoments,opt.data,opt.autoLagsIdx,opt.HOSThirdOrder,opt.HOSFourthOrder,opt.qLag);    
    tmphFunc    = get_hFunc(opt.data,opt.autoLagsIdx,opt.dataMoments,opt.HOSThirdOrder,opt.HOSFourthOrder);    
    tmp_momentsData = [mean(tmphFunc + repmat(opt.dataMoments',opt.numSim,1),1)'];
    setupStep1.Sw            = chol(diag(diag(eye(size(setupStep1.Wstep1)))));
    setup = setupStep1;    
    setup.lowerBoundsValues = setupStep1.lowerBounds;
    setup.upperBoundsValues = setupStep1.upperBounds;
    %AUX = GetAUX(options_.order,nx,ne);
    R = 1000; numSim = 1000; antithetic = 0;
    Z = zeros(size(setupStep1.dataMoments,1),R);
    reverseStr = '';
    for r=1:R
        if rem(r,10) == 0
            msg = sprintf('    Simulate data processed %d/%d', r, R); fprintf([reverseStr, msg]); reverseStr = repmat(sprintf('\b'), 1, length(msg));
                elseif r==R
            msg = sprintf('    Simulate data processed %d/%d \n', r, R); fprintf([reverseStr, msg]); reverseStr = repmat(sprintf('\b'), 1, length(msg));                
        end
        switch opt.distrib
            case 'Gaussian'
                exo_shocks = transpose(mvnrnd(zeros(M_.exo_nbr,1),M_.Sigma_e,numSim));        
            case 'Student_t'
                for i = 1:M_.param_nbr
                    if strcmp(deblank(M_.param_names(i,:)),'df_studt')
                        df_studt = M_.params(i);                
                    end
                end
                V = gamrnd(df_studt/2,2/df_studt,numSim,1);        
                U = mvnrnd(zeros(M_.exo_nbr,1),(df_studt-2)/df_studt*M_.Sigma_e,numSim);
                exo_shocks = transpose(repmat(V.^(-1/2),1,M_.exo_nbr).*U);
        end
        if antithetic
            exo_shocks(:,numSim/2+1:numSim) = -exo_shocks(:,1:numSim/2); % Antithetic shock generation - match skewness
            for i=1:M_.exo_nbr
                exo_shocks(i,:) = (exo_shocks(i,:) - mean(exo_shocks(i,:)'))/std(exo_shocks(i,:)')*sqrt(M_.Sigma_e(i,i)); % Quadratic resampling - match variance
            end
        end
        datz = transpose(simult_(oo_.dr.ys,oo_.dr,exo_shocks',options_.order));
        Z(:,r) = momentsGMMData(datz(2:end,opt.indx_obs),opt.autoLagsIdx,opt.HOSThirdOrder,opt.HOSFourthOrder);
    end
    fprintf('\n');
    tmp_momentsMC = mean(Z,2);    
    if setup.transpar
        % Transform par from model(restricted) to max(unrestricted)
        paramsInputValues = invtrans(setupStep1.params0,setup.lowerBounds,setup.upperBounds,setup.tr);
    else
        paramsInputValues           = setupStep1.params0;
    end
    
    [tmp_out,tmp_momentsModel,tmp_errorMes,tmp_Q] = objectFunc(paramsInputValues,setup);
    format long; disp([full(tmp_momentsData) full(setup.dataMoments) full(tmp_momentsMC) full(tmp_momentsModel)]); format short;

else

    %% Section 6: Step 1 of the GMM estimation
    setupStep1.optimWeightMat    = 0;
    Wstep1                       = getOptimalWeighting(opt.dataMoments,opt.data,opt.autoLagsIdx,opt.HOSThirdOrder,opt.HOSFourthOrder,opt.qLag);    
    setupStep1.Sw                = chol(diag(diag(Wstep1)));
    paramsStep1                  = opt.params0;
    for i=1:setupStep1.numOptimStep1
        [paramsStep1,setupStep1] = runOptimization(paramsStep1,setupStep1);        
    end
    results.Step1                 = getSEGMM(paramsStep1,setupStep1);%The standard errors
    results.paramsStep1 = paramsStep1;
    %% Section 7: Step 2 of the GMM estimation
    [~,modelMoments,~,~]         = objectFunc(paramsStep1,setupStep1);    
    Wopt                         = getOptimalWeighting(modelMoments,opt.data,opt.autoLagsIdx,opt.HOSThirdOrder,opt.HOSFourthOrder,opt.qLag);
    setupStep2                   = setupStep1;
    setupStep2.Sw                = chol(Wopt);
    setupStep2.optimWeightMat    = 1;
    disp(['Rank of Wopt = ', num2str(rank(Wopt)), ' Size of Wopt: ', num2str(size(Wopt))])
    params0Step2                 = paramsStep1;
    for i=1:setupStep2.numOptimStep2
        [paramsStep2,setupStep2] = runOptimization(params0Step2,setupStep2);
        params0Step2             = paramsStep2;
    end
    results.Step2                = getSEGMM(paramsStep2,setupStep2); %The standard errors    
    results.params0Step2         = params0Step2;
end