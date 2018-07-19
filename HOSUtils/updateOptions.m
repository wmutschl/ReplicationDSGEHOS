% By Willi Mutschler, September 26, 2016. Email: willi@mutschler.eu
function opt = updateOptions(opt)
% Number of lagged covariances included in the GMM estimation, e.g. [1 5] to include first and fifth lag
if ~isfield(opt,'autoLagsIdx')
    opt.autoLagsIdx = [1];
end
% Third-order product moments in the GMM estimation,  0: do not include or 1: include
if ~isfield(opt,'HOSThirdOrder')
    opt.HOSThirdOrder = 0;
end
% Fourth-order product moments in the GMM estimation, 0: do not include or 1: include
if ~isfield(opt,'HOSFourthOrder')
    opt.HOSFourthOrder = 0;
end
% Method to compute first two unconditional moments: 0 for Andreasen's, 1 for Mutschler's
if ~isfield(opt,'MomComp')
    opt.MomComp = 1;
end
% AUXILIARY TEST IF EMPIRICAL MOMENTS FIT TO THEORETICAL MOMENTS
if ~isfield(opt,'AUXTEST')
    opt.AUXTEST = 0;
end
% Reduce further dimension in ABCD System
if ~isfield(opt,'dimreduce')
    opt.dimreduce = 0;
end
%distributions supposed in the state space system: 'Gaussian' for multivariate normal distribution, 'Student_t' for multivariate student-t distribution
if ~isfield(opt,'distrib')
    opt.distrib = 'Gaussian';
end
%distribution used for data simulation: 'Gaussian' for multivariate normal distribution, 'Student_t' for multivariate student-t distribution
if ~isfield(opt,'datadistrib')
    opt.datadistrib = opt.distrib;
end

% Seed for simulated sample path
if ~isfield(opt,'seed_nr')
    opt.seed_nr = [];
end
% use antithetic shocks in simulated data to fit moments better
if ~isfield(opt,'antithetic')
    opt.antithetic = 0;
end
% Length of the simulated sample path
if ~isfield(opt,'numSim')
    opt.numSim = 250;
end
%Burn-in simulated data
if ~isfield(opt,'burnin')
    opt.burnin = 1000;
end
%Compare theoretical statistics with simulated data
if ~isfield(opt,'simdata')
    opt.simdata = 0;
end
%Monte Carlo Runs for comparision with simulated data
if ~isfield(opt,'MCruns')
    opt.MCruns = 1000;
end

% 1 for the CMAES, 2 for a gradient based optimizer, 3 for the simple algorithm by Nelder-Mead.
if ~isfield(opt,'optim')
    opt.optim = 2;
end
%Transform parameters
if ~isfield(opt,'transpar')
    opt.transpar = 0;
end

% Number of optimizations in step 1
if ~isfield(opt,'numOptimStep1')
    opt.numOptimStep1 = 2;
end
% Number of optimizations in step 1
if ~isfield(opt,'numOptimStep2')
    opt.numOptimStep2 = 2;
end
% Max number of iterations for the optimizer
if ~isfield(opt,'MaxIter')
    opt.MaxIter = 5000;
end
% Max number of function evaluations for the optimizer
if ~isfield(opt,'MaxEvals')
    opt.MaxEvals = 5000;
end
% Function tolerance at the objective function
if ~isfield(opt,'TolFun')
    opt.TolFun = 1D-6;
end
% Function tolerance for the coefficients
if ~isfield(opt,'TolX')
    opt.TolX = 1D-6;
end
% Number of draws per generation i.e. lambda in the CMAES 
if ~isfield(opt,'PopSize')
    opt.PopSize = 50;
end
% Use sparse matrices in state space representation
if ~isfield(opt,'sparseflag')
    opt.sparseflag = 0;
end
% Number of lags for the Newey-West estimator of W
if ~isfield(opt,'qLag')
    opt.qLag = 20;
end
% stepsize for computing standard errors
if ~isfield(opt,'epsValue')
    opt.epsValue = 1e-6;
end
%methods in solve_gensylv.m: 'disclyap' 'doubling' 'fixed-point', 'dlyap'.
if ~isfield(opt,'sylv_method')
    opt.sylv_method = 'disclyap';
end
%tolerance in solve_gensylv.m
if ~isfield(opt,'sylv_tol')
    opt.sylv_tol = 1e-16;
end
%maximum number of iterations in solve_gensylv.m
if ~isfield(opt,'sylv_max_iter')
    opt.sylv_max_iter = 500;
end
%use sparse matrices in solve_gensylv.m
if ~isfield(opt,'sylv_sparseflag')
    opt.sylv_sparseflag = 0;
end
