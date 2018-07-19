% Original author: Martin M. Andreasen
% Modified by Willi Mutschler, September 26, 2016. Email: willi@mutschler.eu

function [paramsOpt,setup] = runOptimization(params0,setup)
if setup.transpar
    % Transform par from model(restricted) to max(unrestricted)
    params0Values = invtrans(params0,setup.lowerBounds,setup.upperBounds,setup.tr);
else
    params0Values = params0;
end

if setup.optim == 1
%     InsigmaValues      = struc2values(setup.Insigma,setup.selectParams);  %The std deviation in the initial search distributions
%     sigma              = 0.1;                          %The step size
%     opts.SigmaMax      = 1;                            %The maximal value for sigma
%     opts.LBounds       = setup.lowerBoundsValues;      %Lower bound for params
%     opts.UBounds       = setup.upperBoundsValues;      %Upper bound for params
%     opts.MaxIter       = setup.MaxIter;                %The maximum number of iterations
%     opts.MaxFunEvals   = setup.MaxEvals;               %The maximum number of function evaluations
%     opts.PopSize       = setup.PopSize;                %The population size
%     opts.VerboseModulo = 1;                            %Display results after every 10'th iteration
%     opts.TolFun        = setup.TolFun;                 %Function tolerance
%     opts.TolX          = setup.TolX;                   %Tolerance in the parameters
%     opts.Plotting      = 'off';                        %Dislpay plotting or not
%     opts.Saving        = 'off';                        %Saving results
%     paramsOptValuesCMEAS = cmaes_dsgeDisplay(@objectFunc,params0Values,sigma,InsigmaValues,opts,setup,AUX);    
%     paramsOptValues = paramsOptValuesCMEAS;
elseif setup.optim == 2
    options = optimset('Display','iter','MaxIter',setup.MaxIter,'MaxFunEvals',setup.MaxEvals,'TolFun',setup.TolFun,'TolX',setup.TolX);   
    paramsOptValuesGRAD = lsqnonlin(@objectFunc,params0Values,[],[],options,setup);
    paramsOptValues = paramsOptValuesGRAD;
elseif setup.optim == 3
%     options = optimset('Display','iter','MaxIter',setup.MaxIter,'MaxFunEvals',setup.MaxEvals,'TolFun',setup.TolFun,'TolX',setup.TolX);
%     paramsOptValuesNM = fminsearch(@objectFunc,params0Values,options,setup,AUX);
%     paramsOptValues = paramsOptValuesNM;
end

if setup.transpar
    % Transform par from max(unrestricted) to model(restricted)
    paramsOpt = trans(paramsOptValues,setup.lowerBounds,setup.upperBounds,setup.tr);
else
    paramsOpt = paramsOptValues;
end


end