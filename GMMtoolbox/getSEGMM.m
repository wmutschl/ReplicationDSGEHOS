% Original author: Martin M. Andreasen
% Modified by Willi Mutschler, September 26, 2016. Email: willi@mutschler.eu

% This function computes standard errors to the GMM estimates
% using the optimal weighting matrix

function out = getSEGMM(paramsInputValues,setup)
epsValue = setup.epsValue;
% Evaluating the objective function to get modelMoments
[~,modelMoments,info,Q] = objectFunc(paramsInputValues,setup);

% Some dimensions
numMom      = size(modelMoments,1);
dimParams   = size(paramsInputValues,1);
D           = zeros(numMom,dimParams);

for i=1:dimParams
    %Positive step
    paramsInputValues_p_eps      = paramsInputValues;
    paramsInputValues_p_eps(i,1) = paramsInputValues_p_eps(i) + epsValue;
    [~,modelMoments_p_eps,errorMesP] = objectFunc(paramsInputValues_p_eps,setup);

    % Negative step
    paramsInputValues_m_eps      = paramsInputValues;
    paramsInputValues_m_eps(i,1) = paramsInputValues_m_eps(i) - epsValue;
    [~,modelMoments_m_eps,errorMesM] = objectFunc(paramsInputValues_m_eps,setup);

    % The Jacobian: 
    if errorMesP == 0 && errorMesM == 0        
        D(:,i) = (modelMoments_p_eps-modelMoments_m_eps)/(2*epsValue);
    else
        disp('Cannot compute the jacobiant - no standard errors available')
        out = NaN;
        return
    end
end

T        = size(setup.data,1); %Number of observations
if setup.optimWeightMat == 0
    % We do not have the optimal weighting matrix
    WWused        = setup.Sw'*setup.Sw;
    [~,modelMoments] = objectFunc(paramsInputValues,setup);
    WWopt         = getOptimalWeighting(modelMoments,setup.data,setup.autoLagsIdx,setup.HOSThirdOrder,setup.HOSFourthOrder,setup.qLag);    
    S             = WWopt\eye(size(WWopt,1));
    AA            = (D'*WWused*D)\eye(dimParams);
    AVar          = 1/T*AA*D'*WWused*S*WWused*D*AA;    
    % The J-test (see John Cochranes book on asset pricing Ch. 11)
    % Does not seem to work well at the moment, probably because the way
    % tempInv is computed!
    DWD           = D'*WWused*D;
    numMoments    = size(setup.Sw,1);
    temp          = (eye(numMoments)-D*DWD*D'*WWused)*S*(eye(numMoments)-D*DWD*D'*WWused)';
    tempInv       = pinv(temp);   % the a pseudo-inversion
    out.Jtest     = T*(setup.dataMoments-modelMoments)'*tempInv*(setup.dataMoments-modelMoments);
    out.JtestDf   = size(setup.Sw,1)-dimParams;
    out.ProbJtest = chi2cdf(out.Jtest,out.JtestDf,'upper');  
elseif setup.optimWeightMat == 1
    % We have the optimal weigthing matrix
    WW            = setup.Sw'*setup.Sw;
    AVar          = 1/T*((D'*WW*D)\eye(dimParams));
    % The J-test
    out.Jtest     = T*Q;
    out.JtestDf   = size(setup.Sw,1)-dimParams;
    out.ProbJtest = chi2cdf(out.Jtest,out.JtestDf,'upper');  
end
%out.params       = values2struct(paramsInputValues,setup.selectParams);
out.params       = paramsInputValues;

SEvalues         = sqrt(diag(AVar));
%out.paramsSE     = values2struct(SEvalues,setup.selectParams);
out.paramsSE     = SEvalues;
out.Jacobian     = D;
out.Q            = Q;
%out.model        = model;
out.modelMoments = modelMoments;
out.dataMoments  = setup.dataMoments;
out.nameMoments  = setup.nameMoments;


end