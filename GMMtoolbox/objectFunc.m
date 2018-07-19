% Original author: Martin M. Andreasen
% Modified by Willi Mutschler, September 26, 2016. Email: willi@mutschler.eu

function [out,momentsModel,info,Q] = objectFunc(paramsInputValues,setup)
% This function evaluates the objective function for GMM estimation
% info = 0, everything is ok, else for reporting problems)
% Initializing the outputs
global oo_ M_ options_ % get Dynare structures; used to pass them on to resol.m
out = [];
momentsModel = [];

% For efficient starting values when computing unconditional second moments
global Var_z2old Var_z3old 
global C2z0_old_2 C3z0_old_1 C4z0_old_1 C2z0_old_1 C3z0_old_2 C4z0_old_2 C2z0_old_3 C3z0_old_3 C4z0_old_3
if ~exist('Var_z2old','var'); Var_z2old = []; end
if ~exist('Var_z3old','var'); Var_z3old = []; end

if ~exist('C2z0_old_1','var'); C2z0_old_1 = []; end
if ~exist('C3z0_old_1','var'); C3z0_old_1 = []; end
if ~exist('C4z0_old_1','var'); C4z0_old_1 = []; end
if ~exist('C2z0_old_2','var'); C2z0_old_2 = []; end
if ~exist('C3z0_old_2','var'); C3z0_old_2 = []; end
if ~exist('C4z0_old_2','var'); C4z0_old_2 = []; end
if ~exist('C2z0_old_3','var'); C2z0_old_3 = []; end
if ~exist('C3z0_old_3','var'); C3z0_old_3 = []; end
if ~exist('C4z0_old_3','var'); C4z0_old_3 = []; end

if setup.transpar
    % Transform par from max(unrestricted) to model(restricted)
    paramsInputValues = trans(paramsInputValues,setup.lowerBounds,setup.upperBounds,setup.tr);
end

% Check if the lower and upper bounds are violated
if any(paramsInputValues < setup.lowerBounds) || any(paramsInputValues > setup.upperBounds)
    info = 1;
    if setup.optim == 1 || setup.optim == 3
        out = NaN;
    elseif setup.optim == 2
        out = ones(size(setup.dataMoments,1),1)*1D10;
    end
    return;
end

%% set parameter for use in Dynare
M_.params(setup.indx_params0) = paramsInputValues;

% Solve the model
[oo_.dr,info,M_,options_,oo_] = resol(0,M_,options_,oo_);

if any(info >0)
    if setup.optim == 1 || setup.optim == 3
        out = NaN;
    elseif setup.optim == 2
        out = ones(size(setup.dataMoments,1),1)*1D10;
    end
    return;
end


if setup.MomComp == 0
    % % Transformation of the approximated solution
% % We load and unfold the output from Dynare
% f_11 = [oo_.dr.ghx oo_.dr.ghu];
% [f_ss,fv,fvv,fss,fvvv,fssv,fsss,sigma,nv,nu,nx,nz] = Dynare_Unfold_Matrices(setup.orderApp,oo_,M_,f_11);
% 
% % We set up the alternative state space representation
% [gSteadyState_tmp,gv_tmp,gvv_tmp,gss_tmp,gvvv_tmp,gssv_tmp,gsss_tmp,hSteadyState,hv,hvv,hss,hvvv,hssv,hsss,eta,ny_tmp,label_y_tmp,label_v,label_u] = ...
%     StateSpaceDynare_LinearInov(f_ss,fv,fvv,fss,fvvv,fssv,fsss,sigma,nv,nx,nz,M_,oo_);
% 
% % We only use the selected variables in g(v)
% [ny,label_y,gSteadyState,gv,gvv,gss,gvvv,gssv,gsss] = ...
%     reduceControlsDynare(setup.selectY,nv,ny_tmp,label_y_tmp,gSteadyState_tmp,gv_tmp,gvv_tmp,gss_tmp,gvvv_tmp,gssv_tmp,gsss_tmp);
% clear gSteadyState_tmp gv_tmp gvv_tmp gss_tmp gvvv_tmp gssv_tmp gsss_tmp ny_tmp label_y_tmp
    % The moments in the normal distribution
%     ne         = size(model.eta,2);       %number of shocks
%     vectorMom3 = zeros(1,ne);
%     vectorMom4 = 3*ones(1,ne);
%     vectorMom5 = zeros(1,ne);
%     vectorMom6 = 15*ones(1,ne);
% 
%     % Computing moments
%     sig         = 1;
%     maxAuto_lag = max(setup.autoLagsIdx);
%     mom2nd      = UnconditionalMoments_2nd_Lyap(model.gx,model.gxx,model.gss,...
%                   model.hx,model.hxx,model.hss,model.eta,sig,vectorMom3,vectorMom4,maxAuto_lag,Var_z2old);
%     Var_z2old = mom2nd.Var_z;
%     if setup.orderApp == 3
%         mom3rd = UnconditionalMoments_3rd_Lyap(model.gx,model.gxx,model.gss,...
%               model.gxxx,model.gssx,model.gsss,...
%               model.hx,model.hxx,model.hss,model.hxxx,model.hssx,model.hsss,model.eta,sig,...
%               vectorMom3,vectorMom4,vectorMom5,vectorMom6,mom2nd.Var_xs,mom2nd.Var_xfxf,mom2nd.Var_xs_xfxf,maxAuto_lag,Var_z3old);
%           Var_z3old = mom3rd.Var_z;
%     end
%     % The moments implied by the model    
%     if setup.orderApp == 3
%         momentsModel =  momentsGMM(model,mom3rd,setup.autoLagsIdx,setup.selectHOS,0);        
%     else
%         momentsModel =  momentsGMM(model,mom2nd,setup.autoLagsIdx,setup.selectHOS,0);        
%     end
elseif setup.MomComp == 1
    if setup.orderApp == 1
        HOS = HOSrun(setup,[],C2z0_old_1,C3z0_old_1,C4z0_old_1,0);
        C2z0_old_1 = HOS.C2z0;
        if setup.HOSThirdOrder; C3z0_old_1 = HOS.C3z0; end
        if setup.HOSFourthOrder; C4z0_old_1 = HOS.C4z0; end
    elseif setup.orderApp == 2
        oldThirdOrder = setup.HOSThirdOrder; oldFourthOrder = setup.HOSFourthOrder;
        setup.orderApp = 1;
        setup.HOSThirdOrder = 1; setup.HOSFourthOrder = 1;
        HOS1st = HOSrun(setup,[],C2z0_old_1,C3z0_old_1,C4z0_old_1,1);        
        C2z0_old_1 = HOS1st.C2z0; C3z0_old_1 = HOS1st.C3z0; C4z0_old_1 = HOS1st.C4z0;
        setup.orderApp = 2;
        setup.HOSThirdOrder = oldThirdOrder; setup.HOSFourthOrder = oldFourthOrder;
        HOS2nd = HOSrun(setup,HOS1st,C2z0_old_2,C3z0_old_2,C4z0_old_2,0);
        C2z0_old_2 = HOS2nd.C2z0; 
        if setup.HOSThirdOrder; C3z0_old_2 = HOS2nd.C3z0; end
        if setup.HOSFourthOrder; C4z0_old_2 = HOS2nd.C4z0; end
        HOS = HOS2nd;
    elseif setup.orderApp == 3
        oldThirdOrder = setup.HOSThirdOrder; oldFourthOrder = setup.HOSFourthOrder;        
        setup.HOSThirdOrder = 1; setup.HOSFourthOrder = 1;         
        
        setup.orderApp = 1;
        HOS1st = HOSrun(setup,[],C2z0_old_1,C3z0_old_1,C4z0_old_1,1);        
        C2z0_old_1 = HOS1st.C2z0; 
        if setup.HOSThirdOrder; C3z0_old_1 = HOS1st.C3z0; end
        if setup.HOSFourthOrder;C4z0_old_1 = HOS1st.C4z0; end
        
        setup.orderApp = 2;
        HOS2nd = HOSrun(setup,HOS1st,C2z0_old_2,C3z0_old_2,C4z0_old_2,1);        
        C2z0_old_2 = HOS2nd.C2z0; 
        if setup.HOSThirdOrder; C3z0_old_2 = HOS2nd.C3z0; end
        if setup.HOSFourthOrder;C4z0_old_2 = HOS2nd.C4z0; end
        
        setup.orderApp = 3;
        setup.HOSThirdOrder = oldThirdOrder; setup.HOSFourthOrder = oldFourthOrder;
        HOS3rd = HOSrun(setup,HOS2nd,C2z0_old_3,C3z0_old_3,C4z0_old_3,0);        
        C2z0_old_3 = HOS3rd.C2z0;        
        if setup.HOSThirdOrder; C3z0_old_3 = HOS3rd.C3z0; end
        if setup.HOSFourthOrder; C4z0_old_3 = HOS3rd.C4z0; end
        HOS = HOS3rd; 
    end
    momentsModel =  momentsGMM(HOS,setup.autoLagsIdx,setup.HOSThirdOrder,setup.HOSFourthOrder,1);
end
% To Compare results with andreasen
% format long
% [mom2nd.E_y+model.g0 HOS2nd.Ey]
% [mom2nd.Var_y(:) HOS2nd.C2y0]
% [mom2nd.Cov_y(:) HOS2nd.C2yt]
% [mom3rd.E_y+model.g0 HOS3rd.Ey]
% [mom3rd.Var_y(:) HOS3rd.C2y0]
% [mom3rd.Cov_y(:) HOS3rd.C2yt]
% format short


% The moments implied by the model
residuals    = setup.Sw*(setup.dataMoments-momentsModel);
%Q = (dataMoments-momentsModel)'*(Sw'*Sw)*(dataMoments-momentsModel);
Q = residuals'*residuals;
if setup.optim == 1 || setup.optim == 3
    out = residuals'*residuals;
elseif setup.optim == 2
    out = residuals;
end
end

