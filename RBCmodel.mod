%/////////////////////////////////////////////////////////////
%////    Model: RBC model with habit formation and variable labor                  ////
%////    Author: Willi Mutschler                          ////
%////    Email:  willi@mutschler.eu                       ////
%////    Version: August 26, 2016      		   		      ////
%/////////////////////////////////////////////////////////////
%----------------------------------------------------------------
% 0. Specify options: User Settings (For other settings and defaults see updateOptions.m
%----------------------------------------------------------------
% Order of approximation
@#define orderApp = 3
% Which distribution in state space system: 0:for Gaussian, 1:for Student_t
@#define distrib = 0

@#if distrib == 1
    opt.distrib = 'Student_t'; 
@#endif

opt.numSim               = 250;     % Number of observations of simulated data
opt.burnin               = 0;      % Initial burnin periods, will be discarded from simulation
opt.autoLagsIdx          = [1]      % Which autocovariances to include
opt.HOSThirdOrder        = 0;         % Use third-order statistics,  0: do not include or 1: include
opt.HOSFourthOrder       = 0;         % User fourth-order statistics, 0: do not include or 1: include
opt.seed_nr              = 10;
opt.transpar             = 0;

%----------------------------------------------------------------
% 1. Declare variables and parameters
%----------------------------------------------------------------
var k c a iv y la n rk w;
predetermined_variables k;
varexo u_a;
varobs c iv n;
parameters DELTA BETTA B ETAl ETAc THETA ALFA RHOA STDA;
@#if distrib == 1
    parameters df_studt; % Please always use df_studt for degrees of freedom
@# endif

%----------------------------------------------------------------
% 2. Calibrate parameter values for simulated series (i.e. true values
%----------------------------------------------------------------	
 @#if distrib == 1
    df_studt = 13;
@# endif
DELTA           = 0.025;
BETTA           = 0.984;
B               = 0.5;
ETAl            = 1; 
ETAc            = 2; 
THETA           = 3.48;
ALFA            = 0.667;
RHOA            = 0.979;
STDA            = 0.0072;

%----------------------------------------------------------------
% 3. Declare model equations
%----------------------------------------------------------------
model;    
% FOC for consumption
0 = -exp(la) +(exp(c)-B*exp(c(-1)))^(-ETAc) - BETTA*B*(exp(c(+1))-B*exp(c))^(-ETAc);
% Household's FOC for labor
0 = -THETA*(1-exp(n))^-ETAl + exp(la)*exp(w);
% FOC for capital
0 = -exp(la) + BETTA*exp(la(+1))*(exp(rk(+1)) + (1-DELTA));
% Firm's FOC for capital
0 = -exp(a)*(1-ALFA)*exp(k)^(-ALFA)*exp(n)^(ALFA) + exp(rk);
% Firm's FOC for labor
0 = -exp(a)*ALFA*exp(k)^(1-ALFA)*exp(n)^(ALFA-1) + exp(w);
% National income identity
0 = -exp(c) - exp(iv) + exp(y);
% Production function
0 = -exp(y) + exp(a)*exp(k)^(1-ALFA)*exp(n)^(ALFA);
% Law of motion for capital
0 = -exp(k(+1)) + (1-DELTA)*exp(k) + exp(iv);
% Law of motion for technology
0 = -log(exp(a)) + RHOA*log(exp(a(-1))) + STDA*u_a;
end;

%----------------------------------------------------------------
% 4. Specify variance of shock processes depending on distribution
%----------------------------------------------------------------
shocks;
@#if distrib == 1
    % Student-t distribution
    var u_a = df_studt/(df_studt-2);
@#else 
    % Gaussian distribution
    var u_a = 1;        
@#endif
end; 

%----------------------------------------------------------------
% 5. Specify steady-state (either steady_state_model or initval)
%----------------------------------------------------------------
% steady_state_model;
% errorMes = 0;
% A = 1; % The value of technology
% RK = 1/BETTA - (1-DELTA); % The rental rate of capital
% K_O_N = (RK/(A*(1-ALFA)))^(-1/ALFA); % Capital divided by labor
% % if K_O_N <= 0
% %     errorMes = 1;
% % end
% W = A*ALFA*(K_O_N)^(1-ALFA); % The wage level
% IV_O_N = DELTA*K_O_N; % Investment over labor
% Y_O_N = A*K_O_N^(1-ALFA); % Output over labor
% C_O_N = Y_O_N - IV_O_N; % Consumption over labor
% % if C_O_N <= 0
% %     errorMes = 1;
% % end
% %if ETAc == 1 && ETAl == 1
% N = (1-BETTA*B)*(C_O_N*(1-B))^-1*W/THETA/(1+(1-BETTA*B)*(C_O_N*(1-B))^-1*W/THETA);
% %else
% %     % No closed-form solution and we therefore use a fixed-point algorithm
% %     if errorMes == 0
% %         options = optimset('Display','off','TolX',1e-12,'TolFun',1e-12);
% %         N0 = 1/3;
% %         [N,~,exitflag] = fsolve(@findN,N0,options);
% %         if exitflag <= 0
% %             errorMes = 1;
% %         end
% %     else
% %         N = NaN;
% %     end
% % end
% C  = C_O_N*N;
% Y  = Y_O_N*N;
% IV = IV_O_N*N;
% K  = K_O_N*N;
% LA = (C-B*C)^(-ETAc)-BETTA*B*(C-B*C)^(-ETAc);
% k  = log(K);
% c  = log(C);
% a  = log(A);
% iv  = log(IV);
% y   = log(Y);
% la  = log(LA);
% n   = log(N);
% rk  = log(RK);
% w   = log(W);
% end; 

%----------------------------------------------------------------
% 6. Specify parameters for GMM estimation
%    Syntax: parameter, initial value, lower bound, upper bound, unifrom_pdf,,,,,tr
%    Note that uniform_pdf has no meaning, we simply abuse Dynare's syntax to enable 
%    parameter transformations (if selected as an option: opt.transpar = 1, default is opt.transpar = 0)
%    tr is parameter transformation type
%    0: no transformation needed
%    1: [a,b] -> [-1,1] -> [-inf,inf] by z/sqrt(1-z^2)
%    2: [0,inf] -> [-inf,inf] by b + ln(z-a);
%    with a: lower bound and b: upper bound:
%----------------------------------------------------------------
estimated_params;
   %parameter,     initial value, lower bound, upper bound, unifrom_pdf,,,,, tr
    DELTA,         0.025,         0,           1,           uniform_pdf,,,,, 1;
    BETTA,         0.984,         0,           1,           uniform_pdf,,,,, 1;
    B,             0.5,           0,           1,           uniform_pdf,,,,, 1;
    %ETAl,          1,             0,           10,          uniform_pdf,,,,, 1;
    ETAc,          2,             0,           10,          uniform_pdf,,,,, 1;
    ALFA,          0.667,         0,           1,           uniform_pdf,,,,, 1;
    RHOA,          0.979,         0,           1,           uniform_pdf,,,,, 1;
    STDA,          0.0072,        0,           1,           uniform_pdf,,,,, 1;
    %THETA,         3.48,          0,           10,           uniform_pdf,,,,, 1;
end;

%----------------------------------------------------------------
% 7. Computations
%----------------------------------------------------------------
steady; check;
stoch_simul(order=@{orderApp},pruning,noprint,nomoments,irf=0);
[resultsGMM,paramsGMM] = RunGMM(opt);