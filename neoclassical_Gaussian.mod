%/////////////////////////////////////////////////////////////
%////    Model: Neolassical growth model                  ////
%////    Author: Willi Mutschler                         ////
%////    Email:  willi@mutschler.eu                     ////
%////    Version: August 26, 2016      		   		      ////
%/////////////////////////////////////////////////////////////
%----------------------------------------------------------------
% 0. Specify options: User Settings (For other settings and defaults see updateOptions.m
%----------------------------------------------------------------
% Order of approximation
@#define orderApp = 3
% Which distribution: 0:for Gaussian, 1:for Student_t
@#define distrib = 0
% Do GMM estimation: 0 no or 1 yes
@#define GMMestimation   = 0

@#if distrib == 1
    opt.distrib = 'Student_t'; %multivariate student-t distribution ('Gaussian' is default)
@#endif

@#if GMMestimation == 1    
    opt.numSim               = 250;     % Number of observations of simulated data
    opt.autoLagsIdx          = [1]    
    opt.HOSThirdOrder        = 1;         % Use third-order statistics,  0: do not include or 1: include
    opt.HOSFourthOrder       = 0;         % User fourth-order statistics, 0: do not include or 1: include
    opt.AUXTEST              = 0;         % AUXILIARY TEST IF EMPIRICAL MOMENTS FIT TO THEORETICAL MOMENTS    
    opt.seed_nr              = 10;
    opt.transpar             = 0;
@#else
    opt.HOSThirdOrder        = 1;         % Use third-order statistics,  0: do not include or 1: include
    opt.HOSFourthOrder       = 1;         % User fourth-order statistics, 0: do not include or 1: include
    opt.MCruns               = 1000;      % Monte Carlo Runs for comparison
    opt.numSim               = 10000;     % Number of observations of simulated data
    opt.burnin               = 1000;      % Initial burnin periods, will be discarded from simulation
    opt.antithetic           = 1;         % 1: Use antithetic shocks and quadratic resampling to reduce Monte Carlo variation
    opt.simdata              = 1;         % compare theoretical statistics with simulated data, 0: no or 1: yes
@#endif

%----------------------------------------------------------------
% 1. Declare variables and parameters
%----------------------------------------------------------------
var k a c;
varexo u_a;
varobs c;
parameters ALFA RHO SIGMA_A BETTA DELTA GAM;
@#if distrib == 1
    parameters df_studt; % Please always use df_studt for degrees of freedom
@# endif

%----------------------------------------------------------------
% 2. Calibrate parameter values for simulated series (i.e. true values
%----------------------------------------------------------------	
 @#if distrib == 1
    df_studt = 13;
@# endif
ALFA = 0.3; 
RHO = 0; 
BETTA = 0.95; 
DELTA = 1; 
GAM = 2;
SIGMA_A = 1;

%----------------------------------------------------------------
% 3. Declare model equations
%----------------------------------------------------------------
model;    
-exp(k) + (1-DELTA)*exp(k(-1))+exp(a(-1)+ALFA*k(-1))-exp(c)=0;
-exp(-GAM*c) + BETTA * exp(-GAM*c(+1))*(ALFA*exp(a+(ALFA-1)*k)+1-DELTA)=0;
-a + RHO * a(-1) + SIGMA_A*u_a=0;        
end;

%----------------------------------------------------------------
% 4. Specify variance of shock processes depending on distribution
%----------------------------------------------------------------
shocks;
@#if distrib == 1
    % Student-t distribution
    var u_a = df_studt/(df_studt-2)*1^2;
@#else 
    % Gaussian distribution
    var u_a = 1^2;        
@#endif
end; 

%----------------------------------------------------------------
% 5. Specify steady-state (either steady_state_model or initval)
%----------------------------------------------------------------
steady_state_model;	        
a = log(1);
k = log(((1/BETTA+DELTA-1)/ALFA)^(1/(ALFA-1)));
c = log(exp(a+ALFA*k)-DELTA*exp(k));
end; 
    
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
%parameter, initial value, lower bound, upper bound, unifrom_pdf,,,,, tr
    ALFA,    0.1,           0.05,        0.9,         uniform_pdf,,,,, 1; 
    RHO,     0.5,           1e-5,        0.99999,     uniform_pdf,,,,, 1; 
    %BETTA,   0.9,           0.8,         0.99999,     uniform_pdf,,,,, 1; 
    %DELTA,   0.2,          0.001,        1,         uniform_pdf,,,,, 1; 
    GAM,     2,             0.5,         3,           uniform_pdf,,,,, 0; 
    SIGMA_A, 0.5,           1e-5,        5,           uniform_pdf,,,,, 2; 
end;

%----------------------------------------------------------------
% 7. Computations
%----------------------------------------------------------------
steady; check;
stoch_simul(order=@{orderApp},pruning,noprint,nomoments,irf=0);
DispHOS(opt);
@#if GMMestimation
    [resultsGMM,paramsGMM] = RunGMM(opt);
@#endif