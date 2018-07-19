%/////////////////////////////////////////////////////////////
%////    Model:   An and Schorfheide (2007)    		  	  ////
%////    Author:  Willi Mutschler                         ////
%////    Email:   willi@mutschler.eu                      ////
%////    Version: August 26, 2016      		   		      ////
%/////////////////////////////////////////////////////////////

%----------------------------------------------------------------
% 0. Specify options: User Settings (For other settings and defaults see updateOptions.m
%----------------------------------------------------------------
% Order of approximation
@#define orderApp = 2
% Which distribution: 0:for Gaussian, 1:for Student_t
@#define distrib = 1

opt.HOSThirdOrder        = 1;         % Use third-order statistics,  0: do not include or 1: include
opt.HOSFourthOrder       = 1;         % User fourth-order statistics, 0: do not include or 1: include
opt.simdata              = 1;         % compare theoretical statistics with simulated data, 0: no or 1: yes
opt.MCruns               = 1000;      % Monte Carlo Runs for comparison
opt.numSim               = 10000;     % Number of observations of simulated data
opt.burnin               = 1000;      % Initial burnin periods, will be discarded from simulation
opt.antithetic           = 1;         % 1: Use antithetic shocks and quadratic resampling to reduce Monte Carlo variation

@#if distrib == 1
    opt.distrib = 'Student_t'; %multivariate student-t distribution ('Gaussian' is default)
@#endif

%----------------------------------------------------------------
% 1. Declare variables and parameters
%----------------------------------------------------------------
var c dy p y R g z YGR INFL INT;
varexo e_z e_g e_r;

varobs YGR INFL INT;
parameters tau nu kap cyst psi1 psi2 rhor rhog rhoz rrst pist gamst sig_r sig_g sig_z;
@#if distrib == 1
    parameters df_studt; % Please always use df_studt for degrees of freedom
@# endif

%----------------------------------------------------------------
% 2. Calibrate parameter values for simulated series (i.e. true values
%----------------------------------------------------------------	
@#if distrib == 1
    df_studt = 9;
@# endif
tau   = 2.0000; 
nu    = 0.1000;
kap   = 0.3300;
cyst  = 0.8500;
psi1  = 1.5000;
psi2  = 0.1250;
rhor  = 0.7500;
rhog  = 0.9500;
rhoz  = 0.9000;
rrst  = 1.0000;
pist  = 3.2000;
gamst = 0.5500;
sig_r = 0.002;
sig_g = 0.006;
sig_z = 0.003;
    
%----------------------------------------------------------------
% 3. Declare model equations
%----------------------------------------------------------------
model;
% Auxiliary parameters and variables	
#pist2 = exp(pist/400);
#rrst2 = exp(rrst/400);
#bet   = 1/rrst2;
#phi   = tau*(1-nu)/nu/kap/pist2^2;
#gst   = 1/cyst;
#cst   = (1-nu)^(1/tau);
#yst   = cst*gst;

% Euler equation, eq. (21)
1 = exp(-tau*c(+1)+tau*c+R-z(+1)-p(+1));
% Phillips curve, eq. (22)
(1-nu)/nu/phi/(pist2^2)*(exp(tau*c)-1) = (exp(p)-1)*((1-1/2/nu)*exp(p)+1/2/nu) - bet*(exp(p(+1))-1)*exp(-tau*c(+1)+tau*c+dy(+1)+p(+1));
% Equilibrium condition, eq. (23)
exp(c-y) = exp(-g) - phi*pist2^2*gst/2*(exp(p)-1)^2;
% Taylor Rule, eq. (24)    
R = rhor*R(-1) + (1-rhor)*psi1*p + (1-rhor)*psi2*(y-g) + e_r;    
% Fiscal rule, eq. (25)	
g = rhog*g(-1) + e_g;
% Evolution of technology, eq (26)
z = rhoz*z(-1) + e_z;
% Auxiliary equation for output Growth
dy = y - y(-1);
% Measurement equations, eq. (38)
YGR = gamst+100*(dy+z);
INFL = pist+400*p;
INT = pist+rrst+4*gamst+400*R;    
end;

%----------------------------------------------------------------
% 4. Specify variance of shock processes depending on distribution
%----------------------------------------------------------------
shocks;
@#if distrib == 1
    % Student-t distribution
    var e_r = df_studt/(df_studt-2)*sig_r^2;
    var e_g = df_studt/(df_studt-2)*sig_g^2;
    var e_z = df_studt/(df_studt-2)*sig_z^2;
@#else 
    % Gaussian distribution
    var e_r = sig_r^2;
    var e_g = sig_g^2;
    var e_z = sig_z^2;
@#endif
end;

%----------------------------------------------------------------
% 5. Specify steady-state (either steady_state_model or initval)
%----------------------------------------------------------------
steady_state_model;     
y    = 0;
R    = 0;
g    = 0;
z    = 0;
c    = 0;
dy   = 0;
p    = 0;
YGR  = gamst;
INFL = pist;
INT  = pist + rrst + 4*gamst;
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
%   tau,     2,             1e-5,        10,          uniform_pdf,,,,, 2; 
   nu,      0.1,           1e-5,        0.99999,     uniform_pdf,,,,, 1; 
   kap,     0.3,           1e-5,        10,          uniform_pdf,,,,, 2; 
   cyst,    0.85,          1e-5,        0.99999,     uniform_pdf,,,,, 1; 
   psi1,    1,             0.1,         10,          uniform_pdf,,,,, 2; 
%   psi2,    0.5,           1e-5,        10,          uniform_pdf,,,,, 2; 
   rhor,    0.5,           1e-5,        0.99999,     uniform_pdf,,,,, 1; 
   rhog,    0.8,           1e-5,        0.99999,     uniform_pdf,,,,, 1; 
   rhoz,    0.66,          1e-5,        0.99999,     uniform_pdf,,,,, 1; 
%   rrst,    0.8,           1e-5,        10,          uniform_pdf,,,,, 2; 
%   pist,    4,             1e-5,        20,          uniform_pdf,,,,, 2; 
%   gamst,   0.4,           -5,          5,           uniform_pdf,,,,, 0; 
%   sig_r,   0.003, 1e-8, 5,       uniform_pdf,,,,,2; 
%   sig_g,   0.004, 1e-8, 5,       uniform_pdf,,,,,2; 
%   sig_z,   0.004, 1e-8, 5,       uniform_pdf,,,,,2; 
end;

%----------------------------------------------------------------
% 7. Computations
%----------------------------------------------------------------
steady; check;
stoch_simul(order=@{orderApp},pruning,noprint,nomoments,irf=0);
DispHOS(opt);