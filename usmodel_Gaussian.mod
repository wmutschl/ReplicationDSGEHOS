%/////////////////////////////////////////////////////////////
%////    Model:   Smets and Wouters (2007)    		  	  ////
%////    Author:  Willi Mutschler                         ////
%////    Email:   willi@mutschler.eu                      ////
%////    Version: August 26, 2016      		   		      ////
%/////////////////////////////////////////////////////////////

%----------------------------------------------------------------
% 0. Specify options: User Settings (For other settings and defaults see updateOptions.m
%----------------------------------------------------------------
% Which distribution: 0:for Gaussian, 1:for Student_t
@#define distrib = 0

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
// copy van usmodel_hist_dsge_f19_7_71

var   labobs robs pinfobs dy dc dinve dw  ewma epinfma  zcapf rkf kf pkf    cf invef yf labf wf rrf mc zcap rk k pk    c inve y lab pinf w r a  b g qs  ms  spinf sw kpf kp ;    
 
varexo ea eb eg  eqs  em  epinf ew  ;  
 
varobs dy dc dinve labobs pinfobs dw robs;

parameters curvw cgy curvp constelab constepinf constebeta cmaw cmap calfa
czcap cbeta csadjcost ctou csigma chabb ccs cinvs cfc 
cindw cprobw cindp cprobp csigl clandaw 
crdpi crpi crdy cry crr 
crhoa crhoas crhob crhog crhols crhoqs crhoms crhopinf crhow  
ctrend 
conster cg cgamma clandap cbetabar cr cpie crk cw cikbar cik clk cky ciy ccy crkky cwhlc cwly ;
@#if distrib == 1
    parameters df_studt;
@# endif

%----------------------------------------------------------------
% 2. Calibrate parameter values for simulated series (i.e. true values
%----------------------------------------------------------------	
@#if distrib == 1
    df_studt = 5;
@# endif
// fixed parameters
ctou=.025;
clandaw=1.5;
cg=0.18;
curvp=10;
curvw=10;

// estimated parameters initialisation
calfa=.24;
cgamma=1.004;
cbeta=.9995;
csigma=1.5;
cpie=1.005;
cfc=1.5;
cgy=0.51;

csadjcost= 6.0144;
chabb=    0.6361;    
cprobw=   0.8087;
csigl=    1.9423;
cprobp=   0.6;
cindw=    0.3243;
cindp=    0.47;
czcap=    0.2696;
crpi=     1.488;
crr=      0.8762;
cry=      0.0593;
crdy=     0.2347;

crhoa=    0.9977;
crhob=    0.5799;
crhog=    0.9957;
crhols=   0.9928;
crhoqs=   0.7165;
crhoas=1; 
crhoms=0;
crhopinf=0;
crhow=0;
cmap = 0;
cmaw  = 0;

// derived from steady state
clandap=cfc;
cbetabar=cbeta*cgamma^(-csigma);
cr=cpie/(cbeta*cgamma^(-csigma));
crk=(cbeta^(-1))*(cgamma^csigma) - (1-ctou);
cw = (calfa^calfa*(1-calfa)^(1-calfa)/(clandap*crk^calfa))^(1/(1-calfa));
cikbar=(1-(1-ctou)/cgamma);
cik=(1-(1-ctou)/cgamma)*cgamma;
clk=((1-calfa)/calfa)*(crk/cw);
cky=cfc*(clk)^(calfa-1);
ciy=cik*cky;
ccy=1-cg-cik*cky;
crkky=crk*cky;
cwhlc=(1/clandaw)*(1-calfa)/calfa*crk*cky/ccy;
cwly=1-crk*cky;


ctrend=(cgamma-1)*100;
conster=(cr-1)*100;
constepinf=(cpie-1)*100;
constelab=0;


%----------------------------------------------------------------
% 3. Declare model equations
%----------------------------------------------------------------
model(linear); 

//usmodel_stst;

// flexible economy

	      0*(1-calfa)*a + 1*a =  calfa*rkf+(1-calfa)*(wf)  ;
	      zcapf =  (1/(czcap/(1-czcap)))* rkf  ;
	      rkf =  (wf)+labf-kf ;
	      kf =  kpf(-1)+zcapf ;
	      invef = (1/(1+cbetabar*cgamma))* (  invef(-1) + cbetabar*cgamma*invef(1)+(1/(cgamma^2*csadjcost))*pkf ) +qs ;
          pkf = -rrf-0*b+(1/((1-chabb/cgamma)/(csigma*(1+chabb/cgamma))))*b +(crk/(crk+(1-ctou)))*rkf(1) +  ((1-ctou)/(crk+(1-ctou)))*pkf(1) ;
	      cf = (chabb/cgamma)/(1+chabb/cgamma)*cf(-1) + (1/(1+chabb/cgamma))*cf(+1) +((csigma-1)*cwhlc/(csigma*(1+chabb/cgamma)))*(labf-labf(+1)) - (1-chabb/cgamma)/(csigma*(1+chabb/cgamma))*(rrf+0*b) + b ;
	      yf = ccy*cf+ciy*invef+g  +  crkky*zcapf ;
	      yf = cfc*( calfa*kf+(1-calfa)*labf +a );
	      wf = csigl*labf 	+(1/(1-chabb/cgamma))*cf - (chabb/cgamma)/(1-chabb/cgamma)*cf(-1) ;
	      kpf =  (1-cikbar)*kpf(-1)+(cikbar)*invef + (cikbar)*(cgamma^2*csadjcost)*qs ;

// sticky price - wage economy

	      mc =  calfa*rk+(1-calfa)*(w) - 1*a - 0*(1-calfa)*a ;
	      zcap =  (1/(czcap/(1-czcap)))* rk ;
	      rk =  w+lab-k ;
	      k =  kp(-1)+zcap ;
	      inve = (1/(1+cbetabar*cgamma))* (  inve(-1) + cbetabar*cgamma*inve(1)+(1/(cgamma^2*csadjcost))*pk ) +qs ;
          pk = -r+pinf(1)-0*b +(1/((1-chabb/cgamma)/(csigma*(1+chabb/cgamma))))*b + (crk/(crk+(1-ctou)))*rk(1) +  ((1-ctou)/(crk+(1-ctou)))*pk(1) ;
	      c = (chabb/cgamma)/(1+chabb/cgamma)*c(-1) + (1/(1+chabb/cgamma))*c(+1) +((csigma-1)*cwhlc/(csigma*(1+chabb/cgamma)))*(lab-lab(+1)) - (1-chabb/cgamma)/(csigma*(1+chabb/cgamma))*(r-pinf(+1) + 0*b) +b ;
	      y = ccy*c+ciy*inve+g  +  1*crkky*zcap ;
	      y = cfc*( calfa*k+(1-calfa)*lab +a );
	      pinf =  (1/(1+cbetabar*cgamma*cindp)) * ( cbetabar*cgamma*pinf(1) +cindp*pinf(-1) 
               +((1-cprobp)*(1-cbetabar*cgamma*cprobp)/cprobp)/((cfc-1)*curvp+1)*(mc)  )  + spinf ; 
	      w =  (1/(1+cbetabar*cgamma))*w(-1)
               +(cbetabar*cgamma/(1+cbetabar*cgamma))*w(1)
               +(cindw/(1+cbetabar*cgamma))*pinf(-1)
               -(1+cbetabar*cgamma*cindw)/(1+cbetabar*cgamma)*pinf
               +(cbetabar*cgamma)/(1+cbetabar*cgamma)*pinf(1)
               +(1-cprobw)*(1-cbetabar*cgamma*cprobw)/((1+cbetabar*cgamma)*cprobw)*(1/((clandaw-1)*curvw+1))*
               (csigl*lab + (1/(1-chabb/cgamma))*c - ((chabb/cgamma)/(1-chabb/cgamma))*c(-1) -w) 
               + 1*sw ;
	      r =  crpi*(1-crr)*pinf
               +cry*(1-crr)*(y-yf)     
               +crdy*(y-yf-y(-1)+yf(-1))
               +crr*r(-1)
               +ms  ;
	      a = crhoa*a(-1)  + ea;
	      b = crhob*b(-1) + eb;
	      g = crhog*(g(-1)) + eg + cgy*ea;
	      qs = crhoqs*qs(-1) + eqs;
	      ms = crhoms*ms(-1) + em;
	      spinf = crhopinf*spinf(-1) + epinfma - cmap*epinfma(-1);
	          epinfma=epinf;
	      sw = crhow*sw(-1) + ewma - cmaw*ewma(-1) ;
	          ewma=ew; 
	      kp =  (1-cikbar)*kp(-1)+cikbar*inve + cikbar*cgamma^2*csadjcost*qs ;

// measurment equations

dy=y-y(-1)+ctrend;
dc=c-c(-1)+ctrend;
dinve=inve-inve(-1)+ctrend;
dw=w-w(-1)+ctrend;
pinfobs = 1*(pinf) + constepinf;
robs =    1*(r) + conster;
labobs = lab + constelab;

end; 

%----------------------------------------------------------------
% 4. Specify variance of shock processes depending on distribution
%----------------------------------------------------------------
shocks;
@#if distrib == 1
    % Student-t distribution
    var ea = df_studt/(df_studt-2)*0.4618^2;
    var eb = df_studt/(df_studt-2)*1.8513^2;
    var eg = df_studt/(df_studt-2)*0.6090^2;
    var eqs = df_studt/(df_studt-2)*0.6017^2;
    var em = df_studt/(df_studt-2)*0.2397^2;
    var epinf = df_studt/(df_studt-2)*0.1455^2;
    var ew = df_studt/(df_studt-2)*0.2089^2;
@#else
    % Gaussian distribution
    var ea = 0.4618^2;
    var eb = 1.8513^2;
    var eg = 0.6090^2;
    var eqs = 0.6017^2;
    var em = 0.2397^2;
    var epinf = 0.1455^2;
    var ew = 0.2089^2;
@#endif
end;

%----------------------------------------------------------------
% 7. Computations
%----------------------------------------------------------------
steady; check;
stoch_simul(order=1,pruning,noprint,nomoments,irf=0);
DispHOS(opt);