% By Willi Mutschler, September 26, 2016. Email: willi@mutschler.eu
function [ys, info] = RBCmodel_steadystate(ys, exo)
global M_

% read out parameters to access them with their name
DELTA = M_.params(1);
BETTA = M_.params(2);
B = M_.params(3); 
ETAl = M_.params(4);
ETAc = M_.params(5);
THETA = M_.params(6);
ALFA = M_.params(7);
RHOA = M_.params(8);
STDA = M_.params(9);
    info = 0;
    errorMes=0;
    A=1;
    RK = 1/BETTA - (1-DELTA);
    K_O_N = (RK/(A*(1-ALFA)))^(-1/ALFA);
    if K_O_N <= 0
        errorMes = 1;
    end
    W = A*ALFA*(K_O_N)^(1-ALFA);
    IV_O_N = DELTA*K_O_N;
    Y_O_N = A*K_O_N^(1-ALFA);
    C_O_N = Y_O_N - IV_O_N;
    if C_O_N <= 0
        errorMes = 1;
    end
    if ETAc == 1 && ETAl == 1
        N = (1-BETTA*B)*(C_O_N*(1-B))^-1*W/THETA/(1+(1-BETTA*B)*(C_O_N*(1-B))^-1*W/THETA);
    else
    % No closed-form solution and we therefore use a fixed-point algorithm
    if errorMes == 0
        options = optimset('Display','off','TolX',1e-12,'TolFun',1e-12);
        N0 = 1/3;
        [N,~,exitflag] = fsolve(@findN,N0,options);
        if exitflag <= 0
            errorMes = 1;
        end
    else
        N = NaN;
    end
    end
    
    C=C_O_N*N;
    Y=Y_O_N*N;
    IV=IV_O_N*N;
    K=K_O_N*N;
    LA = (C-B*C)^(-ETAc)-BETTA*B*(C-B*C)^(-ETAc);
    ys(1)=log(K);
    ys(2)=log(C);
    ys(3)=log(A);
    ys(4)=log(IV);
    ys(5)=log(Y);
    ys(6)=log(LA);
    ys(7)=log(N);
    ys(8)=log(RK);
    ys(9)=log(W);
    % Auxiliary equations
    check_=0;
    
function error = findN(N)
error = THETA*(1-N)^(-ETAl)*N^ETAc - (1-BETTA*B)*(C_O_N*(1-B))^(-ETAc)*W;
end
end
