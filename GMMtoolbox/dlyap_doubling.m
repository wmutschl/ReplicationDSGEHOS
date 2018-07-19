% Original author: Martin M. Andreasen
% This function solve the discrete time lyaponov system by a doubling
% algorithm. I.e. varx = A*varx*A' + B
% These codes are taken from some work by S. Schmitt-Grohe and M. Uribe.
function [varx,error_mes] = dlyap_doubling(A_input,B_input,varx_old)

% The doubling algorithm
error_mes = 0;
max_iter  = 500;
A         = A_input;
B         = B_input;
if isempty(varx_old) || size(varx_old,1) ~= size(A,1) 
    varx_old  = eye(size(A))*0.001;
end
difference= 0.1;
index1    = 1;
while difference > 1D-16 && index1 < max_iter;    
    varx       = A*varx_old*A'+B;
    difference = max(max(abs(varx-varx_old)));
    if difference > 1e-25
        B          = A*B*A'+B;
        A          = A*A;
        varx_old   = varx;
        index1     = index1 + 1;
    end    
end
if index1 == max_iter
    error_mes = 1;
end
