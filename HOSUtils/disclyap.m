% By Willi Mutschler, September 26, 2016. Email: willi@mutschler.eu
function [X,error_mes] = disclyap(A,C,X_old,sylv_tol,sylv_max_iter,sylv_sparseflag,symmetric,setting)
nA = size(A,1);
error_mes = 0;
index1    = 0;
difference= 1+sylv_tol;

if symmetric==1
    C=0.5*(triu(C)+tril(C)');
    C=triu(C)+triu(C,1)';
end

switch setting
    case 2
        % Solves the equation X = A*X*A' + C
        if isempty(X_old) || size(X_old,1) ~= nA || size(X_old,2) ~= nA            
            if sylv_sparseflag        
                X_old = spalloc(nA,nA,nA^2);
            else
                X_old  = zeros(nA,nA);
            end
        end
        
        while difference > sylv_tol && index1 <= sylv_max_iter;
            X       = A*X_old*transpose(A)+C;
            difference = max(max(abs(X-X_old)));
            if difference > 1e-25
                C          = A*C*transpose(A)+C;
                A          = A*A;
                X_old      = X;
                index1     = index1 + 1;
            end
        end
    case 3
        % Solves the equation X=kron(A,A)*X*A' + C
        if isempty(X_old) || size(X_old,1) ~= nA^2 || size(X_old,2) ~= nA            
            if sylv_sparseflag        
                X_old = spalloc(nA^2,nA,nA^3);
            else
                X_old  = zeros(nA^2,nA);
            end
        end
        while difference > sylv_tol && index1 <= sylv_max_iter;
            A_kron=KronProd({A,A,1},[1,2], [],1);
            X       = A_kron*X_old*transpose(A)+C;
            difference = max(max(abs(X-X_old)));
            if difference > 1e-25
                C          = A_kron*C*transpose(A)+C;
                A          = A*A;
                X_old      = X;
                index1     = index1 + 1;
            end
        end
    case 4
        % Solves the equation X=kron(A,A)*X*kron(A,A)' + C
        if isempty(X_old) || size(X_old,1) ~= nA^2 || size(X_old,2) ~= nA^2            
            if sylv_sparseflag        
                X_old = spalloc(nA^2,nA^2,nA^4);
            else
                X_old  = zeros(nA^2,nA^2);
            end
        end
        while difference > sylv_tol && index1 <= sylv_max_iter;
            A_kron=KronProd({A,A,1},[1,2], [],1);
            X       = A_kron*X_old*transpose(A_kron)+C;
            difference = max(max(abs(X-X_old)));
            if difference > 1e-25
                C          = A_kron*C*transpose(A_kron)+C;
                A          = A*A;
                X_old      = X;
                index1     = index1 + 1;
            end
        end
end

if index1 == sylv_max_iter && difference>sylv_tol
    warning(['Convergence not achieved in doubling algorithm of Sylvester equation after ' int2str(max_iter) ' iterations']);
    error_mes = 1;
end        

if symmetric==1
    X=0.5*(triu(X)+tril(X)');
    X=triu(X)+triu(X,1)';
end
