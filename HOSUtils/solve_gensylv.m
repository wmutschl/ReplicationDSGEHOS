% By Willi Mutschler, September 26, 2016. Email: willi@mutschler.eu
function [X,error_mes] = solve_gensylv(A,B,C,X_old,sylv_tol,sylv_max_iter,sylv_sparseflag,symmetric,sylv_method)
% Solves the equation X = A*X*B + C
if isempty(X_old) || size(X_old,1) ~= size(A,1)    
    if sylv_sparseflag        
        X_old = spalloc(size(A,2),size(B,1),size(A,2)*size(B,1));
    else
        X_old  = zeros(size(A,2),size(B,1));
    end
end
if symmetric==1
    C=0.5*(triu(C)+tril(C)');
    C=triu(C)+triu(C,1)';
end

error_mes = 0;
index1    = 0;
difference= 1+sylv_tol;

switch sylv_method
    case 'doubling'
        if isequal(A,transpose(B))
            equal_flag = 1;
        else
            equal_flag = 0;
        end        
        % Doubling algorithm based on Andreasens and Schmitt-Grohe and Uribe
        while difference > sylv_tol && index1 <= sylv_max_iter;
            X       = A*X_old*B+C;
            difference = max(max(abs(X-X_old)));
            if difference > 1e-25
                C          = A*C*B+C;
                A          = A*A;
                if equal_flag
                    B         =transpose(A);
                else
                    B          = B*B;
                end
                X_old      = X;
                index1     = index1 + 1;
            end
        end
        if index1 == sylv_max_iter && difference>sylv_tol
            warning(['Convergence not achieved in doubling algorithm of Sylvester equation after ' int2str(max_iter) ' iterations']);
            error_mes = 1;
        end        
    case 'fixed-point'
        % Fixed-point method based on Dynare's fastgensylv
        X = X_old;
        Z = A * X * B + C;

        while difference > sylv_tol && index1 <= sylv_max_iter;
            X = Z;    
            Z_old = Z;
            Z = A * X * B + C;
            difference = max(sum(abs(Z-Z_old)));
            index1 = index1 + 1;    
        end

        if index1 == sylv_max_iter && difference>sylv_tol
            warning(['fastgensylv:: Convergence not achieved in fixed point solution of Sylvester equation after ' int2str(maxit) ' iterations']);
            error_mes = 1;
        end
    case 'dlyap'        
        if isequal(A,transpose(B))
            disp('willi')
            X = dlyap(A,C);
        else
            X = dlyap(A,B,C);
        end    
end

if symmetric==1
    X=0.5*(triu(X)+tril(X)');
    X=triu(X)+triu(X,1)';
end
