% By Willi Mutschler, September 26, 2016. Email: willi@mutschler.eu
% Quadruplication Matrix as defined by
% Meijer (2005) - Matrix algebra for higher order moments. Linear Algebra and its Applications, 410,pp. 112–134
%
% Inputs:
%       p: size of vector
% Outputs:
%       QP: quadruplication matrix
%       QPinv: Moore-Penrose inverse of QP
%
function [QP,QPinv] = quadruplication(p,progress,sparseflag)

if nargin <2
    progress =0;
end
if nargin < 3
    sparseflag = 1;
end
reverseStr = ''; counti=1;
np = p*(p+1)*(p+2)*(p+3)/24;

if sparseflag
    QP = spalloc(p^4,p*(p+1)*(p+2)*(p+3)/24,p^4);
else
    QP = zeros(p^4,p*(p+1)*(p+2)*(p+3)/24);
end
if nargout > 1
    if sparseflag
        QPinv = spalloc(p*(p+1)*(p+2)*(p+3)/24,p*(p+1)*(p+2)*(p+3)/24,p^4);
    else
        QPinv = zeros(p*(p+1)*(p+2)*(p+3)/24,p*(p+1)*(p+2)*(p+3)/24);
    end
end

for l=1:p
    for k=l:p
        for j=k:p
            for i=j:p                             
                if progress && (rem(counti,100)== 0)
                    msg = sprintf('    Quadruplication Matrix Processed %d/%d', counti, np); fprintf([reverseStr, msg]); reverseStr = repmat(sprintf('\b'), 1, length(msg));
                elseif progress && (counti==np)
                    msg = sprintf('    Quadruplication Matrix Processed %d/%d\n', counti, np); fprintf([reverseStr, msg]); reverseStr = repmat(sprintf('\b'), 1, length(msg));
                end
                idx = uperm([i j k l]);
                for r = 1:size(idx,1)
                    ii = idx(r,1); jj= idx(r,2); kk=idx(r,3); ll=idx(r,4);
                    n = ii + (jj-1)*p + (kk-1)*p^2 + (ll-1)*p^3;                    
                    m = mue(p,i,j,k,l);
                    QP(n,m)=1;
                    if nargout > 1
                        if i==j && j==k && k==l
                            QPinv(m,n)=1;
                        elseif i==j && j==k && k>l
                            QPinv(m,n)=1/4;
                        elseif i>j && j==k && k==l
                            QPinv(m,n)=1/4;
                        elseif i==j && j>k && k==l
                            QPinv(m,n) = 1/6;
                        elseif i>j && j>k && k==l
                            QPinv(m,n) = 1/12;
                        elseif i>j && j==k && k>l
                            QPinv(m,n) = 1/12;
                        elseif i==j && j>k && k>l
                            QPinv(m,n) = 1/12;
                        elseif i>j && j>k && k>l
                            QPinv(m,n) = 1/24;                    
                        end
                    end
                end
                counti = counti+1;
            end
        end
    end
end
%QPinv = (transpose(QP)*QP)\transpose(QP);

function m = mue(p,i,j,k,l)
     m = i + (j-1)*p + 1/2*(k-1)*p^2 + 1/6*(l-1)*p^3 - 1/2*j*(j-1) + 1/6*k*(k-1)*(k-2) - 1/24*l*(l-1)*(l-2)*(l-3) - 1/2*(k-1)^2*p + 1/6*(l-1)^3*p - 1/4*(l-1)*(l-2)*p^2 - 1/4*l*(l-1)*p + 1/6*(l-1)*p;
     m = round(m);
end


end