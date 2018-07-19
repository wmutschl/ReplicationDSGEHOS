% By Willi Mutschler, September 26, 2016. Email: willi@mutschler.eu
% Triplication Matrix as defined by
% Meijer (2005) - Matrix algebra for higher order moments. Linear Algebra and its Applications, 410,pp. 112–134
%
% Inputs:
%       p: size of vector
% Outputs:
%       TP: triplication matrix
%

function [TP,TPinv] = triplication(p,progress,sparseflag)
if nargin <2
    progress =0;
end
if nargin < 3
    sparseflag = 1;
end


if sparseflag
    TP = spalloc(p^3,p*(p+1)*(p+2)/6,p^2);
    if nargout >1
        TPinv = spalloc(p*(p+1)*(p+2)/6,p*(p+1)*(p+2)/6,p^3);
    end
else
    TP = zeros(p^3,p*(p+1)*(p+2)/6);
    if nargout >1
        TPinv = zeros(p*(p+1)*(p+2)/6,p*(p+1)*(p+2)/6);
    end
end

reverseStr = ''; counti=1;
np = p*(p+1)*(p+2)/6;


for k=1:p
    for j=k:p
        for i=j:p
            if progress && (rem(counti,50)== 0)
                msg = sprintf('    Triplication Matrix Processed %d/%d', counti, np); fprintf([reverseStr, msg]); reverseStr = repmat(sprintf('\b'), 1, length(msg));
            elseif progress && (counti==np)
                msg = sprintf('    Triplication Matrix Processed %d/%d\n', counti, np); fprintf([reverseStr, msg]); reverseStr = repmat(sprintf('\b'), 1, length(msg));
            end
            idx = uperm([i j k]);
            for r = 1:size(idx,1)                
                ii = idx(r,1); jj= idx(r,2); kk=idx(r,3);
                n = ii + (jj-1)*p + (kk-1)*p^2;                                    
                m = mue(p,i,j,k);
                TP(n,m)=1;
                if nargout >1
                    if i==j && j==k
                        TPinv(m,n)=1;
                    elseif i>j && j==k
                        TPinv(m,n)=1/3;
                    elseif i==j && j>k
                        TPinv(m,n)=1/3;
                    elseif i>j && j>k
                        TPinv(m,n)=1/6;                
                    end
                end
            end
            counti = counti+1;
            n=n+1;
        end
    end
end


function m = mue(p,i,j,k)
    m = i+(j-1)*p + 1/2*(k-1)*p^2 - 1/2*j*(j-1) + 1/6*k*(k-1)*(k-2) - 1/2*(k-1)^2*p;
    m = round(m);
end

end