% By Willi Mutschler, September 26, 2016. Email: willi@mutschler.eu
function [DP5,DP5inv] = Q5_plication(p,progress)
if nargin <2
    progress =0;
end
reverseStr = ''; counti=1;
np = p*(p+1)*(p+2)*(p+3)*(p+4)/(1*2*3*4*5);
DP5 = spalloc(p^5,p*(p+1)*(p+2)*(p+3)*(p+4)/(1*2*3*4*5),p^5);

for i1=1:p
    for i2=i1:p
        for i3=i2:p
            for i4=i3:p
                for i5=i4:p
                    if progress && (rem(counti,100)== 0)
                        msg = sprintf('    Q5-plication Matrix Processed %d/%d', counti, np); fprintf([reverseStr, msg]); reverseStr = repmat(sprintf('\b'), 1, length(msg));
                    elseif progress && (counti==np)
                        msg = sprintf('    Q5-plication Matrix Processed %d/%d\n', counti, np); fprintf([reverseStr, msg]); reverseStr = repmat(sprintf('\b'), 1, length(msg));
                    end
                    idx = uperm([i5 i4 i3 i2 i1]);
                    for r = 1:size(idx,1)
                        ii1 = idx(r,1); ii2= idx(r,2); ii3=idx(r,3); ii4=idx(r,4); ii5=idx(r,5);
                        n = ii1 + (ii2-1)*p + (ii3-1)*p^2 + (ii4-1)*p^3  + (ii5-1)*p^4;
                        m = mue(p,i5,i4,i3,i2,i1);
                        DP5(n,m)=1;
                    end
                    counti = counti+1;
                end
            end
        end
    end
end

DP5inv = (transpose(DP5)*DP5)\transpose(DP5);

function m = mue(p,i1,i2,i3,i4,i5)
     m = binom_coef(p,5,1) - binom_coef(p,1,i1+1) - binom_coef(p,2,i2+1) - binom_coef(p,3,i3+1) - binom_coef(p,4,i4+1) - binom_coef(p,5,i5+1);
     m = round(m);
end

function N = binom_coef(p,q,i)
    t = q; r =p+q-i;
    if t==0
        N=1;
    else
        N=1;
        for h = 0:(t-1)
            N = N*(r-h);
        end
        N=N/factorial(t);
    end
end
end