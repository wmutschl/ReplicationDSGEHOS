% By Willi Mutschler, September 26, 2016. Email: willi@mutschler.eu
function y = mom_invGamma(v,k)
    if k > 0
        a =v/2; b=v/2;
        divisor = a-1;
        for ii=2:k
            divisor = divisor*(a-ii);
        end
        y = b^k/divisor;
    elseif k== 0
        y=1;
    else
        error('Something went wrong with the invGamma moments')
    end    
end %mom_inv_gamma end