% By Willi Mutschler, September 26, 2016. Email: willi@mutschler.eu
function y = invtrans(params,lowerBounds,upperBounds,tr)
% INVTRANS(PARA)
%   transforms variables from model(restricted) 
%   to max(unrestricted)
%
% case 0, (-inf,inf) 
%   no adjustment needed
%   
% case 1, [a,b]
%   x/sqrt(1-x^2), x in [0,1]
%
% case 2, (a,inf)
%   ln(x), x in (0,inf)
%


y = params;

for i=1:length(params)
    a = lowerBounds(i);
    b = upperBounds(i);
    c = 1;
    if tr(i)==1;
        cx = 2*(params(i)-(a+b)/2)/(b-a);
        y(i) = (1/c)*cx/sqrt(1-cx^2);
    elseif tr(i)==2;
        y(i) = b + (1/c)*log(params(i)-a);
    end
end
