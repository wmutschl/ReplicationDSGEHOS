% By Willi Mutschler, September 26, 2016. Email: willi@mutschler.eu
function y = trans(params,lowerBounds,upperBounds,tr)
% TRANS(PARA)
%   transforms variables from max(unrestricted) 
%   to model(restricted)
%
% case 0, (-inf,inf) 
%   no adjustment needed
%   
% case 1, [a,b]
%   x/sqrt(1+x^2), x unbounded
%   0.5*[ (a+b) + (b-a)*c*x/sqrt(1+c^2*x^2) ]
%
% case 2, (a,inf)
%   exp(x), x unbounded
%   a + exp(c*(x-b))
%
y = params;
for i=1:length(params)
    a = lowerBounds(i);
    b = upperBounds(i);
    c = 1;

    if tr(i)==1;
        y(i) = (a+b)/2 + 0.5*(b-a)*c*params(i)/sqrt(1+c^2*params(i)^2);
    elseif tr(i)==2;
        y(i) = a + exp(c*(params(i)-b));
    end
end
