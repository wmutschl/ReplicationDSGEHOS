% By Willi Mutschler, September 26, 2016. Email: willi@mutschler.eu
% Duplication Matrix as defined by Magnus and Neudecker (2002), p.49
%
% Inputs:
%       n: size of vector
% Outputs:
%       d: Duplication matrix
%

function [DP,DPinv] = duplication(p,sparseflag)
if nargin < 2
    sparseflag = 1;
end
a = tril(ones(p));
i = find(a);
a(i) = 1:length(i);
a = a + transpose(tril(a,-1));

j = a(:);

m = p*(p+1)/2;
if sparseflag
    DP = spalloc(p*p,m,p^2);
else
    DP = zeros(p*p,m);
end
for r = 1:size(DP,1)
  DP(r, j(r)) = 1;
end

if nargout > 1
    DPinv = (transpose(DP)*DP)\transpose(DP);    
end

