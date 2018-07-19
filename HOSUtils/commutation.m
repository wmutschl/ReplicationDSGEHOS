% By Willi Mutschler, September 26, 2016. Email: willi@mutschler.eu
% Returns Magnus and Neudecker's commutation matrix of dimensions n by m
%
% Inputs: Dimension of original n times m matrix X
%
% Calls: vec for vectorizing a matrix
%
% Outputs: k: commutation matrix that solves k*vec(X)=vec(X')
%
% Original author: Thomas P Minka (tpminka@media.mit.edu), April 22, 2013
%

function k = commutation(n, m,sparseflag)

if nargin < 2
  m = n(2);
  n = n(1);
end

if nargin < 3
    sparseflag = 1;
end

if 0
  % first method
  i = 1:(n*m);
  a = reshape(i, n, m);
  j = vec(transpose(a));
  k = zeros(n*m,n*m);
  for r = i
    k(r, j(r)) = 1;
  end
else
  % second method
  if sparseflag
    k = reshape(kron(vec(speye(n)), speye(m)), n*m, n*m);
  else
    k = reshape(kron(vec(eye(n)), eye(m)), n*m, n*m);
  end
end


function V = vec(A)
    V = A(:); 
end

end