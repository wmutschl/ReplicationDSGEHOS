% By Willi Mutschler, September 26, 2016. Email: willi@mutschler.eu
% Permutation matrix for fourth-order cumulant
function P = PermutMat(nximin,sparseflag)
if nargin < 2
    sparseflag = 1;
end
if sparseflag
    U = spalloc(nximin^3,nximin^3,nximin^3);
else
    U = spalloc(nximin^3,nximin^3);
end
for i=1:nximin^2
    for k=1:nximin
        U((i-1)*nximin+k,(k-1)*nximin^2+i) = 1;        
    end
end
if sparseflag
    P = kron(speye(nximin),U); 
else
    P = kron(eye(nximin),U); 
end
end