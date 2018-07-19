% By Willi Mutschler, September 26, 2016. Email: willi@mutschler.eu
function X1 = GetKronMom2(X0,E_A,E_B)
    % Get noncentral Kronecker Moments E[kron(A,B)] from central Kronecker Moments E[kron(A-E(A),B-E(B))]
    X1 = X0...
        + kron(E_A,E_B);
end