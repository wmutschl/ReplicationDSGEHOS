% By Willi Mutschler, September 26, 2016. Email: willi@mutschler.eu
% Input: 
% Z = E(kron(z,kron(z,z))) with z = z1,...,n_z, note that Z is a vector of all third-order (cross-)product moments
% idx_i : index vector of variables in zi = z(idx_i)
% idx_j : index vecotr of variables in zj = z(idx_j)
% idx_k : index vecotr of variables in zk = z(idx_k)
% n_z: length of vector z
% Output:  X = E(kron(zi,kron(zj,zk))), (n_X times 1)
function X = GetSubMoms3(Z,idx_i,idx_j,idx_k,n_z)
n_X = length(idx_i) + length(idx_j) + length(idx_k);
X=zeros(n_X,1);
counti=1; 
for i1=idx_i
    for i2=idx_j
        for i3=idx_k
            X(counti,1) = Z((i1-1)*n_z^2+(i2-1)*n_z+i3);
            counti=counti+1;                    
        end
    end        
end
end%GetSubMoms end