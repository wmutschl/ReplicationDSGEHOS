% By Willi Mutschler, September 26, 2016. Email: willi@mutschler.eu
% Input: 
% Z = E(kron(z,z)) with z = z1,...,n_z, note that Z is a vector of all second-order (cross-)product moments
% idx_i : index vector of variables in zi = z(idx_i)
% idx_j : index vecotr of variables in zj = z(idx_j)
% n_z: length of vector z
% Output:  X = E(kron(zi,zj)), (n_X times 1)
function X = GetSubMoms2(Z,idx_i,idx_j,n_z)
n_X = length(idx_i) + length(idx_j);
X=zeros(n_X,1);
counti=1; 
for i1=idx_i
    for i2=idx_j
        X(counti,1) = Z((i1-1)*n_z+i2);
        counti=counti+1;                    
    end        
end
end%GetSubMoms end