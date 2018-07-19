% By Willi Mutschler, September 26, 2016. Email: willi@mutschler.eu
% prodmom.m		Date: 4/29/2006
% This Matlab program computes the product moment of X_{i_1}^{nu_1}X_{i_2}^{nu_2}...X_{i_m}^{nu_m},
% where X_{i_j} are elements from X ~ N(0_n,V).  
% V only needs to be positive semidefinite.
% V: variance-covariance matrix of X
% ii: vector of i_j
% nu: power of X_{i_j} 
% Reference: Triantafyllopoulos (2003) On the Central Moments of the Multidimensional
%            Gaussian Distribution, Mathematical Scientist
%            Kotz, Balakrishnan, and Johnson (2000), Continuous Multivariate
%            Distributions, Vol. 1, p.261
% Note that there is a typo in Eq.(46.25), there should be an extra rho in front 
% of the equation.
% Usage: prodmom(V,[i1 i2 ... ir],[nu1 nu2 ... nur])
% Example: To get E[X_2X_4^3X_7^2], use prodmom(V,[2 4 7],[1 3 2])
%
function y = prodmom(V,ii,nu);

if nargin<3
   nu = ones(size(ii));
end
s = sum(nu);
if s==0
   y = 1;
   return
end
if rem(s,2)==1
   y = 0;
   return
end
nuz = nu==0;
nu(nuz) = [];
ii(nuz) = [];
m = length(ii);
s2 = s/2;
%
%  Use univariate normal results
%
if m==1
   ii = find(nuz~=1);
   y = V(ii,ii)^s2*prod([1:2:s-1]);
   return
end
%
%  Use bivariate normal results when there are only two distinct indices
%
if m==2
   ii = find(nuz~=1);
   rho = V(ii(1),ii(2))/sqrt(V(ii(1),ii(1))*V(ii(2),ii(2)));
   y = V(ii(1),ii(1))^(nu(1)/2)*V(ii(2),ii(2))^(nu(2)/2)*bivmom(nu,rho);
   return  
end
%
%  Regular case
%
[nu,inu] = sort(nu,2,'descend');
ii = ii(inu);
V = V(ii,ii);          % Extract only the relevant part of V
x = zeros(1,m);
V = V./2;
nu2 = nu./2;
p = 2;
q = nu2*V*nu2';
y = 0;
for i=1:fix(prod(nu+1)/2)
    y = y+p*q^s2;
    for j=1:m
        if x(j)<nu(j)
           x(j) = x(j)+1;
           p = -round(p*(nu(j)+1-x(j))/x(j));
           q = q-2*(nu2-x)*V(:,j)-V(j,j);
           break
        else
           x(j) = 0;
           if rem(nu(j),2)==1
              p = -p;
           end
           q = q+2*nu(j)*(nu2-x)*V(:,j)-nu(j)^2*V(j,j);
        end
    end
end
y = y/prod([1:s2]);


%
% bivmom.m		Date: 1/11/2004
% This Matlab program computes the product moment of X_1^{p_1}X_2^{p_2},
% where X_1 and X_2 are standard bivariate normally distributed.
% n : dimension of X
% rho: correlation coefficient between X_1 and X_2
% Reference: Kotz, Balakrishnan, and Johnson (2000), Continuous Multivariate
%            Distributions, Vol. 1, p.261
% Note that there is a typo in Eq.(46.25), there should be an extra rho in front 
% of the equation.
% Usage: bivmom(p,rho)
%
function y = bivmom(p,rho);
s1 = p(1);
s2 = p(2);
rho2 = rho^2;
if rem(s1+s2,2)==1
   y = 0;
   return
end
r = fix(s1/2);
s = fix(s2/2);
y = 1;
c = 1;
odd = 2*rem(s1,2);
for j=1:min(r,s)
    c = 2*c*(r+1-j)*(s+1-j)*rho2/(j*(2*j-1+odd));
    y = y+c;
end
if odd
   y = y*rho;
end
y = prod([1:2:s1])*prod([1:2:s2])*y;
end %bivmom end

end % main function end
