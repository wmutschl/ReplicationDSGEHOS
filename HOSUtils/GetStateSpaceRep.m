% By Willi Mutschler, September 26, 2016. Email: willi@mutschler.eu
function [A,B,C,D,c,d,ybar,Fxi,Fz,nx,ny,nu,nximin,id1_xf,id2_xs,id3_xf_xf] = GetStateSpaceRep(setup)
global M_ oo_
sigm = 1; %the perturbation parameter
Sigma = M_.Sigma_e;        
%eta = chol(Sigma)./sigm; % matrix that captures variance and correlation
%sigeta = sigm*eta;
indx_states = setup.indx_states;
indx_obs = oo_.dr.inv_order_var(setup.indx_obs); %note that solution matrices are stored in different ordering than variable declaration!
orderApp = setup.orderApp;
distrib = setup.distrib;
sparseflag = setup.sparseflag;
ybar = oo_.dr.ys(setup.indx_obs); %dr.ys is stored in same ordering as variable declaration

hx = oo_.dr.ghx(indx_states,:);     gx = oo_.dr.ghx(indx_obs,:);
hu = oo_.dr.ghu(indx_states,:);     gu = oo_.dr.ghu(indx_obs,:);
if orderApp >1
    Hxx = oo_.dr.ghxx(indx_states,:);   Gxx = oo_.dr.ghxx(indx_obs,:);
    Hxu = oo_.dr.ghxu(indx_states,:);   Gxu = oo_.dr.ghxu(indx_obs,:);
    Huu = oo_.dr.ghuu(indx_states,:);   Guu = oo_.dr.ghuu(indx_obs,:);
    hss = oo_.dr.ghs2(indx_states,:);   gss = oo_.dr.ghs2(indx_obs,:);
end
if orderApp >2
    Hxxx = oo_.dr.ghxxx(indx_states,:); Gxxx = oo_.dr.ghxxx(indx_obs,:);
    Hxxu = oo_.dr.ghxxu(indx_states,:); Gxxu = oo_.dr.ghxxu(indx_obs,:);
    Hxuu = oo_.dr.ghxuu(indx_states,:); Gxuu = oo_.dr.ghxuu(indx_obs,:);
    Huuu = oo_.dr.ghuuu(indx_states,:); Guuu = oo_.dr.ghuuu(indx_obs,:);
    Hxss = oo_.dr.ghxss(indx_states,:); Gxss = oo_.dr.ghxss(indx_obs,:);
    Huss = oo_.dr.ghuss(indx_states,:); Guss = oo_.dr.ghuss(indx_obs,:);
    Hsss = zeros(length(indx_states),1);Gsss = zeros(length(indx_obs),1); %Hsss and Gsss will be non-zero if the third moment of the shocks is non-zero, Dynare++ solves only for Gaussian shocks. 
end
ny = size(gx,1);
nx = size(hx,2);
nu = size(Sigma,2);

id1_xf       = (1:nx); 
id2_xs       = id1_xf(end)    + (1:nx); 
id3_xf_xf    = id2_xs(end)    + (1:nx^2);
id4_xrd      = id3_xf_xf(end) +(1:nx);
id5_xf_xs    = id4_xrd(end)   + (1:nx^2);
id6_xf_xf_xf = id5_xf_xs(end) + (1:nx^3);

id1_u        = 1:nu;
id2_u_u      = id1_u(end)       + (1:nu^2);    
id3_xf_u     = id2_u_u(end)     + (1:nx*nu);
id4_u_xf     = id3_xf_u(end)    + (1:nx*nu);
id5_xs_u     = id4_u_xf(end)    + (1:nx*nu);    
id6_u_xs     = id5_xs_u(end)    + (1:nx*nu);    
id7_xf_xf_u  = id6_u_xs(end)    + (1:nx^2*nu);
id8_xf_u_xf  = id7_xf_xf_u(end) + (1:nx^2*nu);
id9_u_xf_xf  = id8_xf_u_xf(end) + (1:nx^2*nu);    
id10_xf_u_u  = id9_u_xf_xf(end) + (1:nx*nu^2);
id11_u_xf_u  = id10_xf_u_u(end) + (1:nx*nu^2);
id12_u_u_xf  = id11_u_xf_u(end) + (1:nx*nu^2);
id13_u_u_u   = id12_u_u_xf(end) + (1:nu^3);     


switch distrib
    case 'Gaussian'
        M2u = vec(Sigma);
        M3u = zeros(nu^3,1);
    case 'Student_t'
        M2u = vec(Sigma);
        M3u = zeros(nu^3,1);
end

Fxi = GetFxi(nu,nx,orderApp,sparseflag);

if orderApp == 1
    % First-order state-space representation for all endogenous variables y given states xf and shocks u with E(u_t u_t') = sigma^2*eta*eta' = Sigma
    %   The system reads
    %     xf(t+1) = hx*xf(t) + hu*u(t+1) 
    %     y_(t+1) = ybar + gx*xf(t) + gu*u(t+1)
    % Set ABCD representation:
    % z(t) = c + A*z(t-1) + B*xi(t+1)
    % y(t) = ybar + d + C*z(t-1) + D*xi(t+1)
    % z(t)=xf(t) xi(t+1)=u(t+1)
    nz = nx; nxi = nu; nximin = nu;
    A = zeros(nz,nz);
    B = zeros(nz,nxi);
    C = zeros(ny,nz);
    D = zeros(ny,nxi);
    c = zeros(nz,1);
    d = zeros(ny,1);
    
    A(id1_xf,id1_xf) = hx; 
    B(id1_xf,id1_u) = hu; 
    C(1:ny,id1_xf) = gx;
    D(1:ny,id1_u) = gu;
    
    Fz = [];
elseif orderApp == 2
    % Second-order state-space representation for all endogenous variables y given first-order effects states xf, second-order effects states xs and shocks u with E(u_t u_t') = sigma^2*eta*eta' = Sigma
    %   The system reads
    %     xf(t+1) = hx*xf(t) + hu*u(t+1)
    %     xs(t+1) = hx*xs(t) + 1/2*Hxx*kron(xf(t),xf(t)) + 1/2*Huu*kron(u(t+1),u(t+1)) + Hxu*kron(xf(t),u(t+1)) + 1/2*hss*sig^2
    %     y_(t+1) = ybar + gx*(xf(t)+xs(t)) + gu*u(t+1)
    % Set ABCD representation:
    % z(t) = c + A*z(t-1) + B*xi(t+1)
    % y(t) = ybar + d + C*z(t-1) + D*xi(t+1)
    % z(t)=[xf(t);xs(t);kron(xf(t),xf(t))] 
    % xi(t+1)=[u(t+1);kron(u(t+1),u(t+1))-vec(Sigma);kron(xf(t),u(t+1));kron(u(t+1),xf(t))]
    nz = 2*nx + nx^2;
    nxi = nu+nu^2+2*nx*nu;
    nximin=nu+nu*(nu+1)/2+nu*nx;
    hx_hx = kron(hx,hx); hx_hu = kron(hx,hu); hu_hx = kron(hu,hx); hu_hu = kron(hu,hu);
    
    A = zeros(nz,nz);
    B = zeros(nz,nxi);
    C = zeros(ny,nz);
    D = zeros(ny,nxi);
    c = zeros(nz,1);
    d = zeros(ny,1);
    
    A(id1_xf,id1_xf) = hx;
    A(id2_xs,id2_xs) = hx;
    A(id2_xs,id3_xf_xf) = 0.5*Hxx;
    A(id3_xf_xf,id3_xf_xf) = hx_hx;
        
    B(id1_xf,id1_u) = hu;
    B(id2_xs,id2_u_u) = 1/2*Huu;
    B(id2_xs,id3_xf_u) = Hxu;
    B(id3_xf_xf,id2_u_u) = hu_hu;
    B(id3_xf_xf,id3_xf_u) = hx_hu;
    B(id3_xf_xf,id4_u_xf) = hu_hx;
    
    C(1:ny,id1_xf) = gx;
    C(1:ny,id2_xs) = gx;
    C(1:ny,id3_xf_xf) = 1/2*Gxx;
    
    D(1:ny,id1_u) = gu;
    D(1:ny,id2_u_u) = 1/2*Guu;
    D(1:ny,id3_xf_u) = Gxu;

    c(id2_xs,1) = 1/2*hss*sigm^2 + 1/2*Huu*M2u; 
    c(id3_xf_xf,1) =hu_hu*M2u;

    d(1:ny,1) = 1/2*gss*sigm^2 + 1/2*Guu*M2u;
    
    % Reduce dimensions
    if setup.dimreduce        
        nx2 = nx*(nx+1)/2;
        nzmin = 2*nx + nx2;    

        [DPx,DPxinv] = duplication(nx);
        Inx = eye(nx);
        Fz = zeros(nz,nzmin);    
        Fz(id1_xf,id1_xf) = Inx;
        Fz(id2_xs,id2_xs) = Inx;
        Fz(id3_xf_xf,id2_xs(end)+(1:nx2)) = DPx;

        Fzinv = zeros(nzmin,nz);
        Fzinv(id1_xf,id1_xf) = Inx;
        Fzinv(id2_xs,id2_xs) = Inx;
        Fzinv(id2_xs(end)+(1:nx2),id3_xf_xf) = DPxinv;

        A = Fzinv*A*Fz;
        B = Fzinv*B;
        C = C*Fz;
        c = Fzinv*c;
        disp('DIMENSION REDUCTION 2nd order')
        disp(nz);
        disp(nzmin);
    else
        Fz = [];
    end
    
elseif orderApp == 3
    % Third-order state-space representation for all endogenous variables y given first-order states xf, second-order states xs, third-order states xrd and shocks u with E(u_t u_t') = sigma^2*eta*eta' = Sigma
    % xf(t+1) = hx*xf(t) + hu*u(t+1)
    % xs(t+1) = hx*xs(t) + 1/2*Hxx*kron(xf(t),xf(t)) 
    %           + 1/2*Huu*kron(u(t+1),u(t+1)) + Hxu*kron(xf(t),u(t+1)) + 1/2*hss*sig^2
    % xrd(t+1} = hx*xrd(t) + Hxx*kron(xf(t),xs(t)) + Hxu*kron(xs(t),u_(t+1))
    %             +3/6*Hxss*sigm^2*xf(t) + 3/6*Huss*u(t+1)
    %             +1/6*Hxxx*kron(xf(t),xf(t),xf(t)) + 1/6*Huuu*(kron(u(t+1),u(t+1),u(t+1))-M3u+M3u)
    %             +3/6*Hxxu*kron(xf(t),xf(t),u(t+1)) + 3/6*Hxuu*kron(xf(t),u(t+1),u(t+1))
    %             +1/6*Hsss*sigm^3
    % y(t+1) = ybar + gx*(xf(t)+xs(t)+xrd(t)) + gu*u(t+1)
    %        + 1/2*Gxx*kron(xf(t),xf(t)) + 1/2*Guu*kron(u(t+1),u(t+1)) + Gxu*kron(xf(t),u(t+1)) + 1/2*gss*sig^2
    %        + Gxx*kron(xf(t),xs(t)) + Gxu*kron(xs(t),u_(t+1))
    %        +3/6*Gxss*sigm^2*xf(t) + 3/6*Guss*u(t+1)
    %        +1/6*Gxxx*kron(xf(t),xf(t),xf(t)) + 1/6*Guuu*(kron(u(t+1),u(t+1),u(t+1))-M3u+M3u)
    %        +3/6*Gxxu*kron(xf(t),xf(t),u(t+1)) + 3/6*Gxuu*kron(xf(t),u(t+1),u(t+1))
    %        +1/6*Gsss*sigm^3
    %
    % Set ABCD representation:
    % z(t) = c + A*z(t-1) + B*xi(t+1)
    % y(t) = ybar + d + C*z(t-1) + D*xi(t+1)
    % z(t)=[xf(t);xs(t);kron(xf(t),xf(t));xrd(t);kron(xf(t),xs(t));kron(xf(t),xf(t),xf(t))] 
    %  xi(t+1) = [u(t+1); kron(u(t+1),u(t+1))-M2u; kron(xf(t),u(t+1)); kron(u(t+1),xf(t)); kron(xs(t),u(t+1));kron(u(t+1),xs(t));...
    %       kron(xf(t),xf(t),u(t+1));kron(xf(t),u(t+1),xf(t));kron(u(t+1),xf(t),xf(t));...
    %       kron(xf(t),u(t+1),u(t+1));kron(u(t+1),xf(t),u(t+1)); kron(u(t+1),u(t+1),xf(t));...
    %       kron(u(t+1),u(t+1),u(t+1))-M3u];
    nz = 3*nx + 2*nx^2 +nx^3;
    nxi = nu+nu^2+2*nx*nu+2*nx*nu+3*nx^2*nu+3*nu^2*nx+nu^3;    
    nu2 = nu*(nu+1)/2; nx2 = nx*(nx+1)/2; nu3 = nu2*(nu+2)/3;
    nximin = nu + nu2 + 2*nu*nx + nu*nx2 + nu2*nx + nu3;
    
    hx_hx = kron(hx,hx); hu_hu = kron(hu,hu); hx_hu=kron(hx,hu); hu_hx = kron(hu,hx);
    A = zeros(nz,nz);
    B = zeros(nz,nxi);
    C = zeros(ny,nz);
    D = zeros(ny,nxi);
    c = zeros(nz,1);
    d = zeros(ny,1);
    
    A(id1_xf,id1_xf) = hx;
    A(id2_xs,id2_xs) = hx;
    A(id2_xs,id3_xf_xf) = 1/2*Hxx;
    A(id3_xf_xf,id3_xf_xf) = hx_hx;
    A(id4_xrd,id1_xf) = 3/6*Hxss*sigm^2;
    A(id4_xrd,id4_xrd) = hx;
    A(id4_xrd,id5_xf_xs) = Hxx;
    A(id4_xrd,id6_xf_xf_xf) = 1/6*Hxxx;
    A(id5_xf_xs,id1_xf) = kron(hx,1/2*hss*sigm^2);
    A(id5_xf_xs,id5_xf_xs) = hx_hx;
    A(id5_xf_xs,id6_xf_xf_xf) = kron(hx,1/2*Hxx);
    A(id6_xf_xf_xf,id6_xf_xf_xf) = kron(hx,hx_hx);
    
    B(id1_xf,id1_u) = hu;    
    B(id2_xs,id2_u_u) = 1/2*Huu;
    B(id2_xs,id3_xf_u) = Hxu;    
    B(id3_xf_xf,id2_u_u) = hu_hu;
    B(id3_xf_xf,id3_xf_u) = hx_hu;
    B(id3_xf_xf,id4_u_xf) = hu_hx;    
    B(id4_xrd,id1_u) = 3/6*Huss*sigm^2;
    B(id4_xrd,id5_xs_u) = Hxu;
    B(id4_xrd,id7_xf_xf_u) = 3/6*Hxxu;
    B(id4_xrd,id10_xf_u_u) = 3/6*Hxuu;
    B(id4_xrd,id13_u_u_u) =  1/6*Huuu;    
    B(id5_xf_xs,id1_u) = kron(hu,1/2*hss*sigm^2);
    B(id5_xf_xs,id6_u_xs) =  hu_hx;
    B(id5_xf_xs,id7_xf_xf_u) = kron(hx,Hxu);
    B(id5_xf_xs,id9_u_xf_xf) = kron(hu,1/2*Hxx);
    B(id5_xf_xs,id10_xf_u_u) = kron(hx,1/2*Huu);
    B(id5_xf_xs,id11_u_xf_u) = kron(hu,Hxu);
    B(id5_xf_xs,id13_u_u_u) = kron(hu,1/2*Huu);    
    B(id6_xf_xf_xf,id7_xf_xf_u) =  kron(hx_hx,hu);
    B(id6_xf_xf_xf,id8_xf_u_xf) = kron(hx,hu_hx);
    B(id6_xf_xf_xf,id9_u_xf_xf) = kron(hu,hx_hx);
    B(id6_xf_xf_xf,id10_xf_u_u) = kron(hx_hu,hu);
    B(id6_xf_xf_xf,id11_u_xf_u) = kron(hu,hx_hu);
    B(id6_xf_xf_xf,id12_u_u_xf) = kron(hu_hu,hx);
    B(id6_xf_xf_xf,id13_u_u_u)  = kron(hu,hu_hu);         
    
    C(1:ny,id1_xf) = gx+.5*Gxss*sigm^2;
    C(1:ny,id2_xs) = gx;
    C(1:ny,id3_xf_xf) = 0.5*Gxx;
    C(1:ny,id4_xrd) = gx;
    C(1:ny,id5_xf_xs) = Gxx;
    C(1:ny,id6_xf_xf_xf) = 1/6*Gxxx;
    
    D(1:ny,id1_u) = gu+.5*Guss*sigm^2;
    D(1:ny,id2_u_u) = 0.5*Guu;
    D(1:ny,id3_xf_u) = Gxu;
    D(1:ny,id5_xs_u) = Gxu;    
    D(1:ny,id7_xf_xf_u) = 1/2*Gxxu;
    D(1:ny,id10_xf_u_u) = 1/2*Gxuu;
    D(1:ny,id13_u_u_u) = 1/6*Guuu;
    
    c(id2_xs,1) = 1/2*hss*sigm^2 + 1/2*Huu*M2u;
    c(id3_xf_xf,1) = hu_hu*M2u; 
    c(id4_xrd,1) = 1/6*Huuu*M3u + 1/6*Hsss*sigm^3; 
    c(id5_xf_xs,1) =  kron(hu,1/2*Huu)*M3u; 
    c(id6_xf_xf_xf,1) = kron(hu_hu,hu)*M3u;
    
    d(1:ny,1) = 0.5*gss*sigm^2 + 0.5*Guu*M2u + 1/6*Guuu*M3u + 1/6*Gsss*sigm^3;
    
    % Reduce dimensions
    if setup.dimreduce
        nx2 = nx*(nx+1)/2; nx3 = nx*(nx+1)*(nx+2)/(2*3);
        nzmin = 3*nx + nx2 + nx^2 + nx3;

        [DPx,DPxinv] = duplication(nx);
        [TPx,TPxinv] = triplication(nx);
        Inx = eye(nx); Inx2 = eye(nx^2);
        Fz = zeros(nz,nzmin);    
        Fz(id1_xf,id1_xf) = Inx;
        Fz(id2_xs,id2_xs) = Inx;
        Fz(id3_xf_xf,id2_xs(end)+(1:nx2)) = DPx;
        Fz(id4_xrd,id2_xs(end)+nx2+(1:nx)) = Inx;
        Fz(id5_xf_xs,id2_xs(end)+nx2+nx+(1:nx^2)) = Inx2;
        Fz(id6_xf_xf_xf,id2_xs(end)+nx2+nx+nx^2+(1:nx3)) = TPx;

        Fzinv = zeros(nzmin,nz);
        Fzinv(id1_xf,id1_xf) = Inx;
        Fzinv(id2_xs,id2_xs) = Inx;
        Fzinv(id2_xs(end)+(1:nx2),id3_xf_xf) = DPxinv;
        Fzinv(id2_xs(end)+nx2+(1:nx),id4_xrd) = Inx;
        Fzinv(id2_xs(end)+nx2+nx+(1:nx^2),id5_xf_xs) = Inx2;
        Fzinv(id2_xs(end)+nx2+nx+nx^2+(1:nx3),id6_xf_xf_xf) = TPxinv;

        A = Fzinv*A*Fz;
        B = Fzinv*B;
        C = C*Fz;
        c = Fzinv*c;
        disp('DIMENSION REDUCTION 3rd order')
        disp(nz);
        disp(nzmin);
    else
        Fz =[];
    end
end

if sparseflag
    A = sparse(A);
    B = sparse(B);
    C = sparse(C);
    D = sparse(D);
    c = sparse(c);
    d = sparse(d);
else
    A = full(A);
    B = full(B);
    C = full(C);
    D = full(D);
    c = full(c);
    d = full(d);
end