% By Willi Mutschler, September 26, 2016. Email: willi@mutschler.eu
function Fxi = GetFxi(nu,nx,orderApp,sparseflag)
if nargin < 4
    sparseflag = 1;
end

if orderApp == 1
    nxi = nu;
    nximin = nu;    
    col1_u       = 1:nu;
    row1_u       = 1:nu;
    
    if sparseflag
        Iu = speye(nu);
        Fxi = spalloc(nxi,nximin,nu);
    else
        Iu = eye(nu);
        Fxi = zeros(nxi,nximin);
    end    
    Fxi(row1_u,col1_u) = Iu;

elseif orderApp == 2    
    nu2 = nu*(nu+1)/2;
    nxi = nu + nu^2 + 2*nx*nu;
    nximin = nu + nu2 + nu*nx;
    
    col1_u       = 1:nu;
    col2_u_u     = col1_u(end)   + (1:nu2);
    col3_xf_u    = col2_u_u(end) + (1:nu*nx);
    
    row1_u       = 1:nu;
    row2_u_u     = row1_u(end)    + (1:nu^2);    
    row3_xf_u    = row2_u_u(end)  + (1:nu*nx);
    row4_u_xf    = row3_xf_u(end) + (1:nx*nu);
    
    [DPu,~] = duplication(nu,sparseflag);
    K_u_x = commutation(nu,nx,sparseflag);
    if sparseflag
        Iu = speye(nu); 
        Iux = speye(nu*nx);
        Fxi = spalloc(nxi,nximin,nu+nu^2+2*nu*nx);
    else
        Iu = eye(nu); 
        Iux = eye(nu*nx);
        Fxi = zeros(nxi,nximin);
    end   
    
    Fxi(row1_u,col1_u) = Iu; 
    Fxi(row2_u_u,col2_u_u) = DPu; 
    Fxi(row3_xf_u,col3_xf_u) = Iux; 
    Fxi(row4_u_xf,col3_xf_u) = K_u_x; 
    
elseif orderApp == 3
    nu2 = nu*(nu+1)/2;  nu3=nu2*(nu+2)/3;
    nx2 = nx*(nx+1)/2;  
    nxi = nu + nu^2 + 3*nx*nu + 3*nu*nx^2 + 3*nu^2*nx + nu^3;
    nximin = nu + nu2 + 2*nu*nx + nu*nx2 + nu2*nx + nu3;
    
    col1_u       = 1:nu;
    col2_u_u     = col1_u(end)       + (1:nu2);
    col3_xf_u    = col2_u_u(end)     + (1:nu*nx);
    col4_xs_u    = col3_xf_u(end)    + (1:nu*nx);
    col5_xf_xf_u = col4_xs_u(end)    + (1:nu*nx2);
    col6_xf_u_u  = col5_xf_xf_u(end) + (1:nu2*nx);
    col7_u_u_u   = col6_xf_u_u(end)  + (1:nu3);
    
    row1_u       = 1:nu;
    row2_u_u     = row1_u(end)       + (1:nu^2);    
    row3_xf_u    = row2_u_u(end)     + (1:nx*nu);
    row4_u_xf    = row3_xf_u(end)    + (1:nx*nu);
    row5_xs_u    = row4_u_xf(end)    + (1:nx*nu);
    row6_u_xs    = row5_xs_u(end)    + (1:nx*nu);    
    row7_xf_xf_u = row6_u_xs(end)    + (1:nu*nx^2);
    row8_xf_u_xf = row7_xf_xf_u(end) + (1:nu*nx^2);
    row9_u_xf_xf = row8_xf_u_xf(end) + (1:nu*nx^2);    
    row10_xf_u_u = row9_u_xf_xf(end) + (1:nx*nu^2);
    row11_u_xf_u = row10_xf_u_u(end) + (1:nx*nu^2);
    row12_u_u_xf = row11_u_xf_u(end) + (1:nx*nu^2);
    row13_u_u_u  = row12_u_u_xf(end) + (1:nu^3);    
    
    [DPx,] = duplication(nx,sparseflag);
    [DPu,] = duplication(nu,sparseflag);
    [TPu,~] = triplication(nu,0,sparseflag);    
    K_u_x = commutation(nu,nx,sparseflag);    
    K_u_xx = commutation(nu,nx^2,sparseflag);
    K_u_xu = commutation(nu,nu*nx,sparseflag);    
    K_ux_x = commutation(nu*nx,nx,sparseflag);    
    K_uu_x = commutation(nu^2,nx,sparseflag);
    if sparseflag
        Ix = speye(nx);
        Iu = speye(nu); 
        Iux = speye(nu*nx);
        Fxi = spalloc(nxi,nximin,nu+nu^2+4*nu*nx+3*nx^2*nu+3*nx*nu^2+nu^3);
    else
        Ix = eye(nx);
        Iu = eye(nu);
        Iux = eye(nu*nx);
        Fxi = zeros(nxi,nximin);
    end
    DPx_Iu = kron(DPx,Iu);
    Ix_DPu = kron(Ix,DPu);

    Fxi(row1_u,col1_u) = Iu;
    Fxi(row2_u_u,col2_u_u) = DPu;
    Fxi(row3_xf_u,col3_xf_u) = Iux;
    Fxi(row4_u_xf,col3_xf_u) = K_u_x; 
    Fxi(row5_xs_u,col4_xs_u) = Iux;
    Fxi(row6_u_xs,col4_xs_u) = K_u_x;
    Fxi(row7_xf_xf_u,col5_xf_xf_u) = DPx_Iu;
    Fxi(row8_xf_u_xf,col5_xf_xf_u) = K_ux_x*DPx_Iu;
    Fxi(row9_u_xf_xf,col5_xf_xf_u) = K_u_xx*DPx_Iu;
    Fxi(row10_xf_u_u,col6_xf_u_u) = Ix_DPu;
    Fxi(row11_u_xf_u,col6_xf_u_u) = K_u_xu*Ix_DPu;
    Fxi(row12_u_u_xf,col6_xf_u_u) = K_uu_x*Ix_DPu;
    Fxi(row13_u_u_u,col7_u_u_u) = TPu;
    
    
end