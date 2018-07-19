% By Willi Mutschler, September 26, 2016. Email: willi@mutschler.eu
function createHOSScriptFiles(nu,nx,orderApp,distrib,orderPM)
toolboxpath = fileparts(mfilename('fullpath'));

%% SYMBOLIC DERIVATIONS product moments for extended innovations
% These are model independent, they only dependent on the number of states
% nx and number of shocks nu
[DPu,DPuinv] = duplication(nu);
% Set distributional assumptions
switch distrib
    case 'Gaussian'
        %if u = epsi with epsi ~ N(0,SIGe), then u ~N(0,SIGe)
        SIGe=sym('SIGe_',[nu nu]); 
        Mom.SIGe = tril(SIGe,0) + tril(SIGe,-1).';
        Mom.E_uu = DPuinv*Mom.SIGe(:);
        Mom.E_uuu = sym(zeros(nu*(nu+1)*(nu+2)/6,1));
        auxparu = DPuinv*Mom.SIGe(:);
    case 'Student_t'
        % if u = v^(-1/2)*epsi with epsi~N(0,SIGe) and v ~ Gamma(df/2,2/df), epsi and v independent of each other,
        % then u is multivariate student t with df degrees of freedom and covariance df/(df-2)*SIGe
        SIGe= sym('SIGe_',[nu nu]);        
        Mom.SIGe = tril(SIGe,0) + tril(SIGe,-1).';
        Mom.df = sym('df','positive');
        Mom.E_uu = DPuinv*vec(Mom.df/(Mom.df-2)*Mom.SIGe);
        Mom.E_uuu = sym(zeros(nu*(nu+1)*(nu+2)/6,1));        
        auxparu = [Mom.df;DPuinv*Mom.SIGe(:)];
end


if orderApp == 1        
    nximin = nu;
    % Symbolic variables for shocks
    u = sym('u',[nu 1]); % Create symbolic variables for epsilon
    % Create minimal xi_t vector    
    ximin = u;    
%   Test whether this is correct:    
%    xi = [u];
%    [xi - GetFxi(nu,nx,orderApp,0)*ximin]
    auxpar =  auxparu; % These are the product moments we need to compute GAMMA2XI, GAMMA3XI or GAMMA4XI
elseif orderApp == 2
    nx2= nx*(nx+1)/2; nx3=nx2*(nx+2)/3; nx4=nx3*(nx+3)/4;     
    nu2 = nu*(nu+1)/2;        
    nximin = nu + nu2 + nu*nx;    
    % Symbolic variables for shocks
    u = sym('u',[nu 1]); % Create symbolic variables for epsilon
    xf = sym('xf',[nx 1]); % Symbolic variables for first-order terms xf    
    E_uu   = sym('E_uu_',[nu2 1]); % unique elements for second-order product moments of epsi
    % Create minimal xi_t vector    
    u_u = kron(u,u); 
    xf_u = kron(xf,u);     
    ximin = [u; DPuinv*(u_u)-E_uu; xf_u];
%   Test whether this is correct:
%     xi = [u; 
%          u_u-DPu*E_uu;          
%          kron(xf,u);
%          kron(u,xf)];
%     [xi - GetFxi(nu,nx,orderApp,0)*ximin]
    
    % Create symbolic variables for unique terms of product moments of xf and xs    
    Mom.E_XF{1} = sym('E_XF1_',[nx 1]); %E(xf)
    Mom.combos.E_XF{1} = flipud(allVL1(nx, 1)); % matrix of powers for unique elelments
    Mom.E_XF{2} = sym('E_XF2_',[nx2 1]); % unique values of E(kron(xf,xf))
    Mom.combos.E_XF{2} = flipud(allVL1(nx, 2)); % matrix of powers for unique elelments
    if orderPM > 2
        Mom.E_XF{3} = sym('E_XF3_',[nx3 1]); % unique values of E(kron(xf,xf,xf))
        Mom.combos.E_XF{3} = flipud(allVL1(nx, 3)); % matrix of powers for unique elelments
    end
    if orderPM > 3
        Mom.E_XF{4} = sym('E_XF4_',[nx4 1]); % unique values of E(kron(xf,xf,xf,xf))
        Mom.combos.E_XF{4} = flipud(allVL1(nx, 4)); % matrix of powers for unique elelments
    end
    
    if orderPM == 2
        auxpar =  [auxparu;Mom.E_XF{1}; Mom.E_XF{2}]; % These are the product moments we need to compute GAMMA2XI
    elseif orderPM == 3
        auxpar  = [auxparu;Mom.E_XF{1}; Mom.E_XF{2}; Mom.E_XF{3}]; % These are the product moments we need to compute GAMMA3XI
    elseif orderPM == 4
        auxpar =  [auxparu;Mom.E_XF{1}; Mom.E_XF{2}; Mom.E_XF{3}; Mom.E_XF{4}]; % These are the product moments we need to compute GAMMA4XI
    end
elseif orderApp == 3
    nx2= nx*(nx+1)/2; nx3=nx2*(nx+2)/3; nx4=nx3*(nx+3)/4; 
    nx5=nx4*(nx+4)/5; nx6=nx5*(nx+5)/6; nx7=nx6*(nx+6)/7; nx8=nx7*(nx+7)/8; 
    nu2 = nu*(nu+1)/2; nu3 = nu2*(nu+2)/3;
    nximin = nu + nu2 + 2*nu*nx + nu*nx2 + nu2*nx + nu3;
    [DPx,DPxinv] = duplication(nx);
    [TPu,TPuinv] = triplication(nu);
    I_nu = eye(nu);
    I_nx = eye(nx);
    % Symbolic variables for shocks
    u = sym('u',[nu 1]); % Create symbolic variables for epsilon
    xf = sym('xf',[nx 1]); % Symbolic variables for first-order terms xf
    xs = sym('xs',[nx 1]); % Symbolic variables for second-order terms xf
    E_uu   = sym('E_uu_',[nu2 1]); % unique elements for second-order product moments of epsi
    E_uuu = sym('E_uuu_',[nu3 1]); % unique elements for third-order product moments of epsi    
    % Create minimal xi_t vector    
    u_u = kron(u,u); 
    xf_u = kron(xf,u); 
    xs_u = kron(xs,u);
    xf_xf_u = kron(xf,xf_u);
    xf_u_u = kron(xf,u_u);
    u_u_u = kron(u,u_u);
    ximin = [u; DPuinv*(u_u)-E_uu;   xf_u;   xs_u;   kron(DPxinv,I_nu)*xf_xf_u;   kron(I_nx,DPuinv)*xf_u_u;   TPuinv*u_u_u-E_uuu];
%   Check whether this is corret    
%     xi = [u; 
%         kron(u,u)-DPu*E_uu; 
%         kron(xf,u); 
%         kron(u,xf); 
%         kron(xs,u);
%         kron(u,xs);
%         kron(kron(xf,xf),u); 
%         kron(xf,kron(u,xf));
%         kron(u,kron(xf,xf));         
%         kron(xf,kron(u,u)); 
%         kron(u,kron(xf,u));
%         kron(u,kron(u,xf));
%         kron(u,kron(u,u))-TPu*E_uuu];    
%    [xi - GetFxi(nu,nx,orderApp,0)*ximin]

    % Create symbolic variables for unique terms of product moments of xf and xs
    Mom.E_XF{1} = sym('E_XF1_',[nx 1]); %E(xf)
    Mom.combos.E_XF{1} = flipud(allVL1(nx, 1)); % matrix of powers for unique elelments
    Mom.E_XF{2} = sym('E_XF2_',[nx2 1]); % unique values of E(kron(xf,xf))
    Mom.combos.E_XF{2} = flipud(allVL1(nx, 2)); % matrix of powers for unique elelments
    Mom.E_XF{3} = sym('E_XF3_',[nx3 1]); % unique values of E(kron(xf,xf,xf))
    Mom.combos.E_XF{3} = flipud(allVL1(nx, 3)); % matrix of powers for unique elelments
    Mom.E_XF{4} = sym('E_XF4_',[nx4 1]); % unique values of E(kron(xf,xf,xf,xf))
    Mom.combos.E_XF{4} = flipud(allVL1(nx, 4)); % matrix of powers for unique elelments
    Mom.E_XS{1} = sym('E_XS1_',[nx 1]); % E(xs)
    Mom.E_XS{2} = sym('E_XS2_',[nx2 1]); % unqiue values of E(kron(xs,xs))
    Mom.combos.E_XS{2} = flipud(allVL1(nx, 2)); % matrix of powers for unique elelments
    Mom.E_XF_XS{1,1} = sym('E_XF1_XS1_',[nx*nx 1]);  % unqiue values of E(kron(xf,xs))
    Mom.combos.E_XF_XS{1,1} = GetPowers(eye(nx),eye(nx));% matrix of powers for unique elelments        
    Mom.E_XF_XS{2,1} = sym('E_XF2_XS1_',[nx2*nx 1]); % unqiue values of E(kron(xf,xf,xs))
    Mom.combos.E_XF_XS{2,1} = GetPowers(flipud(allVL1(nx, 2)),eye(nx));% matrix of powers for unique elelments 
    
    if orderPM > 2    
        Mom.E_XF{5} = sym('E_XF5_',[nx5 1]); % unqiue values of E(kron(xf,xf,xf,xf,xf))
        Mom.combos.E_XF{5} = flipud(allVL1(nx, 5)); % matrix of powers for unique elelments   
        Mom.E_XF{6} = sym('E_XF6_',[nx6 1]); % unqiue values of E(kron(xf,xf,xf,xf,xf,xf))
        Mom.combos.E_XF{6} = flipud(allVL1(nx, 6)); % matrix of powers for unique elelments
        Mom.E_XS{3} = sym('E_XS3_',[nx3 1]); % unqiue values of E(kron(xs,xs,xs))
        Mom.combos.E_XS{3} = flipud(allVL1(nx, 3)); % matrix of powers for unique elelments
        Mom.E_XF_XS{1,2} = sym('E_XF1_XS2_',[nx*nx2 1]); % unqiue values of E(kron(xf,xs,xs))
        Mom.combos.E_XF_XS{1,2} = GetPowers(eye(nx),flipud(allVL1(nx, 2)));% matrix of powers for unique elelments         
        Mom.E_XF_XS{3,1} = sym('E_XF3_XS1_',[nx3*nx 1]);  % unqiue values of E(kron(xf,xf,xf,xs))
        Mom.combos.E_XF_XS{3,1} = GetPowers(flipud(allVL1(nx, 3)),eye(nx));% matrix of powers for unique elelments
        Mom.E_XF_XS{2,2} = sym('E_XF2_XS2_',[nx2*nx2 1]);  % unqiue values of E(kron(xf,xf,xs,xs))
        Mom.combos.E_XF_XS{2,2} = GetPowers(flipud(allVL1(nx, 2)),flipud(allVL1(nx, 2)));% matrix of powers for unique elelments 
        Mom.E_XF_XS{4,1} = sym('E_XF4_XS1_',[nx4*nx 1]);    % unqiue values of E(kron(xf,xf,xf,xf,xs))
        Mom.combos.E_XF_XS{4,1} = GetPowers(flipud(allVL1(nx, 4)),eye(nx)); % matrix of powers for unique elelments        
    end
    if orderPM > 3
        Mom.E_XF{7} = sym('E_XF7_',[nx7 1]); % unqiue values of E(kron(xf,xf,xf,xf,xf,xf,xf))
        Mom.combos.E_XF{7} = flipud(allVL1(nx, 7)); % matrix of powers for unique elelments 
        Mom.E_XF{8} = sym('E_XF8_',[nx8 1]); % unqiue values of E(kron(xf,xf,xf,xf,xf,xf,xf,xf))
        Mom.combos.E_XF{8} = flipud(allVL1(nx, 8)); % matrix of powers for unique elelments         
        Mom.E_XS{4} = sym('E_XS4_',[nx4 1]); % unqiue values of E(kron(xs,xs,xs,xs))
        Mom.combos.E_XS{4} = flipud(allVL1(nx, 4)); % matrix of powers for unique elelments
        Mom.E_XF_XS{1,3} = sym('E_XF1_XS3_',[nx*nx3 1]);    % unqiue values of E(kron(xf,xs,xs,xs))
        Mom.combos.E_XF_XS{1,3} = GetPowers(eye(nx),flipud(allVL1(nx, 3))); % matrix of powers for unique elelments       
        Mom.E_XF_XS{3,2} = sym('E_XF3_XS2_',[nx3*nx2 1]);    % unqiue values of E(kron(xf,xf,xf,xs,xs))  
        Mom.combos.E_XF_XS{3,2} = GetPowers(flipud(allVL1(nx, 3)),flipud(allVL1(nx, 2)));% matrix of powers for unique elelments 
        Mom.E_XF_XS{2,3} = sym('E_XF2_XS3_',[nx2*nx3 1]);    % unqiue values of E(kron(xf,xf,xs,xs,xs))
        Mom.combos.E_XF_XS{2,3} = GetPowers(flipud(allVL1(nx, 2)),flipud(allVL1(nx, 3)));% matrix of powers for unique elelments 
        Mom.E_XF_XS{5,1} = sym('E_XF5_XS1_',[nx5*nx 1]);    % unqiue values of E(kron(xf,xf,xf,xf,xf,xs))
        Mom.combos.E_XF_XS{5,1} = GetPowers(flipud(allVL1(nx, 5)),eye(nx));% matrix of powers for unique elelments 
        Mom.E_XF_XS{4,2} = sym('E_XF4_XS2_',[nx4*nx2 1]);    % unqiue values of E(kron(xf,xf,xf,xf,xs,xs))  
        Mom.combos.E_XF_XS{4,2} = GetPowers(flipud(allVL1(nx, 4)),flipud(allVL1(nx, 2)));% matrix of powers for unique elelments
        Mom.E_XF_XS{6,1} = sym('E_XF6_XS1_',[nx6*nx 1]);    % unqiue values of E(kron(xf,xf,xf,xf,xf,xf,xs))   
        Mom.combos.E_XF_XS{6,1} = GetPowers(flipud(allVL1(nx, 6)),eye(nx));% matrix of powers for unique elelments     
    end
    
    if orderPM == 2
        auxpar =  [auxparu;
                    Mom.E_XF{1};Mom.E_XF{2}; Mom.E_XF{3}; Mom.E_XF{4}; 
                    Mom.E_XS{1};Mom.E_XS{2}; 
                    Mom.E_XF_XS{1,1}; Mom.E_XF_XS{2,1}]; % These are the product moments we need to compute GAMMA2XI
    elseif orderPM == 3
        auxpar  = [auxparu;
                Mom.E_XF{1};Mom.E_XF{2}; Mom.E_XF{3}; Mom.E_XF{4}; Mom.E_XF{5}; Mom.E_XF{6}; 
                Mom.E_XS{1};Mom.E_XS{2}; Mom.E_XS{3}; 
                Mom.E_XF_XS{1,1}; Mom.E_XF_XS{2,1}; Mom.E_XF_XS{1,2};
                Mom.E_XF_XS{3,1}; Mom.E_XF_XS{2,2}; 
                Mom.E_XF_XS{4,1}]; % These are the product moments we need to compute GAMMA3XI
    elseif orderPM == 4
        auxpar =  [auxparu;
                Mom.E_XF{1};Mom.E_XF{2}; Mom.E_XF{3}; Mom.E_XF{4}; Mom.E_XF{5}; Mom.E_XF{6}; Mom.E_XF{7}; Mom.E_XF{8};
                Mom.E_XS{1};Mom.E_XS{2}; Mom.E_XS{3}; Mom.E_XS{4};
                Mom.E_XF_XS{1,1}; Mom.E_XF_XS{2,1}; Mom.E_XF_XS{1,2};
                Mom.E_XF_XS{3,1}; Mom.E_XF_XS{2,2}; Mom.E_XF_XS{1,3};
                Mom.E_XF_XS{4,1}; Mom.E_XF_XS{3,2}; Mom.E_XF_XS{2,3};
                Mom.E_XF_XS{5,1}; Mom.E_XF_XS{4,2};
                Mom.E_XF_XS{6,1}]; % These are the product moments we need to compute GAMMA4XI
    end
end

if orderPM == 2; 
    nunique  = nximin*(nximin+1)/2;
elseif orderPM == 3; 
    nunique  = nximin*(nximin+1)*(nximin+2)/6;    
elseif orderPM == 4; 
    nunique  = nximin*(nximin+1)*(nximin+2)*(nximin+3)/24;
end

%% Create and initialize files
positionFiles     = [toolboxpath,'\OrderApp_',num2str(orderApp)];
if orderApp == 1
    nameOfFunction = sprintf('Xi_%s_prodmom%i_orderApp%i_nu%i',distrib,orderPM,orderApp,nu);
else
    nameOfFunction = sprintf('Xi_%s_prodmom%i_orderApp%i_nx%i_nu%i',distrib,orderPM,orderApp,nx,nu);
end

saveFile = [positionFiles,'\',nameOfFunction,'.m'];

if exist(saveFile,'file') > 0; delete(saveFile); end

text = sprintf('function nXI%imin = %s(arg)',orderPM,nameOfFunction);
dlmwrite(saveFile,text,'-append','delimiter','');
for j=1:length(auxpar)
    text = sprintf('%s = arg(%d);',char(auxpar(j)),j);
    dlmwrite(saveFile,text,'-append','delimiter','');
end

text=sprintf('nXI%imin=zeros(%d,1);',orderPM,nunique);
dlmwrite(saveFile,text,'-append','delimiter','');

%% Computations
%Get Combos
nameOfFunctionCombos = sprintf('combos%i_%i.mat',orderPM,nximin);
saveFileCombos = [toolboxpath,'\Combos\',nameOfFunctionCombos];
if exist(saveFileCombos,'file') == 0
    combos   = sparse(flipud(allVL1(nximin, orderPM)));% All integer permutations with sum criteria == orderPM
    save(saveFileCombos, 'combos');
else
    load(saveFileCombos, 'combos'); 
end

reverseStr = '';
for j = 1:nunique
    if (rem(j,25)== 0)
        msg = sprintf('    %i-order product moments processed %d/%d', orderPM, j, nunique); fprintf([reverseStr, msg]); reverseStr = repmat(sprintf('\b'), 1, length(msg));
    elseif (j==nunique)
        msg = sprintf('    %i-order product moments processed %d/%d\n', orderPM, j, nunique); fprintf([reverseStr, msg]); reverseStr = repmat(sprintf('\b'), 1, length(msg));
    end
    xikmin = expand(prod(power(ximin,combos(j,:)')));% evaluate string of product moments of each element in xi_t
    XIkmin = ProdMom_inov(xikmin,nu,nx,Mom,distrib);
    if (XIkmin ~= 0)
        text = sprintf('nXI%imin(%d,1) = %s;',orderPM,j,char(XIkmin));
        dlmwrite(saveFile,text,'-append','delimiter','');
    end        
end