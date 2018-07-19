% By Willi Mutschler, September 26, 2016. Email: willi@mutschler.eu
function xinames = GetNamesXi(orderApprox,indx_states,nx,nu)
global M_ oo_
namesu = deblank(M_.exo_names); 
if orderApprox > 1
    namesxf = deblank(strcat(M_.endo_names(oo_.dr.order_var(indx_states),:),repmat('^f',nx,1)));
    if orderApprox > 2
        namesxs = deblank(strcat(M_.endo_names(oo_.dr.order_var(indx_states),:),repmat('^s',nx,1)));
    end
end

if orderApprox == 1
    xinames = namesu;
elseif orderApprox == 2    
    xinames =[Get_u(namesu);Get_u_u_min(namesu,nu);Get_xf_u(namesxf,nx,namesu,nu)];
elseif orderApprox == 3    
    xinames =[Get_u(namesu);Get_u_u_min(namesu,nu);Get_xf_u(namesxf,nx,namesu,nu);Get_xs_u(namesxf,nx,namesu,nu);...
              Get_xf_xf_u(namesxf,nx,namesu,nu);Get_xf_u_u(namesxf,nx,namesu,nu);Get_u_u_u_min(namesu,nu)];
end


function xinames_u = Get_u(namesu)
    xinames_u = namesu;    
end

function xinames_u_u_min = Get_u_u_min(namesu,nu)
% u_u(idxvechu)-E_u_u
    countuu = 1; xinames_u_u_min = cell(nu*(nu+1)/2,1);
    for i=1:nu
        for j=i:nu
            xinames_u_u_min{countuu,1} = [namesu(i,:), '*', namesu(j,:), '-M2(',namesu(i,:), '*', namesu(j,:),')'];  
            countuu = countuu+1;
        end
    end
end

function xinames_xf_u = Get_xf_u(namesxf,nx,namesu,nu)
    % kron(xf,u)
    countxfu = 1; xinames_xf_u = cell(nx*nu,1);
    for i=1:nx;
        for j=1:nu
            xinames_xf_u{countxfu,1} = [namesxf(i,:), '*', namesu(j,:)];
            countxfu = countxfu+1;             
        end
    end
end

function xinames_xs_u = Get_xs_u(namesxs,nx,namesu,nu)
    % kron(xs,u)    
    countxsu = 1; xinames_xs_u = cell(nx*nu,1);
    for i=1:nx;
        for j=1:nu
            xinames_xs_u{countxsu,1} = [namesxs(i,:), '*', namesu(j,:)];
            countxsu = countxsu+1;             
        end
    end
end

function xinames_xf_xf_u = Get_xf_xf_u(namesxf,nx,namesu,nu)
    %kron(xf,xf,u)
    countxfxfu = 1; xinames_xf_xf_u = cell(nx^2*nu,1);
    for i=1:nx
       for j=1:nx
            for k=1:nu
                xinames_xf_xf_u{countxfxfu,1} = [namesxf(i,:), '*', namesxf(j,:), '*', namesu(k,:)];
                countxfxfu = countxfxfu+1; 
            end
       end
    end
end

function xinames_xf_u_u = Get_xf_u_u(namesxf,nx,namesu,nu)
    %kron(xf,u,u)
    countxfuu = 1; xinames_xf_u_u = cell(nx*nu^2,1);
    for i=1:nx
       for j=1:nu
            for k=1:nu
                xinames_xf_u_u{countxfuu,1} = [namesxf(i,:), '*', namesu(j,:), '*', namesu(k,:)];
                countxfuu = countxfuu+1; 
            end
       end
    end
end

function xinames_u_u_u_min = Get_u_u_u_min(namesu,nu)
    %unique values in kron(u,u,u)
    countuuu = 1; xinames_u_u_u_min = cell(nu*(nu+1)*(nu+2)/6,1);
    for i=1:nu
        for j=i:nu
            for k=j:nu
                xinames_u_u_u_min{countuuu,1} = [namesu(i,:), '*', namesu(j,:), '*', namesu(k,:),'-M3(',namesu(i,:), '*', namesu(j,:), '*', namesu(k,:),')'];  countuuu = countuuu+1; % unique kron(u,u,u)
            end
        end
    end
end


end % Main function end

