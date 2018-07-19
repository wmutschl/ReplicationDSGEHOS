% By Willi Mutschler, September 26, 2016. Email: willi@mutschler.eu
% This function 
% (1) Computes unique elements of the product moments of xi_t symbolically
% (2) Evaluates the unique elements of the product moments of xi_t for a Gaussian or t-distribution and writes them to files


function XImin = ProdMom_inov(xikmin,nu,nx,Mom,distrib)
nu2 = nu*(nu+1)/2; 
nu3 = nu*(nu+1)*(nu+2)/6;

STR = char(xikmin); % get individual entry as string
% find summands: how many, + or -
Sumsplus = strfind(STR,'+'); Sumsminus = strfind(STR,'-'); 
Sums_idx1 = sort([1 (Sumsplus+2) (Sumsminus+2)]);
Sums_idx2 = sort([(Sumsplus-1) (Sumsminus-2) length(STR)]);
numSummands = numel(Sumsplus) + numel(Sumsminus)+1;        
m_sum = sym(0); %initialize
for iSum = 1:numSummands      % look at each summand individually      
    indx_mom_u = zeros(1,nu); % initialize index that determines powers of E(e_1^{i_1}... e_{n_e}^{i_{n_e})
    indx_mom_xf = zeros(1,nx);% initialize index that determines powers of E(xf_1^{i_1}... xf_{n_x}^{i_{n_x} xs_1^{j_1}... xs_{n_x}^{j_{n_x})
    indx_mom_xs = zeros(1,nx);% initialize index that determines powers of E(xs_1^{i_1}... xs_{n_x}^{i_{n_x} xs_1^{j_1}... xs_{n_x}^{j_{n_x})
    indx_uu = zeros(1,nu2);
    indx_uuu = zeros(1,nu3);
    Str = STR(Sums_idx1(iSum):Sums_idx2(iSum));                
    % Find indx_mom_e for powers of epsi recursively
    for iu=nu:-1:1
        str=['u',num2str(iu), '^'];
        check = strfind(Str,str);
        if ~isempty(check)
            indx_mom_u(iu) = str2double(Str(check + length(str)));                
            Str = strrep(Str,['u',num2str(iu), '^',Str(check + length(str))],'');
        else
            str(end) = [];
            check = strfind(Str,str);
            if ~isempty(check)
                indx_mom_u(iu) = 1;
                Str = strrep(Str,['u',num2str(iu)],'');
            end        
        end
    end
    % Find indx for constants of second moments of epsi
    for iuu=nu2:-1:1
        str=['E_uu_',num2str(iuu), '^'];
        check = strfind(Str,str);
        if ~isempty(check)
            indx_uu(iuu) = str2double(Str(check + length(str)));
            Str = strrep(Str,['E_uu_',num2str(iuu), '^',Str(check + length(str))],'');
        else
            str(end) = [];
            check = strfind(Str,str);
            if ~isempty(check)
                indx_uu(iuu) = 1;
                Str = strrep(Str,['E_uu_',num2str(iuu)],'');
            end        
        end        
    end
    % Find indx for constants of third moments of epsi
    for iuuu=nu3:-1:1
        str=['E_uuu_',num2str(iuuu), '^'];
        check = strfind(Str,str);
        if ~isempty(check)
            indx_uuu(iuuu) = str2double(Str(check + length(str)));
            Str = strrep(Str,['E_uuu_',num2str(iuuu), '^',Str(check + length(str))],'');
        else
            str(end) = [];
            check = strfind(Str,str);
            if ~isempty(check)
                indx_uuu(iuuu) = 1;
                Str = strrep(Str,['E_uuu_',num2str(iuuu)],'');
            end        
        end        
    end
    % Find indx_mom_xf for xf
    for ix=nx:-1:1
        str=['xf',num2str(ix), '^'];
        check = strfind(Str,str);
        if ~isempty(check)
            indx_mom_xf(ix) = str2double(Str(check + length(str)));
            Str = strrep(Str,['xf',num2str(ix), '^',Str(check + length(str))],'');
        else
            str(end) = [];
            check = strfind(Str,str);
            if ~isempty(check)
                indx_mom_xf(ix) = 1;
                Str = strrep(Str,['xf',num2str(ix)],'');
            end        
        end        
    end
    % Find indx_mom_xs for xs
    for ix=nx:-1:1
        str=['xs',num2str(ix), '^'];
        check = strfind(Str,str);
        if ~isempty(check)
            indx_mom_xs(ix) = str2double(Str(check + length(str)));
            Str = strrep(Str,['xs',num2str(ix), '^',Str(check + length(str))],'');
        else
            str(end) = [];
            check = strfind(Str,str);
            if ~isempty(check)
                indx_mom_xs(ix) = 1;
                Str = strrep(Str,['xs',num2str(ix)],'');
            end        
        end        
    end
        
    % Find other constants, ie numbers
    Str = strrep(Str,'*','');
    const = strrep(Str,' ','');
    if isempty(str2double(const)) || isnan(str2double(const))
        const = 1;
    else
        const = str2double(const);
    end

    % Evaluate indices
    % Evaluate index for product moments of epsi
    if sum(indx_mom_u)==0
        mom_u=1;
    else
        switch distrib
            case 'Gaussian'
                mom_u = prodmom(Mom.SIGe,1:nu,indx_mom_u);
            case 'Student_t'
                mom_u = mom_invGamma(Mom.df,sum(indx_mom_u(1:nu))/2)*prodmom(Mom.SIGe,1:nu,indx_mom_u);            
        end
    end
    
    % Evaluate index for powers of terms belonging to E_uu
    if sum(indx_uu)==0
        mom_u2min=1; % no terms belonging to powers of E_uu
    else
        mom_u2min=prod(power(Mom.E_uu,indx_uu'));
    end
    
    % Evaluate index for powers of terms belonging to E_uuu
    if sum(indx_uuu)==0
        mom_u3min=1; % no terms belonging to powers of E_uuu
    else
        mom_u3min=prod(power(Mom.E_uuu,indx_uuu'));
    end

   
    % Evaluate index for product moments of xf and xs
    sum_idxf = sum(indx_mom_xf); sum_idxs = sum(indx_mom_xs); 
    indx_mom_xfs = [indx_mom_xf indx_mom_xs]; sum_idxfs = sum_idxf + sum_idxs;
    if sum_idxfs == 0
        mom_xfs=1; % no terms belonging to powers of xf or xs
    elseif sum_idxfs == 1
        if sum_idxf == 1
            mom_xfs = Mom.E_XF{sum_idxf}(find(indx_mom_xf)); % Mean of xf
        elseif sum_idxs == 1
            mom_xfs = Mom.E_XS{sum_idxs}(find(indx_mom_xs)); % Mean of xs
        end
    else
        if sum_idxs==0 %Only powers of xf
            E_XF = Mom.E_XF{sum_idxf}; combo = Mom.combos.E_XF{sum_idxf};
            mom_xfs = E_XF(find(ismember(combo,indx_mom_xf,'rows')));                    
        elseif sum_idxf==0 %Only powers of xs
            E_XS = Mom.E_XS{sum_idxs}; combo = Mom.combos.E_XS{sum_idxs};
            mom_xfs = E_XS(find(ismember(combo,indx_mom_xs,'rows')));
        else
            E_XF_XS = Mom.E_XF_XS{sum_idxf,sum_idxs}; combo = Mom.combos.E_XF_XS{sum_idxf,sum_idxs};
            mom_xfs = E_XF_XS(find(ismember(combo,indx_mom_xfs,'rows')));                    
        end
    end

    
    % Evaluate signs
    if iSum == 1 && strcmp(STR(1),'-')
        sumsign = -1;
    elseif iSum > 1 && strcmp(STR(Sums_idx1(iSum)-2),'-')
        sumsign = -1;
    else
        sumsign = 1;
    end
    % Put everything together    
    m_sum = m_sum + sumsign*const*mom_u*mom_u2min*mom_u3min*mom_xfs;        
end % for end                
XImin = simplify(m_sum); % Simplify symbolic expression





end % main function end
