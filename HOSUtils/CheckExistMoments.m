% By Willi Mutschler, September 26, 2016. Email: willi@mutschler.eu
function CheckExistMoments(orderApp,HOSThirdOrder,HOSFourthOrder,distrib,M_)
if HOSFourthOrder
    minorderMom(1) = 4; minorderMom(2) = 8; minorderMom(3) = 12;
elseif HOSThirdOrder
    minorderMom(1) = 3; minorderMom(2) = 6; minorderMom(3) = 9;
else
    minorderMom(1) = 2; minorderMom(2) = 4; minorderMom(3) = 6;
end

switch distrib
    case 'Student_t'
        for i = 1:M_.param_nbr
            if strcmp(deblank(M_.param_names(i,:)),'df_studt')
                if M_.params(i) <= minorderMom(orderApp)           
                    error('To compute specified moments for a %d-order approximation, one requires at least finite %d moments of the shocks. Increase degrees of freedoms of t distribution.',orderApp,minorderMom(orderApp));
                end                
            end
        end        
end
