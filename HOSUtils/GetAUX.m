% By Willi Mutschler, September 26, 2016. Email: willi@mutschler.eu
function AUX = GetAUX(orderApp,nx,nu,selectHOS)
if orderApp >= 1
    nximin1 = nu;  
    AUX.XI1.Fxi = GetFxi(nu,nx,1);
    load(sprintf('DP_%i',nximin1),'DP'); 
    AUX.XI1.DPxi = DP;
    
    if selectHOS.ThirdOrder || selectHOS.FourthOrder
        load(sprintf('TP_%i',nximin1),'TP'); 
        AUX.XI1.TPxi = TP;    
    end
    
    if selectHOS.FourthOrder
        load(sprintf('QP_%i',nximin1),'QP'); 
        AUX.XI1.QPxi = QP;    
    end
end
if orderApp >= 2
    nximin2 = nu + nu*(nu+1)/2 + nu*nx;
    nameFileDuTriQuadrup2 = sprintf('DuTriQuadrup_%i',nximin2);
    AUX.XI2 = load(nameFileDuTriQuadrup2);
    AUX.XI2.Fxi = GetFxi(nu,nx,2);
end
if orderApp >= 3   
    nx2= nx*(nx+1)/2;    
    nu2 = nu*(nu+1)/2; 
    nu3 = nu2*(nu+2)/3;
    nximin3 = nu + nu2 + 2*nu*nx + nu*nx2 + nu2*nx + nu3;    
    nameFileDuTriQuadrup3 = sprintf('DuTriQuadrup_%i',nximin3);
    AUX.XI3 = load(nameFileDuTriQuadrup3);
    AUX.XI3.Fxi = GetFxi(nu,nx,3);
end

