% By Willi Mutschler, September 26, 2016. Email: willi@mutschler.eu
function [D2,D3,D4] = DiagionalizationMatrix(p,second,third,fourth)
Ip = speye(p);
if second;
    D2 = spalloc(p^2,p,p);
else
    D2=[];
end
if third
    D3 = spalloc(p^3,p,p);
else
    D3=[];
end
if fourth
    D4 = spalloc(p^4,p,p);
else
    D4=[];
end
for i = 1:p
    ei = Ip(:,i);    
    eiei = kron(ei,ei);
    if third || fourth
        eieiei = kron(ei,eiei);
    end
    if fourth
        eieieiei = kron(ei,eieiei);
    end
    
    if second
        D2 = D2 + eiei*transpose(ei);
    end
    if third
        D3 = D3 + eieiei*transpose(ei);
    end
    if fourth
        D4 = D4 + eieieiei*transpose(ei);
    end
end
