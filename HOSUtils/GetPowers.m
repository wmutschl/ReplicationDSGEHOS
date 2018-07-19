% By Willi Mutschler, September 26, 2016. Email: willi@mutschler.eu
function combos = GetPowers(A,B)
    [I,J]=ndgrid(1:size(B,1),1:size(A,1));  
    combos = [A(J(:),:) B(I(:),:)]; % matrix of powers for unique elelments
end