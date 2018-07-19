function nXI4min = Xi_Student_t_prodmom4_orderApp1_nu1(arg)
df = arg(1);
SIGe_1 = arg(2);
nXI4min=zeros(1,1);
nXI4min(1,1) = (3*SIGe_1^2*df^2)/(4*(df/2 - 1)*(df/2 - 2));
