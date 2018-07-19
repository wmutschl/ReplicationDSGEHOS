function nXI2min = Xi_Student_t_prodmom2_orderApp1_nu1(arg)
df = arg(1);
SIGe_1 = arg(2);
nXI2min=zeros(1,1);
nXI2min(1,1) = (SIGe_1*df)/(df - 2);
