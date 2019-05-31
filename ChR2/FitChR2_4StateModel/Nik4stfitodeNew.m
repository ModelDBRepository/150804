% this function defines the 4-state model

function dy = Nik4stfitodeNew(t, y, ton, toff, PP);

tau_ChR = PP(9);

Gd1 = PP(3); Gd2 = PP(4); e12 = PP(5); e21 = PP(6); 
Gr = PP(10); 

S0 = 0.5*(1+tanh(120*(heaviside(t-ton).*heaviside(toff-t) - 0.1))); 

dy = zeros(4,1);

 dy(1) = PP(1)*y(4).*(1-y(1)-y(2)-y(3))-(Gd1+e12)*y(1) + e21*y(2);
 dy(2) = PP(2)*y(4).*y(3) + e12*y(1) - (Gd2+e21)*y(2);
 dy(3) = Gd2*y(2) - (PP(2).*y(4)+Gr)*y(3);
 dy(4) = (S0 - y(4))/tau_ChR;



