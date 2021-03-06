% this function is integrating the 4-state model for 1s optostimulation;
% note that the photocurrent is not normalized

function Inik = Nikolic4stFitNew(t, ton, toff, PP)

g1 = PP(7); gama = PP(8);


% integrating with ode45

[T,y] = ode45(@Nik4stfitodeNew, t, [0 0 0 0], [], ton, toff,  PP);

 I = g1*(y(:,1)+gama*y(:,2));

Inik = I;
%Inik = I./max(I);
