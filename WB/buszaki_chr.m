% This function is designed to integrate the WB model with ChR2 kinetics
% provided by the 4-state model;
%
% Last update: RAS 09/10/2012

function Vmdot = buszaki_chr(t,Vmh,P1,P2);

global Cm phi gNa ENa gK EK gL EL Idc
global Gd1 Gd2 Gr e12 e21 g1 gama tau_ChR


% WB model with ChR2
V = Vmh(1); h = Vmh(2); n = Vmh(3);
y(1) = Vmh(4); y(2) = Vmh(5); y(3) = Vmh(6); y(4) = Vmh(7);

            am = -0.1*(V+35)/(exp(-0.1*(V+35)) - 1);
            ah = 0.07*exp(-(V+58)/20);
            an = -0.01*(V+34)/(exp(-0.1*(V+34)) - 1);

            bm = 4*exp(-(V+60)/18);
            bh = 1/(exp(-0.1*(V+28)) + 1);
            bn = 0.125*exp(-(V+44)/80);

            minf = am/(am+bm);

                    % evaluating the currents
                    INa  = gNa*(minf^3)*h*(V-ENa);
                    IK   = gK*(n^4)*(V-EK);
                    IL   = gL*(V-EL);
                    %IChR = 0;
                    IChR = V*g1*(y(1)+gama*y(2));
                    Itot = -INa - IK - IL - IChR + Idc;


% integration
hdot = phi*(ah*(1-h) - bh*h);
ndot = phi*(an*(1-n) - bn*n);

% 4-state model
S0 = 0.5*(1+tanh(120*((P1>0) - 0.1))); 

 dy(1) = P1*y(4)*(1-y(1)-y(2)-y(3))-(Gd1+e12)*y(1) + e21*y(2);
 dy(2) = P2*y(4)*y(3) + e12*y(1) - (Gd2+e21)*y(2);
 dy(3) = Gd2*y(2) - (P2*y(4)+Gr)*y(3);
 dy(4) = (S0 - y(4))/tau_ChR;
 
 
 Vdot = (Itot/Cm);

 Vmdot = [Vdot hdot ndot dy(1) dy(2) dy(3) dy(4)];
 
 
 
 
 
 
 
 
 
 
 