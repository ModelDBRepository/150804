% This function is designed to integrate the WB model with ChR2 kinetics
% provided by the 3-state model;
%
% Last update: RAS 09/10/2012

function Vmdot = buszaki_chr3st(t,Vmh,P);

global Cm phi gNa ENa gK EK gL EL Idc
global Gd Gr g1


% neural functional parameters
V = Vmh(1); h = Vmh(2); n = Vmh(3);
y(1) = Vmh(4); y(2) = Vmh(5);

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
                    IChR = V*g1*y(1);
                    Itot = -INa - IK - IL - IChR + Idc;


% integration
hdot = phi*(ah*(1-h) - bh*h);
ndot = phi*(an*(1-n) - bn*n);

% 3-state model
 dy(1) = P*(1-y(1)-y(2)) - Gd*y(1);
 dy(2) = Gd*y(1) - Gr*y(2);
 
 
 Vdot = (Itot/Cm);

 Vmdot = [Vdot hdot ndot dy(1) dy(2)];
 
 
 
 
 
 
 
 
 
 
 