% This function is integrating the Golomb neural model with ChR2 kinetics 
% implemented by the 4-state model
%
% Last update: RAS 09/10/2012

function Vmh_dot = golombChR(t,Vmh,P1,P2);

global C phi g_L V_L pms pns g_Na g_NaP g_Kdr g_A g_M V_Na V_K
global teta_m sigma_m teta_p sigma_p teta_h sigma_h t_tauh teta_n sigma_n t_taun teta_a sigma_a
global teta_b sigma_b teta_z sigma_z tau_b tau_z Idc
global Gd1 Gd2 Gr e12 e21 g1 gama tau_ChR


% recovering the variable
V = Vmh(1); 
h = Vmh(2); n = Vmh(3);
b = Vmh(4); z = Vmh(5);
y(1) = Vmh(6); y(2) = Vmh(7); y(3) = Vmh(8); y(4) = Vmh(9);


h_dot = phi*(gammafn(V,teta_h,sigma_h)-h)/(1.0+7.5*gammafn(V,t_tauh,-6.0));
n_dot = phi*(gammafn(V,teta_n,sigma_n)-n)/(1.0+5.0*gammafn(V,t_taun,-15.0));
b_dot = (gammafn(V,teta_b,sigma_b)-b)/tau_b;
z_dot = (gammafn(V,teta_z,sigma_z)-z)/tau_z;


    m_inf = gammafn(V,teta_m,sigma_m);
    p_inf = gammafn(V,teta_p,sigma_p);
    a_inf = gammafn(V,teta_a,sigma_a);
    
            I_Na = g_Na*(m_inf^pms)*h*(V - V_Na);
            I_NaP = g_NaP*p_inf*(V - V_Na);

            I_Kdr = g_Kdr*(n^pns)*(V - V_K);
            I_A   = g_A*(a_inf^3)*b*(V - V_K);
            I_M   = g_M*z*(V - V_K);
            
            I_L = g_L*(V - V_L);
            
            %IChR = 0;
            IChR = V*g1*(y(1)+gama*y(2));
                        
            Itot = -I_L - I_Na - I_NaP - I_Kdr - I_A - I_M - IChR + Idc; 
            
S0 = 0.5*(1+tanh(120*((P1>0) - 0.1))); 

 dy(1) = P1*y(4)*(1-y(1)-y(2)-y(3))-(Gd1+e12)*y(1) + e21*y(2);
 dy(2) = P2*y(4)*y(3) + e12*y(1) - (Gd2+e21)*y(2);
 dy(3) = Gd2*y(2) - (P2*y(4)+Gr)*y(3);
 dy(4) = (S0 - y(4))/tau_ChR;
            

V_dot = Itot/C;
        
Vmh_dot = [V_dot h_dot n_dot b_dot z_dot dy(1) dy(2) dy(3) dy(4)]; 




