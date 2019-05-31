% This code is designed to implement different optostimulation protocols in a
% Golomb (Gol) neuron model when the 4-state model is employed to 
% account for ChR2 kinetics; the protocol is a train of ns = # number of 
% stimuli, each of with ws = 2ms, presented at a frequency f = # ;
%
% A set of previously determined parameters for ChRwt and ChETA are 
% provided in comment text which must be appropriately uncomment when the 
% code is run for the chosen variant;
%
% The code is evaluating the firing success rate (%) for ChRwt and 
% (when chosen) ChETA; runtime is about 1,30h
%
% The significance of other parameters is indicated in comments;
%
% Last update of the code: RAS 09/12/2012.

clear all; clc;

% constant parameters in Gol neuron model
global C phi g_L V_L pms pns g_Na g_NaP g_Kdr g_A g_M V_Na V_K
global teta_m sigma_m teta_p sigma_p teta_h sigma_h t_tauh teta_n sigma_n t_taun teta_a sigma_a 
global teta_b sigma_b teta_z sigma_z tau_b tau_z Idc
global Gd1 Gd2 Gr e12 e21 g1 gama tau_ChR

% constant parameters in Gol neuron model
C = 1;  phi = 10;
g_L = 0.05;
V_L = -70; %-62

pms = 3; pns = 4;

g_Na = 35;
g_NaP = 0.0; % non-zero for bursting regime
g_Kdr = 6;
g_A = 1.4;
g_M = 1;

V_Na = 55;
V_K = -90;

teta_m = -30; sigma_m = 9.5;
teta_p = -47;
sigma_p = 3;
teta_h = -45; sigma_h = -7;
t_tauh = -40.5;
teta_n = -35; sigma_n = 10;
t_taun = -27;
teta_a = -50; sigma_a = 20;
teta_b = -80; sigma_b = -6;
teta_z = -39; sigma_z = 5;

tau_b = 15; tau_z = 75;

Idc = 0.12; % this value provides a resting potential of ~ -70 mV

%%%%%%%%%%%%%%%%% ChR2 PARAMETERS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% % parameters ChRwt model (Berndt)
% PP1 = 0.1243; PP2 = 0.0125;    % normal intensity (42 mW/mm2)
% %PP1 = 0.06807; PP2 = 0.00684; % reduced intensity (23 mW/mm2)
% %PP1 = 0.0198;   PP2 = 0.002;  % low intensity (6.7mW/mm2)
% Gd1 = 0.0105; Gd2 = 0.1181; e12 = 4.3765; e21 = 1.6046; 
% Gr = 1/10700; gama = 0.0157; %g1 = 0.0393;
% g1 = 4.9; 
% tau_ChR = 0.504; 


% parameters ChRET/TC model 
PP1 = 0.1252; PP2 = 0.0176;   % normal intensity (42mW/mm2)
%PP1 = 0.0686; PP2 = 0.00964; % reduced intensity (23 mW/mm2)
%PP1 = 0.02;   PP2 = 0.0028;  % low intensity (6.7 mW/mm2)
Gd1 = 0.0104; Gd2 = 0.1271; e12 = 16.1087; e21 = 1.090; 
Gr = 1/2600; gama = 0.0179; 
g1 = 33.0; 
tau_ChR = 0.3615; 

%%%%%%%%%%%%%% Integration Module %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% integration of the model using Euler/Runge Kutta 4
t(1) = 0; dt = 0.05; 
   
% light protocol;
ff = [3 5 10 20 40 60 80 100]; % the values of frequencies (in Hz) for which the firing success rate will be evaluated

for kk = 1:length(ff)
    f = ff(kk);
T = round(1000*(1/f)); % period of light stimulation (in ms)
TT = round(T/dt); % integrations time coresponding to the period

ws = 2;  % the width of the stimulus in ms;
tws = round(ws/dt); % integration time coresponding to each stimulus

%buiding the light protocol
ns = 40; % number of stimulations
in = 1000; % number of units before and after stimulation protocol (transient)
light = [zeros(1,in)]; % "in" units of integration are zero before the protocol starts
for ii = 1:ns
    light = [light ones(1,tws) zeros(1,TT-2*round(tws/2))];
    c1(ii) = in + (ii-1)*(tws+(TT-2*round(tws/2))) + 1;    % this is the index at the begining of each stimulation pulse
    c2(ii) = in + (ii-1)*(tws+(TT-2*round(tws/2))) + tws;  % this is the index at the end of each stimulation pulse
end

light = [light zeros(1,in)];  % "in" units of integration are zero after the protocol ends
iters = length(light);

% defining the rates of excitation
P1 = PP1*light;
P2 = PP2*light;

%initial conditions
    V_in = -71;
    h_in = 0.9; n_in = 0.02;
    b_in = 0.2; z_in = 0.001; 

    y(1) = 0; y(2) = 0; y(3) =0; y(4) = 0;

Vmh(1,:) = [V_in h_in n_in b_in z_in y(1) y(2) y(3) y(4)];

for ii = 1:iters-1
            
    % using RK4
    K1 = golombChR(t,Vmh(ii,:),P1(ii),P2(ii));
    K2 = golombChR(t+dt/2,Vmh(ii,:)+dt*K1/2,P1(ii),P2(ii)); 
    K3 = golombChR(t+dt/2,Vmh(ii,:)+dt*K2/2,P1(ii),P2(ii));
    K4 = golombChR(t+dt,Vmh(ii,:)+dt*K3,P1(ii),P2(ii));
    
    Vmh(ii+1,:) = Vmh(ii,:) + dt*(K1 + 2*K2 + 2*K3 + K4)/6;
   
    t(ii+1) = t(ii) + dt;
    
end

 
% for each interpulse interval evaluate the nuber of spikes
for ss = 1:ns-1
        % select the time series of the potential for each window of interest 
    Vs = Vmh(c2(ss):c1(ss+1),1);
 % count the number of spikes
V1 = Vs(:,1).*(Vs(:,1)>0); %select only the part of the time series that exceeds a certain threshold (here 0V)
VV = diff(V1); %take the derivative
sn(ss,kk) = 0;
for jj = 1: length(VV)-1
    if (VV(jj)>0)&(VV(jj+1)<0)
        sn(ss,kk) = sn(ss,kk) + 1; % count the number of spikes
    end;
end

end


 clear  Vmh;
 

end

% evaluation of the firing success rate 
fsr = mean(sn>0.1,1)*100;

% plot the firing success rate
figure;
plot(ff,fsr,'b'); hold on;
xlabel('Stimulation Frequency(Hz)'); ylabel('Firing Success Rate(%)');






