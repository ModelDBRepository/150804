% This code can implement different optostimulation protocols in a
% Golomb (Gol) neuron model when the 3-state model is employed to 
% account for ChR2 kinetics; the protocol is a train of ns = # number of 
% stimuli, each of width ws = 2ms, presented at a frequency f = # ;
%
% A set of previously determined parameters for ChRwt and ChRET/TC are 
% provided in comment text which must be appropriately uncomment when the 
% code is run for the chosen variant;
%
% The significance of other parameters is indicated in comments;
%
% Last update of the code: RAS 09/10/2012.

clear all; clc;

% constant parameters in Gol neuron model
global C phi g_L V_L pms pns g_Na g_NaP g_Kdr g_A g_M V_Na V_K
global teta_m sigma_m teta_p sigma_p teta_h sigma_h t_tauh teta_n sigma_n t_taun teta_a sigma_a 
global teta_b sigma_b teta_z sigma_z tau_b tau_z Idc
global Gd Gr g1


% constant parameters in Gol neuron model
C = 1;  phi = 10;
g_L = 0.05;
V_L = -70; %-62

pms = 3; pns = 4;

g_Na = 35;
g_NaP = 0.0;  % non-zero for bursting regime
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

Idc = 0.12;  % this value provides a resting potential of ~ -70 mV

%%%%%%%%%%%%%%%%% ChR2 PARAMETERS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% parameters ChRwt
Gd = 1/11.1; Gr = 1/10700; 
    l1 = 1/9.6;
    Pmax = l1+(Gr*Gd)/(l1-Gr-Gd); 
g1 = 12.6;

% % parameters ChRET/TC
% Gd = 1/8.1; Gr = 1/2600; 
%     l1 = 1/11;
%     Pmax = l1+(Gr*Gd)/(l1-Gr-Gd); 
% g1 = 12;


%%%%%%%%%%%%%% Integration Module %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% integration of the model using Euler/Runge Kutta 4
t = 0; dt = 0.05; 

% light protocol;
f = 40; % frequency (in Hz) of light stimulation
T = round(1000*(1/f)); %period of light stimulation (in ms)
TT = round(T/dt); % integrations time coresponding to the period

ws = 2;  % the width of the stimulus in ms;
tws = round(ws/dt); % integration time coresponding to each stimulus

% building the light stimulation protocol
ns = 60; % number of stimulations
in = 1000; % transient period before the optostimulation protocol
light = [zeros(1,in)]; % transient prior to the optostimulation protocol 
for ii = 1:ns
    light = [light ones(1,tws) zeros(1,TT-2*round(tws/2))];
    c1(ii) = in + (ii-1)*(tws+(TT-2*round(tws/2))) + 1;  % this is the index at the begining of each stimulation pulse
end

iters = length(light); % defining the number of integration steps

% defining the excitation rate
P = Pmax*light;

%initial conditions
    V_in = -71;%-71.962;
    h_in = 0.9; n_in = 0.02;
    b_in = 0.2; z_in = 0.001; 

    y(1) = 0; y(2) = 0;  

Vmh(1,:) = [V_in h_in n_in b_in z_in y(1) y(2)];

for ii = 1:iters-1
    
    % using RK4
    K1 = golombChR3st(t,Vmh(ii,:),P(ii));
    K2 = golombChR3st(t+dt/2,Vmh(ii,:)+dt*K1/2,P(ii)); 
    K3 = golombChR3st(t+dt/2,Vmh(ii,:)+dt*K2/2,P(ii));
    K4 = golombChR3st(t+dt,Vmh(ii,:)+dt*K3,P(ii));
    
    Vmh(ii+1,:) = Vmh(ii,:) + dt*(K1 + 2*K2 + 2*K3 + K4)/6;
   
    t(ii+1) = t(ii) + dt;
    
end

 % plot the membraine potential time series resulting from the applied
 % optostimulation protocol
 figure;
 plot(t,Vmh(:,1),'k','LineWidth',1.5);hold on;
 axis([0 1600 -90 35]);

% plot the optostimulation protocol (one rectangle for each stimulus applied)
 x_light = c1*dt; % the time (the horizontal position) of the each stimulus
 y_light = -87*ones(size(c1)); % the vertical position of each stimulus 
 w_light = 4*ones(size(c1)); % the width of the rectangle representing each stiumulus
 h_light = 10*ones(size(c1)); % the hight of the rectangle representing each stimus

  % plot of the actual train of stimuli represented by rectangles as defined
  % above and presented together with the membraine potential elicited by optostimulation
 for ii = 1:length(c1);
     rectangle('Position',[x_light(ii) y_light(ii) w_light(ii) h_light(ii)],'Facecolor','b','EdgeColor','k');hold on;
 end

xlabel('time(ms)'); ylabel('V(t)');    





