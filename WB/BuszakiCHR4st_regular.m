% This code can implement different optostimulation protocols in a
% Wang-Buszaki (WB) neuron model when the 4-state model is employed to 
% account for ChR2 kinetics; the protocol is a train of ns = # number of 
% stimuli, each of with ws = 2ms, presented at a frequency f = # ;
%
% A set of previously determined parameters for ChRwt and ChETA are 
% provided in comment text which must be appropriately uncomment when the 
% code is run for the chosen variant;
%
% The significance of other parameters is indicated in comments;
%
% Last update of the code: RAS 09/10/2012.

clear all; clc;

% constant parameters in WB neuron model
global Cm phi gNa ENa gK EK gL EL Idc
global Gd1 Gd2 Gr e12 e21 g1 gama tau_ChR


% constant parameters in WB neuron model
Cm = 1;  gNa = 35; ENa = 55; gK = 9; EK = -90; gL = 0.1; EL = -65; 
Idc = -0.51; %for a rest state around -70mV
phi = 5;

%%%%%%%%%%%%%%%%% ChR2 PARAMETERS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% parameters ChRwt model
PP1 = 0.0641; PP2 = 0.06102; Gd1 = 0.4558; Gd2 = 0.0704; e12 = 0.2044; e21 = 0.0090; 
Gr = 1/10700; gama = 0.0305;
g1 = 40;  % this parameter is variable depending on the optostimulation frequency employed
          % see the legend of Figure 5 for more details
tau_ChR =  6.3152; 

% % parameters ChETA model
% PP1 = 0.0661; PP2 = 0.0641; Gd1 = 0.0102; Gd2 = 0.1510; e12 = 10.5128; e21 = 0.0050; 
% Gr = 1/1000; gama = 0.0141; 
% g1 = 70.0; 
% tau_ChR =  1.5855;

%%%%%%%%%%%%%% Integration Module %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% integration parameters
t(1) = 0; 
dt = 0.05;

% light protocol;
f = 80; % frequency (in Hz) of light stimulation
T = round(1000*(1/f)); %period of light stimulation (in ms)
TT = round(T/dt); % integrations time coresponding to the period

ws = 2;  % the width of the stimulus in ms;
tws = round(ws/dt); % integration time coresponding to each stimulus

% building the light stimulation protocol
ns = 20; % number of stimulations
in = 1000; % transient period before the optostimulation protocol
light = [zeros(1,in)]; % transient prior to the optostimulation protocol 
for ii = 1:ns
    light = [light ones(1,tws) zeros(1,TT-2*round(tws/2))];
    c1(ii) = in + (ii-1)*(tws+(TT-2*round(tws/2))) + 1;  % this is the index at the begining of each stimulation pulse
end

iters = length(light); % defining the number of integration steps

% defining the rates of excitation
P1 = PP1*light;
P2 = PP2*light;

% initial conditions
V(1) = -80; h(1) = 0.1; n(1) = 0.01;
y(1) = 0; y(2) = 0; y(3) =0; y(4) = 0;
Vmh(1,:) = [V(1) h(1) n(1) y(1) y(2) y(3) y(4)];

% system integration
for ii = 1:iters
    
      %using RK4
    K1 = buszaki_chr(t,Vmh(ii,:),P1(ii),P2(ii));
    K2 = buszaki_chr(t+dt/2,Vmh(ii,:)+dt*K1/2,P1(ii),P2(ii)); 
    K3 = buszaki_chr(t+dt/2,Vmh(ii,:)+dt*K2/2,P1(ii),P2(ii));
    K4 = buszaki_chr(t+dt,Vmh(ii,:)+dt*K3,P1(ii),P2(ii));
         
    Vmh(ii+1,:) = Vmh(ii,:) + dt*(K1 + 2*K2 + 2*K3 + K4)/6;
    
    t(ii+1) = t(ii)+dt;
    
end

% plot the membraine potential time series resulting from the applied
% optostimulation protocol
 figure;
 plot(t,Vmh(:,1),'k','LineWidth',1.5);hold on;
 axis([20 320 -95 50]);
 
 
 % plot the optostimulation protocol (one rectangle for each stimulus applied)
 x_light = c1*dt;  % the time (the horizontal position) of the each stimulus
 y_light = -87*ones(size(c1)); % the vertical position of each stimulus 
 w_light = 4*ones(size(c1));  % the width of the rectangle representing each stiumulus
 h_light = 10*ones(size(c1)); % the hight of the rectangle representing each stimus

 % plot of the actual train of stimuli represented by rectangles as defined
 % above and presented together with the membraine potential elicited by optostimulation
 for ii = 1:length(c1);
     rectangle('Position',[x_light(ii) y_light(ii) w_light(ii) h_light(ii)],'Facecolor','b','EdgeColor','b');hold on;
 end
xlabel('time(ms)'); ylabel('V(t)');    
    
    
    



