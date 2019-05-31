% This code can implement different optostimulation protocols in a
% Wang-Buszaki (WB) neuron model when the 3-state model is employed to 
% account for ChR2 kinetics; the protocol is a train of ns = # number of 
% stimuly, each of with ws = 2ms, presented at a frequency f = # ;
%
% A set of previously determined parameters for ChRwt and ChETA are 
% provided in comment text which must be appropriately uncomment when the 
% code is run for the chosen variant;
%
% Other parameter which may take different values (depending on the addressed issue)
% are indicated in comments;
%
% Last update of the code: RAS 09/10/2012.

clear all; clc;

% constant parameters in WB neuron model
global Cm phi gNa ENa gK EK gL EL Idc
global Gr Gd g1 


% other parameters in WB neuron model
Cm = 1;  gNa = 35; ENa = 55; gK = 9; EK = -90; gL = 0.1; EL = -65; 
Idc = -0.51;  %for a rest state around -70mV
phi = 5;

%%%%%%%%%%%%%%%%% ChR2 PARAMETERS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% % parameters ChRwt model
% Gd = 1/9.8; Gr = 1/10700; 
%     l1 = 1/55.5; 
%     Pmax = l1+(Gr*Gd)/(l1-Gr-Gd); 
% g1 = 4.8; % this parameter is variable depending on the optostimulation frequency employed
%           % see the legend of Figure 5 for more details

% parameters ChETA model
Gd = 1/5.2; Gr = 1/1000; 
    l1 = 1/15;
    Pmax = l1+(Gr*Gd)/(l1-Gr-Gd); 
g1 = 2.2; % this parameter is variable depending on the optostimulation frequency employed
          % see the legend of Figure 5 for more details

%%%%%%%%%%%%%% Integration Module %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% integration parameters
t(1) = 0; 
dt = 0.05;


% light protocol;
f = 80; % frequency (in Hz) of light stimulation
T = round(1000*(1/f)); %period of light stimulation (in ms)
TT = round(T/dt); % integrations time coresponding to the period

ws = 2;  % the width of each stimulus in ms;
tws = round(ws/dt); % integration time coresponding to each stimulus  

% buiding the light stimulation protocol
ns = 20; % number of stimulations
in = 1000; % transient period before the optostimulation protocol
light = [zeros(1,in)]; % transient prior to the optostimulation protocol 
for ii = 1:ns
    light = [light ones(1,tws) zeros(1,TT-2*round(tws/2))];
    c1(ii) = in + (ii-1)*(tws+(TT-2*round(tws/2))) + 1;  % this is the index at the begining of each stimulation pulse
end

iters = length(light); % defining the number of integration steps

% defining the rate of excitation
P = Pmax*light;


% initial conditions
V(1) = -80; h(1) = 0.1; n(1) = 0.01;
y(1) = 0; y(2) = 0; 
Vmh(1,:) = [V(1) h(1) n(1) y(1) y(2)];

% system integration
for ii = 1:iters
    
      %using RK4
    K1 = buszaki_chr3st(t,Vmh(ii,:),P(ii));
    K2 = buszaki_chr3st(t+dt/2,Vmh(ii,:)+dt*K1/2,P(ii)); 
    K3 = buszaki_chr3st(t+dt/2,Vmh(ii,:)+dt*K2/2,P(ii));
    K4 = buszaki_chr3st(t+dt,Vmh(ii,:)+dt*K3,P(ii));
         
    Vmh(ii+1,:) = Vmh(ii,:) + dt*(K1 + 2*K2 + 2*K3 + K4)/6;
    
    t(ii+1) = t(ii)+dt;
    
end

% plot the membraine potential time series resulting from the applied
% optostimulation protocol
 figure;
 plot(t,Vmh(:,1),'r','LineWidth',1.5);hold on;
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

    
    
    



