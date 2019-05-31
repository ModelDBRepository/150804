% This code is designed to implement different optostimulation protocols in a
% Wang-Buszaki (WB) neuron model when the 4-state model is employed to 
% account for ChR2 kinetics; the protocol is a train of ns = # number of 
% stimuli, each of with ws = 2ms, presented at a frequency f = # ;
%
% A set of previously determined parameters for ChRwt and ChETA are 
% provided in comment text which must be appropriately uncomment when the 
% code is run for the chosen variant;
%
% The code is evaluating the number of extra spikes and the plateau
% potential for ChRwt and (when chosen) ChETA; runtime is about 1,30h
%
% The significance of other parameters is indicated in comments;
%
% Last update of the code: RAS 09/12/2012.

clear all; clc;

% constant parameters in WB model
global Cm phi gNa ENa gK EK gL EL Idc
global Gd1 Gd2 Gr e12 e21 g1 gama tau_ChR


% constant parameters in WB model
Cm = 1;  gNa = 35; ENa = 55; gK = 9; EK = -90; gL = 0.1; EL = -65; 
Idc = -0.51;  %for a rest state around -70mV
phi = 5;


%%%%%%%%%%%%%%%%% FUNCTIONAL PARAMETERS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% parameters ChRwt model
PP1 = 0.0641; PP2 = 0.06102; Gd1 = 0.4558; Gd2 = 0.0704; e12 = 0.2044; e21 = 0.0090; 
Gr = 1/10700; gama = 0.0305; 
g1 = 20.0; % the value of this parameter is different in other implementations 
tau_ChR = 0.5; % the value of this parameter is different in other implementations 


% % parameters ChETA model 
% PP1 = 0.0661; PP2 = 0.0641; Gd1 = 0.0102; Gd2 = 0.1510; e12 = 10.5128; e21 = 0.0050; 
% Gr = 1/1000; gama = 0.0141; 
% g1 = 70.0; 
% tau_ChR = 1.5855; 

%%%%%%%%%%%%%% Integration Module %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% integration parameters
t(1) = 0; 
dt = 0.05;

%light protocol;
ff = [5 10 20 40 60 80 100 125 150 200]; % the values of frequencies for which the number of extra spikes
                                         % and plateau potential will be evaluated
for kk = 1:length(ff)
    f = ff(kk);
T = round(1000*(1/f)); %period of light stimulation (in ms)
TT = round(T/dt); % integrations time coresponding to the period

ws = 2;  % the width of the stimulus in ms;
tws = round(ws/dt); % integration time coresponding to each stimulus

%buiding the light protocol
ns = 40; %number of stimulations
light = [zeros(1,1000)]; % transient period before the optostimulation protocol
for ii = 1:ns
    light = [light ones(1,tws) zeros(1,TT-2*round(tws/2))];
end

iters = length(light); % defining the number of integration steps

% defining the rates of excitation
P1 = PP1*light;
P2 = PP2*light;

% initial conditions
V(1) = -80; h(1) = 0.1; n(1) = 0.01;
y(1) = 0; y(2) = 0; y(3) = 0; y(4) = 0;
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

% evaluate the number of spikes 
V1 = Vmh(:,1).*(Vmh(:,1)>0); %select only the part of the time series that exceeds a certain threshold (here 0V)
VV = diff(V1); %take the derivative
ss(kk) = 0;
for jj = 1: length(VV)-1
    if (VV(jj)>0)&(VV(jj+1)<0)
        ss(kk) = ss(kk) + 1; % count the number of spikes
    end;
end

% evaluation of the plato potential
m1 = min(Vmh(1500:end-TT-500,1));
sm(kk) = abs(Vmh(1000,1) - m1);

clear  Vmh;

end

%ploting the number of extra spikes
figure;hold on;
plot(ff,(ss-ns).*((ss-ns)>=0),'k');
%plot(f,ss-ns,'r')
title('Extra Spikes'); xlabel('Light Pulse Frequency(Hz)'); ylabel('No Extra Spikes');
axis([0 200 0 120])

%ploting the plateau potential
figure;hold on;
plot(ff,sm,'k');
title('Plateau Potential');xlabel('Light Pulse Frequency(Hz)'); ylabel('Plateau Potential');
axis([0 200 0 20]);
    
    
    
    
    



