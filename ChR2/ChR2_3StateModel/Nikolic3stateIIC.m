% This code is design to implement the 3-state model to account for the
% kinetics of different ChR2 variants in Ideal Initial Conditions (IIC)
% when exposed to optostimulation of 1s and 2ms;
%
% It has been used to generate the time series in Figure 1 in the manuscrip
% and Supplementary Information, panels A2,3 and B2,3
%
% To obtain the desired figure, uncomment the appropriate set of
% parameters, light stimulation protocol and plotting section 

clear all
clc

global Gd Gr

% parameters ChR2 variant 
Gd = 1/9.8; Gr = 1/10700; l1 = 1/55.5; Pmax = l1 + (Gd*Gr)/(l1-Gd-Gr) ; V = -100; g1 = 0.07; % WT Gunaydin
%Gd = 1/5.2; Gr = 1/1000; l1 =1/15; Pmax = l1 + (Gd*Gr)/(l1-Gd-Gr) ; V = -100; g1 = 0.03314;  % ETA Gunaydin

%Gd = 1/11.1; Gr = 1/10700; l1 = 1/9.6; Pmax = l1 + (Gd*Gr)/(l1-Gd-Gr) ; V = -75; g1 = 0.03256; % WT Berndt
%Gd = 1/8.1; Gr = 1/2600; l1 = 1/11; Pmax = l1 + (Gd*Gr)/(l1-Gd-Gr) ; V = -75; g1 = 0.06097; % ETC Berndt
           
%integration parameters
dt = 0.05;
t(1) = 0;

% light stimulation protocol
P = [zeros(1,500/dt) Pmax*ones(1,1000/dt) zeros(1,500/dt)]; % for 1s ( = 1000ms) optostimulation
%P = [zeros(1,500/dt) Pmax*ones(1,2/dt) zeros(1,500/dt)]; % for 2ms optostimulation

% evaluation of the number of integrations steps
iters = length(P);
         
% ideal initial conditions
O(1) = 0.0; D(1) = 0.0; C(1) = 1-O(1) - D(1); %ideal conditions

Out3st(1,:) = [O(1) D(1)];

% integration of the 3-state model
for ii = 1:iters-1
                   
    %RG4
    K1 = Nik3stSIIC(t,Out3st(ii,:), P(ii));
    K2 = Nik3stSIIC(t+dt/2,Out3st(ii,:)+dt*K1/2, P(ii));
    K3 = Nik3stSIIC(t+dt/2,Out3st(ii,:)+dt*K2/2, P(ii));
    K4 = Nik3stSIIC(t+dt,Out3st(ii,:)+ dt*K3, P(ii));
    
    Out3st(ii+1,:) = Out3st(ii,:) + dt*(K1 + 2*K2 + 2*K3 + K4)/6;
    
end

I = V*g1*Out3st(:,1);

t = 0:dt:(length(I)-1)*dt;

%%%%%%%%%%%%%%%% plotting section %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

figure;  % uncomment this section for optostimulation of 1s
plot(t,I,'r');hold on;
axis([0 2000 -0.95 0.1]); % adjust the limits if necessary when different variants are used
xlabel('time(ms)'); ylabel('Photocurrent (nA)');


% figure;  % uncomment this section for optostimulation of 2ms
% plot(t,I./(-min(I)),'r'); hold on;
% axis([450 600 -1.02 0.1]); % adjust the limits if necessary when different variants are used
% xlabel('time(ms)'); ylabel('Normalized Photocurrent (nA)');



