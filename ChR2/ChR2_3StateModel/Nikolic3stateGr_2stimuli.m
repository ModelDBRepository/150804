% This code is design to evaluate the recovery of ChRwt (modeled based on the 
% experimental data provided by Gunaydin et al.); a stimulation protocol of
% 1s(light on) + 1s(light off) + 1s(light on) is presented when the ChRwt is
% in ideal initial conditions (IIC) or special initial condition (SIC)
%
% The recovery is evaluated as the ratio r = 100*l2/l1 where l2 is the
% distance between the peak and the steady state of the photocurrent
% response to the first 1s stimulation and l1 is the same quantity
% evaluated for the photocurrent ellicited by the second 1 s
% optostimulation
%
% The code has been used to generate the time series in Figure 3 in the manuscrip
%
% To obtain the desired figure, uncomment the appropriate set of
% parameters, light stimulation protocol and plotting section 
%
% To obtain the recovery curves presented in Figure 3 panel A and C in the
% manuscript, uncomment the indicated lines in the script; the runtime for
% one curve is about 16h
%
% Last modification RAS 09/19/2012

clear all
clc

global Gd Gr

% parameters ChRwt (Gunaydin et al.)
Gd = 1/9.8; Gr = 1/10700; l1 = 1/55.5; Pmax = l1 + (Gd*Gr)/(l1-Gd-Gr) ; V = -100; g1 = 0.07; % WT Gunaydin IIC
%Gd = 1/9.8; Gr = 1/10700; l1 = 1/55.5; Pmax = l1 + (Gd*Gr)/(l1-Gd-Gr) ; V = -100; g1 = 3.687; % WT Gunaydin SIC

for jj = 1%:50 %uncomment here for the recovery curve 

    %integration parameters
    t(1) = 0; dt = 0.05;

light1 = ones(1,1000/dt); light2 = light1;
betweenl = zeros(1,jj*1000/dt);

% light input  
%P = [Pmax*light1 betweenl  Pmax*light2 zeros(1,500/dt)]; % for ideal initial conditions
P = [Pmax*light1 betweenl  Pmax*light2 zeros(1,500/dt)]; % for special initial conditions

%integration parameters
iters = length(P);

%initial conditions
%O(1) = 0.0023; D(1) = 0.9845; C(1) = 1- O(1) - D(1); % special initial conditions
O(1) = 0.0; D(1) = 0.0; C(1) = 1-O(1) - D(1); % ideal initial conditions

Out3st(1,:) = [O(1) D(1)];

% integration of the 3-state model
for ii = 1:iters-1
    
    %RG4
    K1 = Nik3stSIIC(t,Out3st(ii,:), P(ii));
    K2 = Nik3stSIIC(t+dt/2,Out3st(ii,:)+dt*K1/2, P(ii));
    K3 = Nik3stSIIC(t+dt/2,Out3st(ii,:)+dt*K2/2, P(ii));
    K4 = Nik3stSIIC(t+dt,Out3st(ii,:)+ dt*K3, P(ii));
    
    Out3st(ii+1,:) = Out3st(ii,:) + dt*(K1 + 2*K2 + 2*K3 + K4)/6;
    t(ii+1) = t(ii)+dt;
    
end

 % evaluating the photocurrent 
 I = V*g1*Out3st(:,1);
 I = [zeros(1,500/dt) I']; % adjustment necessary for a nice graph 
 t = [(-500:dt:0-dt) t];
 
 % ploting the photocurrent elicited by the optostimulation protocol
    % comment this section for the evaluation of the recovery curve
 figure;
 plot(t,I,'k');hold on;
 xlabel('time(ms)'); ylabel('Photocurrent(nA)');

 % evaluation of the distance from peak to steady state durin the two
 % periods of optostimulation
 I = -I;
 l1 = max(I(500/dt:500/dt+length(light1))) - I(500/dt+length(light1) - 10);
 l2 = max(I( 500/dt+length(light1)+length(betweenl):500/dt+length(light1)+length(betweenl)+length(light2) ) ) -...
     I(500/dt+length(light1)+length(betweenl)+length(light2) - 10);
 
 % evaluation of the recovery ratio
 r(jj) = 100*l2/l1; 

 k(jj) = jj;
 
 %clear t I Out3st; % uncomment here for the recovery curve
end
 
%%%%%%% uncomment the section below for the recovery curve %%%%%%%%%

% % plotting the recovery ratio
% figure;
% plot(k,r,'k'); hold on;
% xlabel('time(s)'); ylabel('Rec(%)');






 