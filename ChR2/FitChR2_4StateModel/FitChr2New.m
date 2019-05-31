%function E = FitChr2New(x)
%function [E Inik1 Iemp1 Inik2 Iemp2] = FitChr2New(x)
%
% This function is designed to evaluate the ChR2 photocurrent for 1s and
% 2ms optostimulation and compare it with the empirical profile developed
% based on the experimental data provided by Gunaydin et al. 
%
% It can be called by searching functions in Matlab to optimize the
% parameters or it can be called to provide the photocurrents time series
% and the error relative to their profiles.
%
% The function calls multiple other functions to integrate the 4-state
% system accounting for ChR2 kinetics; the integration is performed with
% ode45 with a time step dt=0.01 in Matlab; the stiffness of the problem is 
% not significant; the use of ode15s instead of ode45 leads to the same results. 
% The choice of this integrations function is also motivated by the
% significantly shorter time of integration than Runge Kutta 4 provides.
%
% In comments are provided the set of parameteres found by optimization
% paradigms for ChRwt and the ChETA variants investigated in the manuscript;
%
% To obtain a figure like Figure 1, panel A4 from the manuscript,
% comment the very first line, uncomment the appropriate parameter vector x, 
% the code lines for figure ploting and make sure the appropriate code is 
% chosen in the functions called by this script.  
%
% A similar code has been used for ChRwt and the fast variant ChR ET/TC
% following the experimental data presented in Berndt et al.2011
%
% Last modification RAS 09/19/2012

 clear all; clc;
   
% % parameters parameters for ChRwt (Gunaydin et al.)
% x = [0.0641    0.06102    0.4558    0.0704    0.2044    0.0090   11.3589 0.0305    6.3152];  % final for paper
% Gr = 1/10700;

% parameters for ChETA (Gunaydin et al.)
x = [0.0661    0.0641    0.0102    0.1510   10.5128    0.0050   87.5898   0.0141    1.5855]; % final for paper
Gr = 1/1000;
  
%%% VERY IMPORTANT: the parameter g1 hear is actually g1*V 
   
% the parameters in the 4-state model are:
P1 = x(1); P2 = x(2);
Gd1 = x(3); Gd2 = x(4); e12 = x(5); e21 = x(6); 
g1 = x(7);  gama = x(8); tau_ChR = x(9);


% genertating the vector to be passed to other functions called below
PP = [P1 P2 Gd1 Gd2 e12 e21 g1 gama tau_ChR Gr]; 
 
% first condition: the photocurrent must match the profile for long time
% stimulation (1s) 
dt1 = 0.01; % integration time step
time1 = 0:dt1:2000; % total time of integration (1.4s = 1400ms)
ton1= 500; toff1 = 1500; % begining and end of optostimulation (total optostimulation time 1s = 1000ms)
 
Inik1 = Nikolic4stFitNew(time1, ton1, toff1, PP); % evaluate the photocurrent generated  
                                                  % by the 4-state model (set of differentinal equation) 
Iemp1 = monoexpNew1(time1, ton1, toff1); % evaluate the empirical photocurrent profile

% second condition: the photocurrent must match the profile for short time
% stimulation (2ms)
dt2 = 0.01;  % integration time step
time2 = 0:dt2:1000; % total time of integration (1s = 1000ms)
ton2 = 500; toff2 = 502; % begining and end of optostimulation (total optostimulation time 2ms)

Inik2 = Nikolic4stFitNew2ms(time2, ton2, toff2, PP); % evaluate the photocurrent generated  
                                                     % by the 4-state model (set of differentinal equation) 
Iemp2 = monoexpNew2(time2, ton2, toff2 ); % evaluate the empirical photocurrent profile


E1 = sqrt(mean((Inik1' - Iemp1).^2)); % RMSD for 1s optostimulation
E2 = sqrt(mean((Inik2' - Iemp2).^2)); % RMSD for 2 ms
E3 = abs( max(Inik1)-max(Iemp1) )/max(Iemp1); % PE for the photocurrent peak
E = 100*(E1+E2+E3); % total error
 
figure; hold on;
subplot(1,2,1); hold on;plot(time1,-Inik1,'k'); plot(time1,-Iemp1,'b');box on;
    xlabel('time(ms)'); ylabel('Photocurrent(nA)');
subplot(1,2,2); hold on;plot(time2,-Inik2,'k'); plot(time2,-Iemp2,'b');axis([450 600 -1.02 0.1]); box on;
    xlabel('time(ms)'); ylabel('Normalized Photocurrent(nA)');




