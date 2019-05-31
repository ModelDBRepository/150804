function Iemp = monoexpETC();

% This function provides the experimental profile for several ChR2 variants
% to be used in comparison with the theoretical solution

dt = 0.01;
t = [0:dt:2000];

ton = 500; toff = 1500;

%%%%%%%%%%%%%%%%% uncomment this section for ChRwt %%%%%%%%%%%%%%%% %%%%%%%%%%%%%%%%%%%%%%%%

% the parameters for ChRwt (Berndt)
Ipeak = 1;  % normalized peak of the measured photocurrent
Isteady = 0.27; % steady-state photocurrent
ro = Isteady/Ipeak; % steady-state to peak photocurrent ratio

tau_deact = 9.6;
tau_act = 0.23; 
tau_off = 11.1; 
tp = ton + 2.65; 

% the profile (empirical) function
Iemp = Ipeak*(1-exp(-(t-ton)/tau_act)).*heaviside(t-ton).*heaviside(tp-t) +...  % this is the rise current of the activation phase
       Ipeak*(ro + (1-ro)*exp(-(t-tp)/tau_deact)).*heaviside(t-tp).*heaviside(toff-t)+... % this is the decay current (from Ipeak to Isteady)
       Isteady*exp(-(t-toff)/tau_off).*heaviside(t-toff);   % this is the light off current
 
   for jj = 1:length(Iemp)
       if isnan(Iemp(jj)) == 1;
           Iemp(jj)=0;
       end
   end

   %%%%%%%%%%%%%%%%% uncomment this section for ChETA %%%%%%%%%%%%%%%% %%%%%%%%%%%%%%%%%%%%%%%%

% % the parameters for ChETA
% Ipeak = 1;  % normalized peak of the measured photocurrent
% Isteady = 0.6; % steady-state photocurrent
% ro = Isteady/Ipeak; % steady-state to peak photocurrent ratio
% 
% tau_deact = 15;
% tau_act = 0.08;
% tau_off = 5.2;
% tp = ton + 0.9;
% 
% %the profile (empirical) function
% Iemp = Ipeak*(1-exp(-(t-ton)/tau_act)).*heaviside(t-ton).*heaviside(tp-t) +...  % this is the rise current of the activation phase
%        Ipeak*(ro + (1-ro)*exp(-(t-tp)/tau_deact)).*heaviside(t-tp).*heaviside(toff-t)+... % this is the decay current (from Ipeak to Isteady)
%        Isteady*exp(-(t-toff)/tau_off).*heaviside(t-toff);   % this is the light off current
% 
%  
%    for jj = 1:length(Iemp)
%        if isnan(Iemp(jj)) == 1;
%            Iemp(jj)=0;
%        end
%    end

%%%%%%%%%%%%%%%%% uncomment this section for ChR ET/TC %%%%%%%%%%%%%%%% %%%%%%%%%%%%%%%%%%%%%%%%

% % the parameters for ChR ET/TC
% Ipeak = 1;  % normalized peak of the measured photocurrent
% Isteady = 0.31; % steady-state photocurrent
% ro = Isteady/Ipeak; % steady-state to peak photocurrent ratio
% 
% tau_deact = 11;
% tau_act = 0.18;
% tau_off = 8.1;
% tp = ton + 2.17;
% 
% % the profile (empirical) function
% Iemp = Ipeak*(1-exp(-(t-ton)/tau_act)).*heaviside(t-ton).*heaviside(tp-t) +...  % this is the rise current of the activation phase
%        Ipeak*(ro + (1-ro)*exp(-(t-tp)/tau_deact)).*heaviside(t-tp).*heaviside(toff-t)+... % this is the decay current (from Ipeak to Isteady)
%        Isteady*exp(-(t-toff)/tau_off).*heaviside(t-toff);   % this is the light off current
%  
%    for jj = 1:length(Iemp)
%        if isnan(Iemp(jj)) == 1;
%            Iemp(jj)=0;
%        end
%    end

