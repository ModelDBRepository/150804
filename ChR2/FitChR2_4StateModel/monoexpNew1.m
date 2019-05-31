function Iemp = monoexpNew1(t, ton, toff);

% this function is designed to provide the empirical profile function of
% the photocurrent elicited by 1s optostimulation based on the experimental 
% data presented in Gunaydin et al. 2010

%%%%%%%%%%%%%%%%% uncomment this section for ChRwt %%%%%%%%%%%%%%%% %%%%%%%%%%%%%%%%%%%%%%%%

% % the parameters for ChR2 WT
% Ipeak = 0.848; % peak of the measured photocurrent
% ro = 0.4; % steady-state to peak photocurrent ratio
% Isteady = ro*Ipeak; % steady-state photocurrent
% 
% tau_deact = 55.5; % the time constant for the decay from peak to steady-state
% tau_act = 0.2; % the time constant of the photocurrent raise from zero to peak
% tau_off = 9.8; % the time constant of the photocurrent decay from steady state to zero
% tp = ton + 2.4; % the time necessary for the photocurrent to reach the peak
% 
% % the profile (empirical) function 
% Iemp = Ipeak*(1-exp(-(t-ton)/tau_act)).*heaviside(t-ton).*heaviside(tp-t) +...  % this is the rise current of the activation phase
%        Ipeak*(ro + (1-ro)*exp(-(t-tp)/tau_deact)).*heaviside(t-tp).*heaviside(toff-t)+... % this is the decay current (from Ipeak to Isteady)
%        Isteady*exp(-(t-toff)/tau_off).*heaviside(t-toff);   % this is the light off current
% 
% % in some older versions of Matlab, Heaviside function can generate a Nan value at the point of 
% % inflection which may induce failure when search/optimization algorithms are employed
% % to eliminate this issue we set these Nan values to zero
%    for jj = 1:length(Iemp)
%        if isnan(Iemp(jj)) == 1;
%            Iemp(jj)=0;
%        end
%    end
 

   %%%%%%%%%%%%%%%%% uncomment this section for ChETA %%%%%%%%%%%%%%%% %%%%%%%%%%%%%%%%%%%%%%%%
   
% the parameters for ChETA (Gunaydin et al.)
Ipeak = 0.645; % peak of the measured photocurrent
ro = 0.6; % steady-state to peak photocurrent ratio
Isteady = ro*Ipeak; % steady-state photocurrent

tau_deact = 15;  % the time constant for the decay from peak to steady-state
tau_act = 0.08;  % the time constant of the photocurrent raise from zero to peak
tau_off = 5.2;  % the time constant of the photocurrent decay from steady state to zero
tp = ton + 0.9; % the time necessary for the photocurrent to reach the peak

%the profile (empirical) function 
Iemp = Ipeak*(1-exp(-(t-ton)/tau_act)).*heaviside(t-ton).*heaviside(tp-t) +...  % this is the rise current of the activation phase
       Ipeak*(ro + (1-ro)*exp(-(t-tp)/tau_deact)).*heaviside(t-tp).*heaviside(toff-t)+... % this is the decay current (from Ipeak to Isteady)
       Isteady*exp(-(t-toff)/tau_off).*heaviside(t-toff);   % this is the light off current

% in some older versions of Matlab, Heaviside function can generate a Nan value at the point of 
% inflection which may induce failure when search/optimization algorithms are employed
% to eliminate this issue we set these Nan values to zero 
   for jj = 1:length(Iemp)
       if isnan(Iemp(jj)) == 1;
           Iemp(jj)=0;
       end
   end





