% This code is dedicated to find an analythical solution to the 4 state
% model for different ChR2 variants; 
%
% The code implements the steps described in the Supplementary Information
% to evaluate the photocurrent for Light On condition and uses the
% theoretical solution provided in Nikolic et al. 2009 to evaluate the
% photocurrent for Light Off condition
%
% For Light On condition, a visual comparison between the normalized 
% photocurrent evaluated using theoretical solution and the empirical 
% profile can be displayed for the ChR2 variant of choice modeled based on 
% the experimental data
%
% The slight missmatch between the two photocurrents is mainly due to the
% small error introduced by the optimizations parametric search employed to
% find the 4-state model parameters
%
% Last modification RAS 09/18/2012


% declaration of symbolic variables
clear all; clc;
syms P1 P2 Gd1 Gd2 e12 e21 Gr x lamda1 lamda2 lamda3 t a b c d 
syms A B Z P yp yc y0
syms C1 C2 C3 y10 y20 y30 
syms A11  A12 A13 A21 A22 A23 A31 A32 A33 yp1 yp2 yp3


%%%%%%%%%%%%%%%%%%%%%%  LIGHT ON CONDITION  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% the eigenvalues are given by the characteristic (cubic) equation:
xsol = solve('-(P1+Gd1+e12+x)*(Gd2+e21+x)*(P2+Gr+x) - P1*e12*Gd2 + P2*Gd2*(P1+Gd1+e12+x) + e12*(e21-P1)*(P2+Gr+x)=0',...
'IgnoreAnalyticConstraints',false,x);   

% the solutions are:
lamda1 = xsol(1);
lamda2 = xsol(2);
lamda3 = xsol(3);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  ChR2 PARAMETERS   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% parameters ChRwt (Berndt)
P1 = 0.1243; P2 = 0.0125; Gd1 = 0.0105 ; Gd2 = 0.1181; e12 = 4.3765; e21 = 1.6046; 
Gr = 1/10700; gama = 0.0157; 
g1 = 0.098; 

% % parameters ChETA (Gunaydin)
% P1 = 0.0661; P2 = 0.0641; Gd1 = 0.0102 ; Gd2 = 0.1510; e12 = 10.5128; e21 = 0.0050; 
% Gr = 1/1000; gama = 0.0141; 
% g1 = 0.8759;  

% % parameters ChR model ETC as measured in Berndt 2011 (fitting Nikolic 4 state model)
% P1 = 0.1252; P2 = 0.0176; Gd1 = 0.0104 ; Gd2 = 0.1271; e12 = 16.1087; e21 = 1.0900; 
% Gr = 1/2600; gama = 0.0179; 
% g1 = 0.5599;  

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% although the eigenvalues are real, Matlab provides them in a complex
% format with a null imaginary component; we must keep only the real part
lamda1 = real(eval(lamda1)); lamda2 = real(eval(lamda2)); lamda3 = real(eval(lamda3));

% simplifications of notations=
a = P2+Gr; 
b = P1+Gd1+e12; 
c = e21-P1;
d = P1*Gd2;

% the complementary matrix
A = [exp(lamda1*t) (a+lamda1)*(b+lamda1)*exp(lamda1*t)/(c*(a+lamda1)-d) Gd2*(b+lamda1)*exp(lamda1*t)/(c*(a+lamda1)-d) ;...
     exp(lamda2*t) (a+lamda2)*(b+lamda2)*exp(lamda2*t)/(c*(a+lamda2)-d) Gd2*(b+lamda2)*exp(lamda2*t)/(c*(a+lamda2)-d) ;...
     exp(lamda3*t) (a+lamda3)*(b+lamda3)*exp(lamda3*t)/(c*(a+lamda3)-d) Gd2*(b+lamda3)*exp(lamda3*t)/(c*(a+lamda3)-d) ]';

% evaluation of the inverse of the complementary matrix 
B = A^(-1);

% evaluation of the particular solution
Z = B*[P1;0;0]; P = int(Z,t);

yp = A*P;
yc = A*[C1;C2;C3];
y0 = [y10;y20;y30];

t=0;
yp = real(eval(yp)); yc = real(eval(yc)); 

% evaluation of the coefficients using the initial conditions
AA = real(eval(A)); 
A11 = AA(1,1); A12 = AA(1,2); A13 = AA(1,3);
A21 = AA(2,1); A22 = AA(2,2); A23 = AA(2,3);
A31 = AA(3,1); A32 = AA(3,2); A33 = AA(3,3);

yp1 = yp(1); yp2 = yp(2); yp3 = yp(3);

C1 = (-yp1-A12*C2-A13*C3)/A11;
C2 = (A21*yp1-A11*yp2-C3*(A23*A11-A21*A13))/(A22*A11-A21*A12);
s = solve('A31*(-yp1-A12*(A21*yp1-A11*yp2-C3*(A23*A11-A21*A13))/(A22*A11-A21*A12)-A13*C3)/A11 + A32*(A21*yp1-A11*yp2-C3*(A23*A11-A21*A13))/(A22*A11-A21*A12) + A33*C3 + yp3=0',...
    'IgnoreAnalyticConstraints',false,C3);

C3 = real(eval(s));
C2 = real(eval(C2));
C1 = real(eval(C1));

% visual comparison between the theoretical solution for Light On condition
% and the empirical profile
tt = 0:0.01:1000;
O1 = C1*exp(lamda1*tt) +C2*exp(lamda2*tt)+C3*exp(lamda3*tt) + yp(1);
O2 = A21*C1*exp(lamda1*tt) + A22*C2*exp(lamda2*tt) + A23*C3*exp(lamda3*tt) + yp(2);
V = -100;

I = V*g1.*(O1+gama.*O2);
Iplot = I./max(-I);
figure; plot([zeros(1,500/0.01) Iplot],'r'); hold on;
Iemp = monoexpETC(); % change this function to match the variant
plot(-Iemp);

% the maximum amplitude of the photocurrent component associated with each
% eigenvalue and steady state
Ic1 = V*g1*C1*(1+gama*A21);
Ic2 = V*g1*C2*(1+gama*A22);
Ic3 = V*g1*C3*(1+gama*A23);
Ip = V*g1*(yp(1) + gama*yp(2));

l = [abs(lamda1) abs(lamda2) abs(lamda3)]; tau_on = 1./l
Ic = [Ic1 Ic2 Ic3 Ip]; Aon = Ic./(V*g1)

% the photocurrent component associated with each eigenvalue
Il1 = V*g1*C1*(1+gama*A21)*exp(lamda1*tt);
Il2 = V*g1*C2*(1+gama*A22)*exp(lamda2*tt);
Il3 = V*g1*C3*(1+gama*A23)*exp(lamda3*tt);

%%%%%%%%%%%%%%%%%%%%%% LIGHT OFF CONDITION   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% evaluation of the eigenvalues
b = (Gd1+Gd2+e12+e21)/2;
c = sqrt(b^2 - (Gd1*Gd2+Gd1*e21+Gd2*e12));

L1 = b-c; L2 = b+c;
L = [L1 L2]; tau_off = 1./L
O10 = O1(end); O20 = O2(end);

% check O10 and O20 with the steady state from light on condition
syms C2s O1s O2s

C2s = solve('P1*( 1-C2s -((Gd2+e21)*(P2+Gr)/Gd2 - P2)*C2s/e12 - (P2+Gr)*C2s/Gd2) - (Gd1+e12)*((Gd2+e21)*(P2+Gr)/Gd2 - P2)*C2s/e12 +e21*(P2+Gr)*C2s/Gd2 = 0',C2s);

O10bis = eval( ((Gd2+e21)*(P2+Gr)/Gd2 - P2 )*C2s/e12 );
O20bis = eval( (P2+Gr)*C2s/Gd2 );

%evaluation of the slow and fast components
Aslow = ( (L2-(Gd1+(1-gama)*e12))*O10 +((1-gama)*e21+gama*(L2-Gd2))*O20 )/(L2-L1);

Afast =  ( (Gd1+(1-gama)*e12-L1)*O10 +(gama*(Gd2-L1)-(1-gama)*e21)*O20 )/(L2-L1);

Islow = V*g1*Aslow*exp(-L1*tt); Ifast = V*g1*Afast*exp(-L2*tt);

Aoff = [Aslow Afast]










