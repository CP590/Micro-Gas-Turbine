clear all
Pa = 101325;
Ta = 288;
k = 1.333;
Cp = 1148;
R = 287;
P02 = rc*Pa;
delPb = P02*0.04;
P01 = Pa*rc*(1-delPb);
T03 = 1200;
P03 = Pa;

etat = 0.75;
r1h = 3e-3;
r2 = 20e-3;
beta1h = 30;
Z = 12;
md = 20e-3;

N = 140e3;
%20e3:220e3;

%alpha2 = sqrt(acosd(0.63*pi/(2*Z)));
beta2 = -acosd(1-(0.63*pi)/Z);
Sw = etat(1-(P01/P03)^((1-k)/k));
T01 = T03*(1-Sw);
rho01 = P01/(R*T01);
rho03 = P03/(R*T03);
a01 = sqrt(k*R*T01);
U2 = a01*sqrt((1/(k-1))-Sw/(cosd(beta2)));
r2 = U2*60/(2*pi*N);
C2 = sqrt(a01^2*(Sw/(k-1))*(2*cosd(beta2))/(1+cosd(beta2)));
alpha2 = atand(sind(beta2)/(cosd(beta2)-1));
T02 = T03;
T2 = T02 - C2^2/(2*Cp);
M2 = C2/sqrt(k*R*T2);
Cm2 = C2*cos(alpha2);
Ctheta2 = C2*sin(alpha2);
W2 = sqrt(Cm2^2 + (U2-Ctheta2)^2);
%Wout = m*(U2*Ctheta2 
%Wnet = 

M3srel = 0.5;
M3 = M3srel*cos(beta3s);


