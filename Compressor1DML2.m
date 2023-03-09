clear all
P01 = 101325;
T01 = 288;
k = 1.4;
Cp = 1005;
R = 287;
mu = 18.27e-6;
ka = 1e-5;

rcd = 2.15;
r1h = 3e-3;
r2 = 20e-3;
alpha2m = 65;
beta1h = 30;
Z = 12;
md = 20e-3;
betab2d = 15;

N = 150e3;

sigma = 1-sqrt(cosd(betab2d))/Z;
lamda = sigma/(1-(tand(betab2d)/tand(alpha2m)));
U1h = 2*pi*N*r1h/60;
Cm1 = U1h/tand(beta1h);
T1 = T01-Cm1^2/(2*Cp);
P1 = P01/(T01/T1)^(k/(k-1));
rho1 = P1/(R*T1);
A1 = md/(rho1*Cm1);
r1s = sqrt(A1/pi+r1h^2);
U1s = 2*pi*N*r1s/60;
U2 = 2*pi*N*r2/60;
Ctheta2m = lamda*U2;
Cm2 = Ctheta2m/tand(alpha2m);
P02 = rcd*P01;
a01 = sqrt(k*R*T01);
Mu = U2/a01;
T02 = T01*(1+(k-1)*lamda*Mu^2);
T2 = T02-Cm2^2/(2*Cp);
P2 = P02/(T02/T2)^(k/(k-1));
rho2 = P2/(R*T2);
A2 = md/(rho2*Cm2);
b2 = A2/(2*pi*r2);
etac = (rcd^((k-1)/k)-1)/(Mu^2*lamda*(k-1));

W1h = U1h/sind(beta1h);
W1s = sqrt(Cm1^2+U1s^2);
W2 = sqrt(Cm2^2 + (U2-Ctheta2m)^2);
a02 = sqrt(k*R*T02);
M1h = U1h/a01;
M1s = U1s/a01;
M2rel = W2/a02;

lc = 20e-3;
B0h = [0; r1h];
B1h = [0.2*lc; 0.2*r2];
B2h = [lc; 0.5*r2];
B3h = [lc; r2];
B0s = [0; r1s];
B1s = [0.2*lc; r1s];
B2s = [0.8*lc; 0.5*r2];
B3s = [(lc-b2); r2];

t = 0:0.01:1;
Ph = kron((1-t).^3, B0h) + kron(3*(1-t).^2.*t, B1h) + kron(3*(1-t).*t.^2, B2h) + kron(t.^3, B3h);
Ps = kron((1-t).^3, B0s) + kron(3*(1-t).^2.*t, B1s) + kron(3*(1-t).*t.^2, B2s) + kron(t.^3, B3s);

thh = 500e-3;
ths = 500e-3;
Phihsp = atand(0.5.*thh.*Ph(2, :));
Phissp = atand(0.5.*ths.*Ps(2, :));
%{
Phi0h
Phi1h
Phi2h
Phi3h
%}
Phih = 10*(t.^1);
Phis = 10*(t.^1);
n = 12;
Yhp = Ph(2, :).*cosd(Phih + 360*n/Z - Phihsp);
Yh = Ph(2, :).*cosd(Phih + 360*n/Z);
Yhs = Ph(2, :).*cosd(Phih + 360*n/Z + Phihsp);
Ysp = Ps(2, :).*cosd(Phis + 360*n/Z - Phissp);
Ys = Ps(2, :).*cosd(Phis + 360*n/Z);
Yss = Ps(2, :).*cosd(Phis + 360*n/Z + Phissp);

Zhp = Ph(2, :).*sind(Phih + 360*n/Z - Phihsp);
Zh = Ph(2, :).*sind(Phih + 360*n/Z);
Zhs = Ph(2, :).*sind(Phih + 360*n/Z + Phihsp);
Zsp = Ps(2, :).*sind(Phis + 360*n/Z - Phissp);
Zs = Ps(2, :).*sind(Phis + 360*n/Z);
Zss = Ps(2, :).*sind(Phis + 360*n/Z + Phissp);

plot3(Ph(1, :), Zh, Yh)
grid minor
hold on
plot3(Ph(1, :), Zhp, Yhp)
plot3(Ph(1, :), Zhs, Yhs)
plot3(Ps(1, :), Zs, Ys)
plot3(Ps(1, :), Zsp, Ysp)
plot3(Ps(1, :), Zss, Yss)
pbaspect([1 1 1])