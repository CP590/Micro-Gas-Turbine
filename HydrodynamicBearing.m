r = 3e-3;
d = 2*r;
l = 5e-3;
cd = 0.01*r;
cr = 0.5*cd;
mu = 30e-3;
N = 200e3;
n = N/60;
U = pi*n*d;

theta = 30;

eps = 0;
e = eps*cr;
phi = 30;
On = 4*pi*eps;
Keps = (eps*(pi^2*(1-eps^2) + 16*eps^2)^0.5)/(4*(1-eps^2)^2);
F = Keps*4*pi*mu*d*n*l^3/cd^2;

z = 0:l/100:0.5*l;
for i = 1:length(z)
P(i) = mu*U*(0.25*l^2 - z(i)^2)*(3*eps*sind(theta))/(r*cr^2*(1+eps*cosd(theta))^3);
end

Ts = 4*pi^2*mu*r^3*l*n/(cr*(1-eps^2)^0.5);
Tr = Ts + F*e*sind(phi);
Pl = 2*pi*Tr*n;
