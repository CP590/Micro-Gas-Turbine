rho = 1.204;
v = 40;
l= 10e-3;
Cp = 992.532;
k = 0.0262;
alpha = 0.003661;
mu = 0.019;
theta = 0;
g = 9.81;
Tw = 293.15;
Tfl = 1000;

Re = rho*v*l/mu;
Pr = mu*Cp/k;
Gr = l^3*g*rho^2*alpha*abs(Tw-Tfl)/mu^2;
Ra = Pr*Gr;