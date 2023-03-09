function[k, R, Cp] = AirProperties(T, P)
k = 1.42592 - 8.03974e-5*T;
h = 0.919848*T^1.01457;
R = 286.99;
Cp = R*k/(k-1);
Cv = Cp - R;
rho = P/(R*T);