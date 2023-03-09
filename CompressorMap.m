function [rcOD, etacOD, A2, m, mc, r1s, b2] = CompressorMap(P01, T01, k, Cp, R, mu, ka, etac, r1h, r2, alpha2m, beta1h, Z, N, betab2d, md)
m = 0e-4:1e-4:30e-3;

for j = 1:length(N)
[rc, r1s, b2, A1, Ns, beta1s, U2, M1s, M2, M2rel, sigma] = CDesignPoint(P01, T01, k, Cp, R, etac, r1h, r2, alpha2m, beta1h, Z, md, N(j), -betab2d);    
for i = 1:length(m)
[rcOD(j, i), etacOD(j, i) A2, T02OD(j, i)] = COffDesign(P01, T01, k, Cp, R, mu, ka, etac, r1h, r1s, r2, b2, A1, alpha2m, beta1h, Z, m(i), N(j), betab2d, sigma);
mc(j, i) = m(i)*sqrt(T02OD(j, i)/T01)/rcOD(j, i);
end

%{
plot(mc, rcOD)
grid on
grid minor
hold on
%}
end