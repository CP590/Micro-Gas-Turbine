function [rc, r1s, b2, A1, Ns, beta1s, U2, M1s, M2, M2rel, sigma] = CBacksweepEval(P01, T01, k, Cp, R, etac, r1h, r2, alpha2m, beta1h, Z, md, N)
betab2 = 0:1:70;
for j = 1:length(N)
for i = 1:length(betab2)
[rc(i), r1s, b2(i), A1, Ns, beta1s, U2, M1s, M2, M2rel, sigma] = CDesignPoint(P01, T01, k, Cp, R, etac, r1h, r2, alpha2m, beta1h, Z, md, N(j), -betab2(i));
end
plot(betab2, b2)
grid on
grid minor
hold on
end
end