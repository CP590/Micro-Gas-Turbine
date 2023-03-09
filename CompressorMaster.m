clear all
P01 = 101325;
T01 = 288;
[k, R, Cp] = AirProperties(T01, P01);
%{
k = 1.4;
Cp = 1005;
R = 287;
%}
mu = 18.27e-6;
ka = 1e-5;

etac = 0.73;
r1h = 3e-3;
r2 = 20e-3;
alpha2m = 65;
beta1h = 30;
Z = 12;
md = 20e-3;
betab2d = 15;

N = 110e3:10e3:190e3;

%[rc, r1s, b2, A1, Ns, beta1s, U2, M1s, M2, M2rel, sigma] = CBacksweepEval(P01, T01, k, Cp, R, etac, r1h, r2, alpha2m, beta1h, Z, md, N);
[rcOD, etacOD, A2, m, mc, r1s, b2] = CompressorMap(P01, T01, k, Cp, R, mu, ka, etac, r1h, r2, alpha2m, beta1h, Z, N, betab2d, md);
[Xh, Yh, Xs, Ys, t] = MeridionalDefinition(r1s, r1h, r2, b2);
[Bh, Bs, Phih, Phis] = CamberDefinition(t, Yh, Ys);
[dh, ds] = ThicknessDefinition(t);
[Yhn, Yhs, Yhp, Ysn, Yss, Ysp, Zhn, Zhs, Zhp, Zsn, Zss, Zsp] = ThreeDDefinition(Xh, Xs, Yh, Ys, Z, Phih, Phis, dh, ds);



%{
for j = 1:length(N)
for i = 1:length(betab2)
[rc(i), r1s, b2(i), A1, Ns, beta1s, U2, M1s, M2, M2rel, sigma] = DesignPoint(P01, T01, k, Cp, R, etac, r1h, r2, alpha2m, beta1h, Z, md, N(j), betab2(i));
end
%plot(betab2, rc)
%grid on
%grid minor
%hold on
%end

%
for i = 1:length(m)
[rcOD(i), A2] = OffDesign(P01, T01, k, Cp, R, etac, r1h, r1s, r2, b2, A1, alpha2m, beta1h, Z, m(i), N(j), betab2d, sigma);
end

plot(m, rcOD)
grid on
grid minor
hold on
end
%

%[mcc, rcc] = ChokeLine(P01, T01, k, R, etac, A2);
%plot(mcc, rcc)
%}