Wnet = 1e3;
delh = 43100e3;
etat = 0.7;
etac = 0.7;
etam = 0.9;
etab = 0.99;
etahex = 0.75;
k1 = 1.4;
k2 = 1.333;
Cpa = 1005;
Cpg = 1148;
T3 = 800:10:1200;
T1 = 288.15;
rc = 2.15;
P1 = 1;
P4 = P1;
P2 = rc*P1;
delPb = P2*0.04;
delPhg = P1*0.03;
delPha = P2*0.03;
P3 = P2;
rt = P3/P4;
T2s = T1*rc^((k1-1)/k1);
T4s = T3*(1/rt)^((k2-1)/k2);
T2 = (T2s -T1)/etac + T1;
T4 = T3*(1-(etat*(1-(1/rt)^((k2-1)/k2))));
T5 = etahex*(T4-T2) + T2;
wnet = Cpg*(T3 - T4) - ((1/etam)*Cpa*(T2 - T1));
etas = (etab*((Cpg*(T3 - T4)) - ((1/etam)*Cpa*(T2 - T1))).*(etab*delh - Cpg.*T3))./(delh.*(Cpg.*T3 - Cpa*T2));
etas2 = (etab*((Cpg*(T3 - T4)) - ((1/etam)*Cpa*(T2 - T1))).*(etab*delh - Cpg.*T3))./(delh.*(Cpg.*T3 - Cpa*T5));
f = (Cpg*T3 - Cpa*T2)/(etab*delh - Cpg*T3);
mdot = Wnet./wnet;
Wnetv = min(mdot).*wnet;

plot(Wnetv, etas)
hold on
grid minor
plot(Wnetv, etas2)

ma = 5e-3:0.1e-4:40e-3;
mf = 0.2e-3;
%T03 = (ma.*Cpa*T2 + mf*delh)./((ma+mf).*Cpg);
T032 = (ma.*Cpa*(T2*(1-etahex))+(mf*delh))./(((ma+mf).*Cpg).*(1-(ma.*etahex*(1-(etat*(1-(1/rt)^((k2-1)/k2)))))/(ma+mf)));
T04 = T032*(1-etat*(1-(1/rc)^((k2-1)/k2)));
T05 = T2*(1-etahex) + etahex*T04;
f2 = (Cpg*T032 - Cpa*T2)/(delh - Cpg*T032);
Wdnet = ((ma+mf).*(T032-(T032*(1-(etat*(1-(1/rt)^((k2-1)/k2))))))) - (ma.*(T2-T1));
error =  ma.*Cpa.*T05 + mf*delh - (ma+mf).*Cpg.*T032;
