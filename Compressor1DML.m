P01 = 101325;
T01 = 288;
k = 1.4;
Cp = 1005;
R = 287;

etac = 0.73;
r1h = 3e-3;
r2 = 20e-3;
alpha2m = 65;
beta1h = 30;
Z = 12;
m = 20e-3;

N = 160e3;
betab2 = 10;

sigma = 1 - (sqrt(cosd(betab2)))/Z;
%sigma= 1 - (cosd(betab2))/Z^0.7;
lamda = sigma/(1-(tand(betab2)/tand(alpha2m)));
U1h = 2*pi*N*r1h/60;
Cm1 = U1h/tand(beta1h);
T1 = T01 - (Cm1.^2/(2*Cp));
a1 = sqrt(k*R*T1);
rho01 = P01/(R*T01);
U2 = 2*pi*N*r2/60;
Mu = U2/a1;
rc = (1+(k-1)*etac*lamda*Mu^2)^(k/(k-1));
P02 = P01*rc;
Ctheta2 = lamda*U2;
Cm2 = Ctheta2/tand(alpha2m);
%Cm2 = sigma*U2/(tand(alpha2m)+sigma*tand(betab2));
%Ctheta2 = sigma*(U2-Cm2*tand(betab2));
T02 = T01*(1+(k-1)*lamda*Mu^2);
rho02 = P02/(R*T02);
delh = (U2*Ctheta2);
W1h = U1h/sind(beta1h);
M1h = U1h/a1;
P1 = P01/((T01./T1)^(k/(k-1)));
rho1 = P1/(R*T1);
A1 = m/(rho1*Cm1);
r1s = sqrt((A1/pi) + r1h^2);
U1s = 2*pi*N*r1s/60;
W1s = sqrt(Cm1^2 + U1s^2);
M1s = U1s/a1;
beta1s = atand(U1s/Cm1);
C2 = sqrt(Ctheta2^2 + Cm2^2);
T2 = T02 - C2^2/(2*Cp);
a2 = sqrt(k*R*T2);
W2 = sqrt(Cm2^2 + (U2 - Ctheta2)^2);
M2rel = W2/a2;
M2 = C2/a2;
P2 = P02/((T02./T2)^(k/(k-1)));
rho2 = P2/(R*T2);
A2 = m/(rho2*Cm2);
b2 = A2/(2*pi*r2)
Ns = 2*pi*N*sqrt(m/(0.5*(rho01+rho02)))/(60*delh^0.75)
DR2 = W1s/W2;
MR2 = M1s/M2rel;


A1 = pi*(r1s^2 - r1h^2);
%{
for j = 1:1
    if j == 1
        ms(j) = m;
    end
for i = 1:10
    if i == 1 && j == 1
    T1s(i) = T01;
    elseif i == 1 && j > 1
    T1s(i) = T1(j-1);
    else
    T1s(i) = T1s(i-1);
    end
    P1s(i) = P01/(T01/T1s(i))^(k/(k-1));
    rho1s(i) = P1s(i)/(R*T1s(i));
    Cm1s(i) = ms(j)/(A1*rho1s(i));
    T1s(i) = T01 - Cm1s(i)^2/(2*Cp);
    if i == 10    
    T1(j) = T1s(i);
    P1(j) = P1s(i);
    rho1(j) = rho1s(i);
    Cm1(j) = Cm1s(i);
    end
    end

    
a1(j) = sqrt(k*R*T1(j));
a01 = sqrt(k*R*T01);
U1h = 2*pi*N*r1h/60;
alpha1h(j) = atand(U1h/Cm1(j));
W1h(j) = sqrt(U1h^2 + Cm1(j)^2);
M1h(j) = W1h(j)/a1(j);
U1s = 2*pi*N*r1s/60;
alpha1s(j) = atand(U1s/Cm1(j));
W1s(j) = sqrt(U1s^2 + Cm1(j)^2);
M1s(j) = W1s(j)/a1(j);
U1 = ((r1s - r1h)/2 + r1h)*2*pi*N/60;
alpha1(j) = atand(U1/Cm1(j));
W1(j) = real(sqrt(U1^2 + Cm1(j)^2));
M1(j) = W1(j)/a1(j);

A2 = 2*pi*r2*b2;
U2 = 2*pi*N*r2/60;

    Cm2(j) = 2*Cm1(j);
    Cptheta2(j) = U2 - Cm2(j)*tand(betab2);
    Cp2(j) = sqrt(Cptheta2(j)^2 + Cm2(j)^2);
    W2(j) = sqrt(Cm2(j)^2 + (U2 - Cptheta2(j))^2);
    Cslip = U2*(1 - sigma);
    Ctheta2(j) = Cptheta2(j) - Cslip;
    C2(j) = sqrt(Ctheta2(j)^2 + Cm2(j)^2);
    W2(j) = sqrt(Cm2(j)^2 + (U2 - Ctheta2(j))^2);
    alpha2m(j) = atand(Ctheta2(j)/Cm2(j));
    beta2(j) = atand((U2 - Ctheta2(j))/Cm2(j));
    delh(j) = U2*Ctheta2(j);
    rc(j) = (1+((k-1)/a01^2)*etac*delh(j))^(k/(k-1));
    P02(j) = P01*rc(j);
    T02(j) = T01*((P02(j)/P01)^((k-1)/k)-1)/etac + T01;
    M2(j) = (U2*Ctheta2(j))/((a01*U2*sind(alpha2m(j)))*sqrt(1+(k-1)*(U2/a01)^2*(Ctheta2(j)/U2)*(1-(Ctheta2(j)/(2*U2*(sind(alpha2m(j)))^2)))));
    P2(j) = P02(j)/((1 + 0.5*(k-1)*M2(j)^2)^(k/(k-1)));
    T2(j) = T02(j)/((P02(j)/P2(j))^(k/(k-1)));
    rho2(j) = P2(j)/(R*T2(j));
    a2(j) = sqrt(k*R*T2(j));
    M2rel(j) = W2(j)/a2(j);
    ms(j+1) = rho2(j)*A2*Cm2(j);
end
%}



for i = 1:10
    if i == 1
    T1s(i) = T01;
    else
    T1s(i) = T1s(i-1);
    end
    P1s(i) = P01/(T01/T1s(i))^(k/(k-1));
    rho1s(i) = P1s(i)/(R*T1s(i));
    Cm1s(i) = m/(A1*rho1s(i));
    T1s(i) = T01 - Cm1s(i)^2/(2*Cp);
    end

A2 = 2*pi*r2*b2;
U2 = 2*pi*N*r2/60;
for j = 1:10
    if i == 10    
    T1(j) = T1s(i);
    P1(j) = P1s(i);
    rho1(j) = rho1s(i);
    Cm1(j) = Cm1s(i);
    a1(j) = sqrt(k*R*T1(j));
a01 = sqrt(k*R*T01);
U1h = 2*pi*N*r1h/60;
alpha1h(j) = atand(U1h/Cm1(j));
W1h(j) = sqrt(U1h^2 + Cm1(j)^2);
M1h(j) = W1h(j)/a1(j);
U1s = 2*pi*N*r1s/60;
alpha1s(j) = atand(U1s/Cm1(j));
W1s(j) = sqrt(U1s^2 + Cm1(j)^2);
M1s(j) = W1s(j)/a1(j);
U1 = ((r1s - r1h)/2 + r1h)*2*pi*N/60;
alpha1(j) = atand(U1/Cm1(j));
W1(j) = real(sqrt(U1^2 + Cm1(j)^2));
M1(j) = W1(j)/a1(j);
    end
    if j == 1
    Cm2(j) = 2*Cm1(j);
    else
    Cm2(j) = Cm2(j-1);
    end
    Cptheta2(j) = U2 - Cm2(j)*tand(betab2);
    Cp2(j) = sqrt(Cptheta2(j)^2 + Cm2(j)^2);
    W2(j) = sqrt(Cm2(j)^2 + (U2 - Cptheta2(j))^2);
    Cslip = U2*(1 - sigma);
    Ctheta2(j) = Cptheta2(j) - Cslip;
    C2(j) = sqrt(Ctheta2(j)^2 + Cm2(j)^2);
    W2(j) = sqrt(Cm2(j)^2 + (U2 - Ctheta2(j))^2);
    alpha2m(j) = atand(Ctheta2(j)/Cm2(j));
    beta2(j) = atand((U2 - Ctheta2(j))/Cm2(j));
    delh(j) = U2*Ctheta2(j);
    rc(j) = (1+((k-1)/a01^2)*etac*delh(j))^(k/(k-1));
    P02(j) = P01*rc(j);
    T02(j) = T01*((P02(j)/P01)^((k-1)/k)-1)/etac + T01;
    M2(j) = (U2*Ctheta2(j))/((a01*U2*sind(alpha2m(j)))*sqrt(1+(k-1)*(U2/a01)^2*(Ctheta2(j)/U2)*(1-(Ctheta2(j)/(2*U2*(sind(alpha2m(j)))^2)))));
    P2(j) = P02(j)/((1 + 0.5*(k-1)*M2(j)^2)^(k/(k-1)));
    T2(j) = T02(j)/((P02(j)/P2(j))^(k/(k-1)));
    rho2(j) = P2(j)/(R*T2(j));
    a2(j) = sqrt(k*R*T2(j));
    M2rel(j) = W2(j)/a2(j);
    %ms(j) = rho2(j)*A2*Cm2(j);
    Cm2(j) = m/(rho2(j)*A2);
end
