function [rcOD, etacOD, A2, T02OD] = COffDesign(P01, T01, k, Cp, R, mu, ka, etac, r1h, r1s, r2, b2, A1, alpha2m, beta1h, Z, m, N, betab2d, sigma)

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
    Cptheta2(j) = U2 - Cm2(j)*tand(betab2d);
    Cp2(j) = sqrt(Cptheta2(j)^2 + Cm2(j)^2);
    W2(j) = sqrt(Cm2(j)^2 + (U2 - Cptheta2(j))^2);
    Cslip = U2*(1 - sigma);
    Ctheta2(j) = Cptheta2(j) - Cslip;
    C2(j) = sqrt(Ctheta2(j)^2 + Cm2(j)^2);
    W2(j) = sqrt(Cm2(j)^2 + (U2 - Ctheta2(j))^2);
    alpha2m(j) = atand(Ctheta2(j)/Cm2(j));
    beta2(j) = atand((U2 - Ctheta2(j))/Cm2(j));
    delh(j) = U2*Ctheta2(j);
    rc(j) = real((1+((k-1)/a01^2)*etac*delh(j))^(k/(k-1)));
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
    if j == 10
        rcOD = rc(j);
        T02OD = T02(j);
        T02sOD = T01*rcOD^((k-1)/k);
        etacOD = (T02sOD - T01)/(T02OD-T01);
        ua = 0.5*(U1+U2);
        nua = 0.5*mu*((1/rho1(j)) + (1/rho2(j)));
        ea = 0.5*ka/b2;
        Rea = 2*b2*ua/nua;
    for k = 1:11
    if k == 1
    f(k) = 0.0001;
    else
    f(k) = f(k-1);
    end
    f(k) = (-2*log((ea/3.7) + (2.51/(Rea*sqrt(f(k))))))^-2;
    if k == 10
        fa = f(k);
    end
    end
    c = 0.25*r2/b2;
    %etacOD = 1 - (c*fa*U2/Ctheta2(j));
end
end