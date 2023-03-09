%{
        ua = 0.5*(U1+ U2);
        nua = 0.5*mu*(1/rho1(j) + 1/rho2(j));
        ea = 0.5*ka/b2;
        Rea = 2*b2*ua/nua;
    for k = 1:10
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
    etacc = 1 - (c*fa*U2/Ctheta2);
    %}