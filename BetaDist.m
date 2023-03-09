function[Bh, Bs, thetah, thetas, funh, funs] = BetaDist(t, Yh, Ys)
B0h = 30;
B1h = 0;
B2h = 0;
B3h = -15;

B0s = 61;
B1s = 0;
B2s = 0;
B3s = -15;

Bh = B0h.*(1-t).^3 + 3*B1h.*t.*(1-t).^2 + 3*B2h.*t.^2.*(1-t) + B3h.*t.^3;
Bs = B0s.*(1-t).^3 + 3*B1s.*t.*(1-t).^2 + 3*B2s.*t.^2.*(1-t) + B3s*t.^3;

%{
for i = 1:length(t)
    if i == 1
        Bh(i) = B0h;
        Bs(i) = B0s;
    else
u = 0:0.05:(i-1)*t/(length(t)-1)
Bh(i) = B0h.*(1-u).^3 + 3*B1h.*u.*(1-u).^2 + 3*B2h.*u.^2.*(1-u) + B3h.*u.^3;
Bs(i) = B0s.*(1-u).^3 + 3*B1s.*u.*(1-u).^2 + 3*B2s.*u.^2.*(1-u) + B3s.*u.^3;
thetah(i) = trapz(u, tand(Bh(1:i))/Yh(1:i));
thetas(i) = trapz(u, tand(Bs(1:i))/Ys(1:i));
    end

%}

for i = 1:length(t)
funh(i) = tand(Bh(i))/Yh(i);
funs(i) = tand(Bs(i))/Ys(i);
end



for i = 2:length(t)
    u = 0:0.05:(i-1)/(length(t)-1)
    thetah(i) = trapz(u, funh(1:i))
    thetas(i) = trapz(u, funs(1:i))
end

plot(t, thetah)
hold on
plot(t, thetas)

%{
u = 0:0.05:1;
Bh = B0h.*(1-u).^3 + 3*B1h.*u.*(1-u).^2 + 3*B2h.*u.^2.*(1-u) + B3h.*u.^3;
Bs = B0s.*(1-u).^3 + 3*B1s.*u.*(1-u).^2 + 3*B2s.*u.^2.*(1-u) + B3s.*u.^3;
for i = 5:length(t)
    i
    tand(Bh(2:i))/Yh(2:i)
thetah(i) = trapz(u, tand(Bh(2:i))/Yh(2:i));
thetas(i) = trapz(u, tand(Bs(2:i))/Ys(2:i));
end

%}

end

