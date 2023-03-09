function[Bh, Bs, Phih, Phis] = CamberDefinition(t, Yh, Ys)

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

for i = 1:length(t)
funh(i) = tand(Bh(i))/Yh(i);
funs(i) = tand(Bs(i))/Ys(i);
end

for i = 2:length(t)
    u = 0:0.01:(i-1)/(length(t)-1);
    Phih(i) = trapz(u, funh(1:i));
    Phis(i) = trapz(u, funs(1:i));
end

%plot(t, Phih)
%plot(t, Phis)

end

