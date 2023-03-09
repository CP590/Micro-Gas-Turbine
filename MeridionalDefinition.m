function[Xh, Yh, Xs, Ys, t] = MeridionalDefinition(r1s, r1h, r2, b2)

x0h = 0;
x1h = 2*r1s;
x2h = 0.9*r2;
x3h = x2h;

x0s = 0;
x1s = 0.9*r2;
x2s = x3h - 1.2*b2;
x3s = x3h - b2;

y0h = r1h;
y1h = 1.1*r1h;
y2h = 0.3*r2;
y3h = r2;

y0s = r1s;
y1s = r1s;
y2s = 0.8*r2;
y3s = r2;

Ah = x3h - 3*x2h + 3*x1h - x0h;
Bh = 3*x2h - 6*x1h + 3*x0h;
Ch = 3*x1h - 3*x0h;
Dh = x0h;

As = x3s - 3*x2s + 3*x1s - x0s;
Bs = 3*x2s - 6*x1s + 3*x0s;
Cs = 3*x1s - 3*x0s;
Ds = x0s;

Eh = y3h - 3*y2h + 3*y1h - y0h;
Fh = 3*y2h - 6*y1h + 3*y0h;
Gh = 3*y1h - 3*y0h;
Hh = y0h;


Es = y3s - 3*y2s + 3*y1s - y0s;
Fs = 3*y2s - 6*y1s + 3*y0s;
Gs = 3*y1s - 3*y0s;
Hs = y0s;

t = 0:0.01:1;
Xs = (((As.*t) + Bs).*t + Cs).*t + Ds;
Ys = (((Es.*t) + Fs).*t + Gs).*t + Hs;
dXsdt = 3*As.*t.^2 + 2*Bs.*t + Cs;
dYsdt = 3*Es.*t.^2 + 2*Fs.*t + Gs;
thetas = atand(dYsdt./dXsdt);

Xh = (((Ah.*t) + Bh).*t + Ch).*t + Dh;
Yh = (((Eh.*t) + Fh).*t + Gh).*t + Hh;
dXhdt = 3*Ah.*t.^2 + 2*Bh.*t + Ch;
dYhdt = 3*Eh.*t.^2 + 2*Fh.*t + Gh;
thetah = atand(dYhdt./dXhdt);
thetah2 = atand((Xs - Xh)./(Ys - Yh));

%{
plot(Xh, Yh)
hold on
plot(Xs, Ys)
plot(x0h, y0h, '.')
plot(x1h, y1h, '.')
plot(x2h, y2h, '.')
plot(x3h, y3h, '.')
plot(x0s, y0s, '*')
plot(x1s, y1s, '*')
plot(x2s, y2s, '*')
plot(x3s, y3s, '*')
%}
end