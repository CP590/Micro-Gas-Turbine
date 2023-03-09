function[dh, ds] = ThicknessDefinition(t)

d0h = 0;
d1h = 10e-3;
d2h = 1e-4;
d3h = 0;

d0s = 0;
d1s = 10e-3;
d2s = 1e-4;
d3s = 0;

dh = d0h.*(1-t).^3 + 3*d1h.*t.*(1-t).^2 + 3*d2h.*t.^2.*(1-t) + d3h.*t.^3;
ds = d0s.*(1-t).^3 + 3*d1s.*t.*(1-t).^2 + 3*d2s.*t.^2.*(1-t) + d3s*t.^3;
%{
plot(t, dh)
hold on
plot(t, ds)
%}
end