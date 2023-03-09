function[Yhn, Yhs, Yhp, Ysn, Yss, Ysp, Zhn, Zhs, Zhp, Zsn, Zss, Zsp] = ThreeDDefinition(Xh, Xs, Yh, Ys, Z, Phih, Phis, dh, ds)

n = Z;

Phihsp = atand(0.5.*dh.*Yh);
Phissp = atand(0.5.*ds.*Ys);

Yhs = Yh.*cosd(Phih + Phihsp + (360*n/Z));
Yhn = Yh.*cosd(Phih + (360*n/Z));
Yhp = Yh.*cosd(Phih - Phihsp + (360*n/Z));
Yss = Ys.*cosd(Phis + Phissp + (360*n/Z));
Ysn = Ys.*cosd(Phis + (360*n/Z));
Ysp = Ys.*cosd(Phis - Phissp + (360*n/Z));
Zhs = Yh.*sind(Phih + Phihsp + (360*n/Z));
Zhn = Yh.*sind(Phih + (360*n/Z));
Zhp = Yh.*sind(Phih - Phihsp + (360*n/Z));
Zss = Ys.*sind(Phis + Phissp + (360*n/Z));
Zsn = Ys.*sind(Phis + (360*n/Z));
Zsp = Ys.*sind(Phis - Phissp + (360*n/Z));

theta = 0:pi/100:2*pi;
x1 = min(Yh)*cos(theta);
y1 = min(Yh)*sin(theta);
z1 = x1*0;
x2 = max(Yh)*cos(theta);
y2 = max(Yh)*sin(theta);
z2 = ones(1, length(x2))*max(Xh);

plot3(Xh, Yhn, Zhn)
%%plot(Xh, Zhn)
grid on
hold on
plot3(z1, x1, y1)
plot3(z2, x2, y2)
%%plot(Xh, Zhp)
plot3(Xh, Yhs, Zhs)
plot3(Xh, Yhp, Zhp)
plot3(Xs, Ysn, Zsn)
%%plot(Xs, Zsn)
plot3(Xs, Yss, Zss)
plot3(Xs, Ysp, Zsp)

axis square

end