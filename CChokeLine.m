function [mcc, rcc] = CChokeLine(P01, T01, k, R, etac, A2)

rcc = 1:0.01:2.6;
for j = 1:length(rcc)
P02(j) = P01*rcc(j);
T02(j) = T01*((P02(j)/P01)^((k-1)/k)-1)/etac + T01;
mc(j) = A2*P02(j)*sqrt(k/R)*(0.5*(k-1))^(-(k+1)/(2*(k-1)))/sqrt(T02(j));
mcc(j) = mc(j)*sqrt(T02)*P01/(sqrt(T01)*P02);
end