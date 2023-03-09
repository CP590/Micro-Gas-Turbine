clear all
R = 8.314;
Po = 101325;
Pc = 3*Po;


T1CH4 = 200.000:1:1000;
T2CH4 = 1000:1:6000.000;
MCH4 = 16.04276;
cCH4 = [1.63552643E+00 1.00842795E-02 -3.36916254E-06 5.34958667E-10 -3.15518833E-14 -1.00056455E+04 9.99313326E+00 5.14987613E+00 -1.36709788E-02 4.91800599E-05 -4.84743026E-08 1.66693956E-11 -1.02466476E+04 -4.64130376E+00 -8.97226656E+03];
H1CH4 = (cCH4(1) + cCH4(2).*T1CH4/2 + cCH4(3).*T1CH4.^2/3 + cCH4(4).*T1CH4.^3/4 + cCH4(5).*T1CH4.^4/5 + cCH4(6)./T1CH4)*R.*T1CH4;     
H2CH4 = (cCH4(8) + cCH4(9).*T2CH4/2 + cCH4(10).*T2CH4.^2/3 + cCH4(11).*T2CH4.^3/4 + cCH4(12).*T2CH4.^4/5 + cCH4(13)./T2CH4)*R.*T2CH4; 

T1O2 = 200:1:1000;
T2O2 = 1000:1:6000;
MO2 = 31.99880;
cO2 = [3.66096083E+00 6.56365523E-04 -1.41149485E-07 2.05797658E-11 -1.29913248E-15 -1.21597725E+03 3.41536184E+00 3.78245636E+00 -2.99673415E-03 9.84730200E-06 -9.68129508E-09 3.24372836E-12 -1.06394356E+03 3.65767573E+00 0.00000000E+00];
H1O2 = (cO2(1) + cO2(2).*T1O2/2 + cO2(3).*T1O2.^2/3 + cO2(4).*T1O2.^3/4 + cO2(5).*T1O2.^4/5 + cO2(6)./T1O2)*R.*T1O2;
H2O2 = (cO2(8) + cO2(9).*T2O2/2 + cO2(10).*T2O2.^2/3 + cO2(11).*T2O2.^3/4 + cO2(12).*T2O2.^4/5 + cO2(13)./T2O2)*R.*T2O2;


T1H2O = 200:1:1000;
T2H2O = 1000:1:6000;
MH2O = 18.01528;
cH2O = [2.67703787E+00 2.97318329E-03 -7.73769690E-07 9.44336689E-11 -4.26900959E-15 -2.98858938E+04 6.88255571E+00 4.19864056E+00 -2.03643410E-03 6.52040211E-06 -5.48797062E-09 1.77197817E-12 -3.02937267E+04 -8.49032208E-01 -2.90848168E+04];
H1H2O = (cH2O(1) + cH2O(2).*T1H2O/2 + cH2O(3).*T1H2O.^2/3 + cH2O(4).*T1H2O.^3/4 + cH2O(5).*T1H2O.^4/5 + cH2O(6)./T1H2O)*R.*T1H2O;
H2H2O = (cH2O(8) + cH2O(9).*T2H2O/2 + cH2O(10).*T2H2O.^2/3 + cH2O(11).*T2H2O.^3/4 + cH2O(12).*T2H2O.^4/5 + cH2O(13)./T2H2O)*R.*T2H2O;


T1CO = 200:1:1000;
T2CO = 1000:1:6000;
MCO = 28.01040;
cCO = [3.04848583E+00 1.35172818E-03 -4.85794075E-07 7.88536486E-11 -4.69807489E-15 -1.42661171E+04 6.01709790E+00 3.57953347E+00 -6.10353680E-04 1.01681433E-06 9.07005884E-10 -9.04424499E-13 -1.43440860E+04 3.50840928E+00 -1.32936276E+04];
H1CO = (cCO(1) + cCO(2).*T1CO/2 + cCO(3).*T1CO.^2/3 + cCO(4).*T1CO.^3/4 + cCO(5).*T1CO.^4/5 + cCO(6)./T1CO)*R.*T1CO;
H2CO = (cCO(8) + cCO(9).*T2CO/2 + cCO(10).*T2CO.^2/3 + cCO(11).*T2CO.^3/4 + cCO(12).*T2CO.^4/5 + cCO(13)./T2CO)*R.*T2CO;

T1CO2 = 200:1:1000;
T2CO2 = 1000:1:6000;
MCO2 = 44.00980;
cCO2 = [4.63659493E+00 2.74131991E-03 -9.95828531E-07 1.60373011E-10 -9.16103468E-15 -4.90249341E+04 -1.93534855E+00 2.35677352E+00 8.98459677E-03 -7.12356269E-06 2.45919022E-09 -1.43699548E-13 -4.83719697E+04 9.90105222E+00 -4.73281047E+04];
H1CO2 = (cCO2(1) + cCO2(2).*T1CO2/2 + cCO2(3).*T1CO2.^2/3 + cCO2(4).*T1CO2.^3/4 + cCO2(5).*T1CO2.^4/5 + cCO2(6)./T1CO2)*R.*T1CO2;
H2CO2 = (cCO2(8) + cCO2(9).*T2CO2/2 + cCO2(10).*T2CO2.^2/3 + cCO2(11).*T2CO2.^3/4 + cCO2(12).*T2CO2.^4/5 + cCO2(13)./T2CO2)*R.*T2CO2;


%Actual equation with dissociated products: CH4 + nO2 -> aCO2 + 2H2O + cCO + dO2
% 1) CO + 0.5O2 <-> CO2 + 0.5O2
m = 6;
z1 = sym('z1');
D1 = [1 m 0; -z1 -m*z1 z1; 1-z1 m*(1-z1) z1]; % [Initial no. of moles; Change in no. of moles; No. of moles at equilibrium]
%   CO O2 CO2

eqn1 = sum(D1(3, :)); %No. of moles in dissociation equation x at equilibrium in terms of zx

%Actual equation with dissociated products in terms of z: CH4 + nO2 -> z1CO2 + 2H2O + (1-z1)CO + 0.5(1-z1)O2

yCO = D1(3, 1)/eqn1; %Equilibrium mole fractions of each product
yO2 = D1(3, 2)/eqn1;
yCO2 = D1(3, 3)/eqn1;

K1 = (yCO2)/(yCO*yO2^0.5) * (Pc/Po)^-0.5; %Law of mass action
CK1 = [2.00E+08, 0.7, 12000];
T = sym('T');
KA = CK1(1)*T^CK1(2)*exp(-CK1(3)/(R*T));

Tp(1) = 1050;

for n = 2:2
    Ta = Tp(n-1);
    for j = 1:10000
        if j == 1
            Ta = Tp(n-1);
        else
            Ta = Temperature(j);
        end
        eqna = K1 == subs(KA, Ta);
        solz1(j) = vpasolve(eqna, z1);
        if Ta > 1000
        HP(j) = 2*H2H2O(find(T2H2O == Ta)) + subs(D1(3, 1), solz1(j))*H2CO(find(T2CO == Ta)) + subs(D1(3, 2), solz1(j))*H2O2(find(T2O2 == Ta)) + subs(D1(3, 3), solz1(j))*H2CO2(find(T2CO2 == Ta));
        else
        HP(j) = 2*H1H2O(find(T1H2O == Ta)) + subs(D1(3, 1), solz1(j))*H1CO(find(T1CO == Ta)) + subs(D1(3, 2), solz1(j))*H1O2(find(T1O2 == Ta)) + subs(D1(3, 3), solz1(j))*H1CO2(find(T1CO2 == Ta));
        end
        HR = 1*(H1CH4(find(T1CH4 == 298))) + m*(H1O2(find(T1O2 == 298)));
        if HP(j) > HR && (j == 1 || HP(j-1) > HR)
    %fprintf('HP too high at %d K, try lower temperature \n', Ta);
    Temperature(j+1) = Ta - 1;
    elseif HP(j) < HR && (j == 1 || HP(j-1) < HR)
    %fprintf('HP too low at %d K, try higher temperature \n', Ta);   
    Temperature(j+1) = Ta + 1;
    elseif (HP(j) < HR && HP(j-1) > HR) || (HP(j) > HR && HP(j-1) < HR)
    fprintf('Values found either side of HR, interpolating \n'); 
    Tp(n) = ((HR - HP(j-1)))/(HP(j) - HP(j-1)) * (Temperature(j) - Temperature(j-1)) + Temperature(j-1);
    fprintf('Solution found, Tc = %d K \n', Tp(n));
    break
        end
    end
end
