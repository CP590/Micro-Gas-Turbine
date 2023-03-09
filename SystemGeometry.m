syms x1




%Turbine casing
rTV = 15; %270°
tTV = 5;
arTV = 70;
alTV = 30; %Must be less than arTV
lTV = 15;
tHSBH = 5;
toTV = 5;



%Turbine impeller
rTI1 = 5;
rTI2 = 30;
lTI1 = 22.65;
lTI2 = 2;
y1 = rTI1 + 1.2^(x1-5);
y2 = rTI2;
gS1 = 3;

%Heat shield
lHS = lTV - lTI2 - gS1;
tHS = 3;

%Air gap 2
a1AG2 = 20; %Must be less than alTV
a2AG2 = 10; %Must be less than a1AG2 and more than rS1

%Bearing housing assembly 1 (Bearing housing, shaft, oil, housing-bearing
%gap, bearing)
rS1 = 5;
lS1 = 20;
lB = 10;
rB = 3.2; %Must be more than rS2
tB = rS1 - rB;
cS1BHA = 0.2;
rBHA = rS1 + cS1BHA; %Must be more than rS1
tBHO = 2;
dhO = 2;
l1O1 = 2;
l2O1 = 6;
rlO1 = 4;
riO1 = rBHA + tBHO;
roO1 = riO1 + rlO1;
lBHTV = 5;

%Bearing housing assembly 2 (Bearing housing, shaft, oil, air gap)
rS2 = 3; %Must be less than rS1 and rB
lS2 = 10;
aiBH = 25; %Must be less than alTV + tHSBH
aoBH = 50; %Probably greater than alTV + tHSBH + toTV
lUC = 5; %Must be less than lS2

%Bearing housing assembly 3 (Bearing housing, shaft, oil, housing-bearing
%gap, bearing)
lS3 = 20;
eBP = 5;
tBP1 = 5;
tBP2 = 5;
tBP3 = 5;

%Compressor impeller
rCI2 = rTI2;
rCI1 = rTI1;
lCI1 = lTI1;
lCI2 = lTI2;

%Compressor case
arCV = 60;
alCV = 30; %Must be less than arCV
rCV = 10; %270°
tCV = 5;
eCV = 10;

AiTV = 0.75*4*pi^2*rTV*arTV; %CV R from exhaust gas
AoTV = 0.75*4*pi^2*(rTV+tTV)*arTV + 2*pi*(alTV+tHSBH+toTV)*(tHS+lBHTV); %CV R to ambient
VTV = (0.75*2*pi^2*arTV*((rTV+tTV)^2-rTV^2)) + pi*lTV*((alTV+toTV+tHSBH)^2-alTV^2) + pi*(tHS + lBHTV)*((alTV+tHSBH)^2-alTV^2);


dy1 = diff(y1, x1);
y1s = (y1*sqrt(1+(dy1)^2));
AsTI1 = 2*pi*int(y1s, x1, 0, lTI1);
AsTI2 = 2*pi*rTI2*lTI2;
AsTI3 = pi*rTI1^2;
AsTI4 = pi*(rTI2^2 - rS1^2);
AsTI = vpa(AsTI1 + AsTI2 + AsTI3 + AsTI4); %CV R from exhaust gas
VTI1 = pi*int(y1^2, x1, 0, lTI1);
VTI2 = pi*y2^2*lTI2;
VTI = vpa(VTI1 + VTI2);
ATIS1 = pi*rS1^2; %CD TI to S1

AHS = pi*(alTV^2 - rBHA^2); %CV R from exhaust gas
ATVHS = 2*pi*(alTV+tHSBH)*tHS; %CD TV to HS
VHS = AHS*tHS + (lHS-tHS)*(alTV^2 - (alTV-tHS)^2) + tHS*((alTV+tHSBH)^2 - alTV^2);

AHSAG = pi*((alTV-tHS)^2-rS1^2) + 2*pi*(alTV-tHS)*lHS; % CV R HS to air gap
ABH1AG = 2*pi*0.5*lHS*(a2AG2 + a1AG2) + pi*(a1AG2^2-a2AG2^2) + pi*((alTV-tHS)^2-a1AG2^2); %CV R air gap to BH1
VAG = pi*(a1AG2^2-a2AG2^2)*0.5*lHS + pi*((alTV-tHS)^2-a1AG2^2)*lHS;

ABH1Am = pi*((alTV+tHSBH)^2 - aiBH^2); %CV R BH1 to ambient
ABH1TV = 2*pi*(alTV+tHSBH)*lBHTV; %CD TV to BH1
ABH1HS = pi*((alTV+tHSBH)^2 - (alTV-tHS)^2) + pi*(a2AG2^2 - rBHA^2); %CD HS to BH1;
ABH1ABH1B1 = 2*pi*rBHA*lB - 4*pi*(0.5*dhO)^2; %CV R BH1 to air BH1-B1
AB1ABH1B1 = 2*pi*rS1*lB - 4*pi*(0.5*dhO)^2; %CV R air BH1-B1 to B1
AS1O1 = 2*pi*rS2*lB; %CV R S1 to O1
AO1B1 = 2*pi*rB*lB - 4*pi*(0.5*dhO)^2 + 4*pi*dhO*tB; %CV R O1 to B1
ABH1O1 = pi*dhO*tBHO + 2*pi*rlO1*(l1O1 + l2O1); %CV R BH1 to O1
AS1 = 2*pi*rS1*gS1; %CV R exhaust gas to S1
VS1 = pi*rS1^2*(lS1 - lB) + pi*rS2*lB;
VB = pi*(rS1^2 - rB^2)*lB - 4*pi*(0.5*dhO)^2*tB;
VBH1 = pi*(a2AG2^2-rBHA^2)*(lS1-gS1-tHS) + pi*(a1AG2^2-rBHA^2)*(lS1-gS1-tHS-0.5*lHS) + pi*((alTV-tHS)^2-rBHA^2)*lBHTV - pi*rlO1^2*(l1O1+l2O1+dhO) - pi*(0.5*dhO)^2*tBHO;
VABHB = pi*(rBHA^2-rS1^2)*lB - pi*(0.5*dhO)^2*(rBHA-rS1);
VO1 = pi*(rB^2-rS2^2)*lB + 4*pi*(0.5*dhO)^2*(tB+tBHO) + pi*rlO1^2*(l1O1+l2O1+dhO);

AS1S2 = pi*rS2^2; %CD S1 to S2;
AS2AS2BH2 = AS1S2*lS2; %CV R S2 to AS2BH2
AAS2BH2B = pi*(rS1^2 - rB^2); %CV R AS2BH2 to B (1 & 2)
AAS2BH2BH2 = pi*rBHA^2*lS2; %CV R AS2BH2 to BH2
ABH1BH2 = pi*(aiBH^2-rBHA^2) - pi*rlO1^2; %CD BH1 to BH2
ABH2O2 = pi*rlO1^2*lS2; %CV R BH2 to O2
ABH2Am = 2*pi*aiBH*lUC + 2*pi*aoBH*(lS2-lUC) + pi*(aoBH^2-aiBH^2); %CV R BH2 to ambient
VS2 = pi*rS2*lS2;
VAS2BH2 = pi*(rBHA^2-rS2^2)*lS2;
VO2 = pi*rlO1^2*lS2;
VBH2 = pi*(aiBH^2-rBHA^2)*lS2 + pi*(aoBH^2-aiBH^2)*(lS2-lUC) - VO2;

AS2S3 = pi*rS2^2; %CD S2 to S3;
ABH2BH3 = pi*aoBH^2-rBHA^2 - pi*rlO1^2; %CD BH2 to BH3
ABH3ABH3B2 = 2*pi*rBHA*lB - 4*pi*(0.5*dhO)^2; %CV R BH3 to air BH3-B2
AB2ABH3B2 = 2*pi*rS1*lB - 4*pi*(0.5*dhO)^2; %CV R air BH3-B2 to B2
AS3O3 = 2*pi*rS2*lB; %CV R S3 to O3
AO3B2 = 2*pi*rB*lB - 4*pi*(0.5*dhO)^2 + 4*pi*dhO*tB; %CV R O3 to B2
ABH3O3 = pi*dhO*tBHO + 2*pi*rlO1*(l1O1 + l2O1); %CV R BH3 to O3
ABH3Am = 2*pi*aoBH*(lS3-tBP2-lUC) + 2*pi*aiBH*lUC + pi*(aoBH^2-aiBH^2) + pi*((aoBH+eBP)^2-aiBH^2); %CV R BH2 to ambient
ABH3Air = pi*((rCI2+tBP3)^2-rCI2^2); %CV R BH3 to discharge air
VO3 = VO1;
VS3 = pi*rS2^2*lS3;
VBH3 = pi*(rBHA^2-rB^2)*(lS3-lB) + pi*(aiBH^2-rBHA^2)*lS3 + pi*(aoBH^2-aiBH^2)*(lS3-tBP2-lUC) + pi*(rCI2^2-aiBH^2)*tBP2 + pi*((rCI2+tBP3)^2-rCI2^2)*(tBP2+tBP1) + pi*(aoBH^2-(rCI2+tBP3)^2)*tBP2 - pi*rlO1^2*(l1O1+l2O1+dhO) - pi*(0.5*dhO)^2*tBHO; 

ABH3CV = 2*pi*(rCI2+tBP3)*tBP1 + 2*pi*(rCI2+tBP3+eBP)*tBP2 + pi*((rCI2+tBP3+eBP)^2-(rCI2+tBP3)^2); %CD BH3 to CV
AoCV = 0.75*4*pi^2*(rCV+tCV)*arCV; %CV R CV to ambient
AiCV = 0.75*4*pi^2*rCV*arCV; %CV R CV to discharge air
VCV = (0.75*2*pi^2*arCV*((rCV+tCV)^2-rCV^2)) + pi*((aoBH+eBP+eCV)^2-(aoBH+eBP)^2)*tBP2;

AsCI1 = AsTI1;
AsCI2 = AsTI2;
AsCI3 = AsTI3;
AsCI4 = pi*(rTI2^2 - rS2^2);
AsCI = vpa(AsCI1 + AsTI2 + AsCI3 + AsCI4); %CV R to discharge air
VCI = VTI;
ACIS3 = pi*rS2^2; %CD S3 to CI