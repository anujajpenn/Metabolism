% This script takes the output of 3D_cell.m and calculates optimum cell
% shape with and without considering interfacial energy as well as 
% adenosine nucleotide concentrations, and other quantities of interest 

% Quantities ending in "1" represent quantities calculated without
% considering interfacial energy. For example, the variable "TotalEnergy1"
% is the sum of all energy contributions without interfacial energy (or in
% other words, the volumetric part of the metabolic potential), whereas the
% variable "TotalEnergy" = "TotalEnergy1" + interfacial energy.
AR = 1:1:15;
Density = logspace(-0.5,1,30);
RefinedAR = 1:0.02:15;
RefinedDensity = logspace(-0.5,1,120);
g = 0.18; 
gamma0 = 2.95e-4; 
TotalEnergy1 = Energy;
Elongation = zeros(length(AR), 1);
cell_surface_area = zeros(length(AR),1);
Volume = 5e-15;

for i = 1:length(AR)
    b = (Volume/(4*pi/3)/AR(i))^(1/3);
    a = b * AR(i);
    c = b;
    S = 4*pi*(((a*b)^1.6075+(a*c)^1.6075+(b*c)^1.6075)/3)^(1/1.6075);
    circ = (4*pi^2*a*b)/((pi*(3*(a+b)-sqrt((3*a+b)*(a+3*b))))^2);
    cell_surface_area(i) = S;
    Elongation(i) = AR(i)/circ;
end

RefinedElongation = zeros(length(RefinedAR), 1);
RefinedSpreadArea = zeros(length(RefinedDensity), length(RefinedAR));
for i = 1:length(RefinedAR)
    b = (Volume/(4*pi/3)/RefinedAR(i))^(1/3);
    a = b * RefinedAR(i);
    circ = (4*pi^2*a*b)/((pi*(3*(a+b)-sqrt((3*a+b)*(a+3*b))))^2);
    RefinedElongation(i) = RefinedAR(i)/circ;
    RefinedSpreadArea(:, i) = pi*a*b;
end

% Refine quantities through interpolation
[AM, DM] = meshgrid(Elongation, Density);
[RAM, RDM] = meshgrid(RefinedElongation, RefinedDensity);

RefinedTotalEnergy1 = interp2(AM, DM, TotalEnergy1, RAM, RDM, 'spline');
RefinedRho = interp2(AM, DM, Rho, RAM, RDM, 'spline');
RefinedS = interp2(AM, DM, Stressvol, RAM, RDM, 'spline');
RefinedUmech = interp2(AM, DM, Umech, RAM, RDM, 'spline');
RefinedUbinding = interp2(AM, DM, Ubinding, RAM, RDM, 'spline');
RefinedUmotor = interp2(AM, DM, Umotor, RAM, RDM, 'spline');
RefinedUATP = interp2(AM, DM, UATP, RAM, RDM, 'spline');

%Find optimum path w.r.t. cell AR and mtrix density
MinAR1 = zeros(length(RefinedDensity), 1);
MinDS1 = zeros(length(RefinedAR), 1);
for i = 1:length(RefinedDensity)
    curmin = 100;
    curminpos = 0;
    for j = 1:length(RefinedAR)
        if RefinedTotalEnergy1(i, j) < curmin
            curminpos = j;
            curmin = RefinedTotalEnergy1(i, j);
        end
    end
    MinAR1(i) = curminpos;
end
for i = 1:length(RefinedAR)
    curmin = 100;
    curminpos = 0;
    for j = 1:length(RefinedDensity)
        if RefinedTotalEnergy1(j, i) < curmin
            curminpos = i;
            curmin = RefinedTotalEnergy1(j, i);
        end
    end
    MinDS1(i) = curminpos;
end

% Find quantities of interest using the optimum path without adding 
% interfacial energy
Contractility1 = zeros(length(MinAR1));
for i = 1:length(MinAR1)
    Contractility1(i) = RefinedRho(i, MinAR1(i));
end
Stress1 = zeros(length(MinAR1));
for i = 1:length(MinAR1)
    Stress1(i) = RefinedS(i, MinAR1(i));
end
SpreadArea1 = zeros(length(MinAR1));
for i = 1:length(MinAR1)
    % Converting to microns here
    SpreadArea1(i) = (RefinedSpreadArea(i, MinAR1(i)))/((1e-6)^2);
end
MechEn1 = zeros(length(MinAR1));
for i = 1:length(MinAR1)
    MechEn1(i) = (RefinedUmech(i, MinAR1(i)));
end
BindingEn1 = zeros(length(MinAR1));
for i = 1:length(MinAR1)
    BindingEn1(i) = (RefinedUbinding(i, MinAR1(i)));
end
ATPHydEn1 = zeros(length(MinAR1));
for i = 1:length(MinAR1)
    ATPHydEn1(i) = (RefinedUATP(i, MinAR1(i)));
end
MotorEn1 = zeros(length(MinAR1));
for i = 1:length(MinAR1)
    MotorEn1(i) = (RefinedUmotor(i, MinAR1(i)));
end

figure
% Plot landscape and optimum path based on volumetric part of the metabolic
% potential
hold on;
contourf(RefinedDensity, RefinedElongation, transpose(RefinedTotalEnergy1), 16, 'k');
plot(RefinedDensity, RefinedElongation(MinAR1));
hold off;

figure
hold on;
% Plot stress and contractility over optimum path
yyaxis right;
plot(RefinedDensity, Contractility1(:,1)/Volume);
yyaxis left;
plot(RefinedDensity, Stress1(:,1)/Volume);
hold off;


% Now repeat calculations considering interfacial energy
TotalEnergy = Energy;
InterfacialEnergy = Energy;
for i = 1:length(AR)
    b = (Volume/(4*pi/3)/AR(i))^(1/3);
    a = b * AR(i);
    c = b;
    S = 4*pi*(((a*b)^1.6075+(a*c)^1.6075+(b*c)^1.6075)/3)^(1/1.6075);
    circ = (4*pi^2*a*b)/((pi*(3*(a+b)-sqrt((3*a+b)*(a+3*b))))^2);
    Elongation(i) = AR(i)/circ;
    for j = 1:length(Density)
        InterfacialEnergy(j, i) = S*gamma0*Density(j)^g;
        TotalEnergy(j, i) =  Energy(j, i) +  S*gamma0*Density(j)^g;
    end
end

RefinedEnergy = interp2(AM, DM, Energy, RAM, RDM, 'spline');
RefinedTotalEnergy = interp2(AM, DM, TotalEnergy, RAM, RDM, 'spline');
RefinedInterfacialEnergy = interp2(AM, DM, InterfacialEnergy, RAM, RDM, 'spline');
RefinedCellStrain = interp2(AM, DM, Cell_Strain, RAM, RDM, 'spline');
RefinedMatStrain = interp2(AM, DM, Mat_Strain, RAM, RDM, 'spline');
RefinedMatUmech = interp2(AM, DM, Mat_Umech, RAM, RDM, 'spline');

% Find optimum path now taking interfacial energy into account
MinAR = zeros(length(RefinedDensity), 1);
MinDS = zeros(length(RefinedAR), 1);
for i = 1:length(RefinedDensity)
    curmin = 100;
    curminpos = 0;
    for j = 1:length(RefinedAR)
        if RefinedTotalEnergy(i, j) < curmin
            curminpos = j;
            curmin = RefinedTotalEnergy(i, j);
        end
    end
    MinAR(i) = curminpos;
end
for i = 1:length(RefinedAR)
    curmin = 100;
    curminpos = 0;
    for j = 1:length(RefinedDensity)
        if RefinedTotalEnergy(j, i) < curmin
            curminpos = i;
            curmin = RefinedTotalEnergy(j, i);
        end
    end
    MinDS(i) = curminpos;
end

Contractility = zeros(length(MinAR));
for i = 1:length(MinAR)
    Contractility(i) = RefinedRho(i, MinAR(i));
end
Stress = zeros(length(MinAR));
for i = 1:length(MinAR)
    Stress(i) = RefinedS(i, MinAR(i));
end
SpreadArea = zeros(length(MinAR));
for i = 1:length(MinAR)
    SpreadArea(i) = (RefinedSpreadArea(i, MinAR(i)))/((1e-6)^2);
end
MechEn = zeros(length(MinAR));
for i = 1:length(MinAR)
    MechEn(i) = (RefinedUmech(i, MinAR(i)));
end
BindingEn = zeros(length(MinAR));
for i = 1:length(MinAR)
    BindingEn(i) = (RefinedUbinding(i, MinAR(i)));
end
ATPHydEn = zeros(length(MinAR));
for i = 1:length(MinAR)
    ATPHydEn(i) = (RefinedUATP(i, MinAR(i)));
end
MotorEn = zeros(length(MinAR));
for i = 1:length(MinAR)
    MotorEn(i) = (RefinedUmotor(i, MinAR(i)));
end
SpreadArea = zeros(length(MinAR));
for i = 1:length(MinAR)
    SpreadArea(i) = (RefinedSpreadArea(i, MinAR(i)))/((1e-6)^2);
end
CStrain = zeros(length(MinAR));
for i = 1:length(MinAR)
    CStrain(i) = RefinedCellStrain(i, MinAR(i));
end
MStrain = zeros(length(MinAR));
for i = 1:length(MinAR)
    MStrain(i) = RefinedMatStrain(i, MinAR(i));
end
MUmech = zeros(length(MinAR));
for i = 1:length(MinAR)
    MUmech(i) = RefinedMatUmech(i, MinAR(i));
end



% Find ATP, ADP, AMP, and AMPK values based on the volumetric stress found
% in optimum shaped cells at different matrix densities.

% these constants come from the fit of our model to experimental data in 
% 
C = 0.00000006776; 
gamma_prime = 0.3; 
n = 0.5284; 
m = 2*n; % Exponent relating strength of calcium signaling (m = n + l, 
% assuming l is twice n).
x = zeros(length(RefinedDensity),length(RefinedAR));
ATP = zeros(length(RefinedDensity),1);
ADP = zeros(length(RefinedDensity),1);
AMP = zeros(length(RefinedDensity),1);
Aratio = zeros(length(RefinedDensity),1);
AMPK = zeros(length(RefinedDensity),1);
for i = 1:length(RefinedDensity)
    % Note the following concs are normalized to [Atot]. In order to find
    % actual concentrations, a value for [Atot] would have to be assumed.
    x(i) = C * (1 + gamma_prime)^2 / (Stress(i,1)^n);
    ATP(i) = (x(i) + 1 - sqrt(x(i)^2 + 2 * x(i)));
    ADP(i) = ((1 - ATP(i))/(1 + gamma_prime));
    AMP(i) = 1 - ATP(i) - ADP(i);
    Aratio(i) = ATP(i) / ADP(i);
    AMPK(i) = AMP(i) * Stress(i,1)^m; % Normalized to [Atot] and [AMPKtot]
end

IntEn = zeros(length(MinAR));
for i = 1:length(MinAR)
    IntEn(i) = RefinedInterfacialEnergy(i, MinAR(i));
end

figure
% Plot energy landscape and optimum path
hold on;
contourf(RefinedDensity, RefinedElongation([1, 187],1), transpose(RefinedTotalEnergy(:,[1, 187])), 25, 'k');
plot(RefinedDensity, RefinedElongation(MinAR));
xlim([0.3,7]);
ylim([1,10]);
hold off;

Elong = RefinedElongation(MinAR);
% Elongation vs collagen density
figure
hold on;
plot([1,2,3,4],[Elong(17),Elong(41),Elong(65),Elong(96)]);
Elongexp = [1.20 1.90 2.07 2.84 7.60; 1.35 1.93 4.51 14.51 29.24; 1.20 1.67 4.11 8.91 19.63; 1.13 1.38 1.41 1.56 7.16];
Elongexp = Elongexp(:, [1 2 2 3 4 4 5]);
boxplot(Elongexp.', 'Whisker', inf);
hold off;

figure
% Stress and contractility in optimum shapes
hold on;
yyaxis right;
plot(RefinedDensity, Contractility(:,1)/Volume);
yyaxis left;
plot(RefinedDensity, Stress(:,1)/Volume);
hold off;


% AMPK - uncomment to plot - this assumes the experimentally measured 
% values for pAMPK have been loaded into a variable called p_AMPK_3D
%figure
%hold on;
%boxplot(p_AMPK_3D);
%plot([1,2,3], [AMPK(41)/AMPK(96),AMPK(79)/AMPK(96),AMPK(96)/AMPK(96)]);
%hold off;

% Plot ATP:ADP ratio 
figure
hold on;
plot([1,2,3,4], [Aratio(17,1), Aratio(41,1), Aratio(65,1), Aratio(96,1)]);
ATR_norm = 0.619;
ATR = [0.155 0.462 0.619 0.812 1.108; 0.203 0.377 0.614 1.039 1.712; 0.214 0.577 0.724 1.076 2.259; 0.366 0.614 0.849 1.148 2.037]/ATR_norm;
ATR = ATR(:, [1 2 2 3 4 4 5]);
boxplot(ATR.', 'Whisker', inf);
%bar([1,2,3,4], [0.62/0.62,0.7/0.62,0.89/0.62,0.94/0.62]); % Wu et al., Biophys J (2021) Fig 4E*
%errorbar([1,2,3,4], [0.62/0.62,0.7/0.62,0.89/0.62,0.94/0.62], [0.297,0.551,0.495,0.475], [0.193,0.339,0.365,0.375]);
hold off;

% AMPK vs glucose uptake (data from Zanotelli et al., Mol Bio Cell (2018), Figure 3)
figure
hold on;
plot([1,2,3], [AMPK(41,1)/AMPK(96,1), AMPK(79,1)/AMPK(96,1), AMPK(96,1)/AMPK(96,1)]);
bar([1,2,3], [0.7197,0.8252,1]);
errorbar([1,2,3], [0.7197,0.8252,1], [0,0,0], [0.03615,0.05273,0.05399]);
hold off;

% AMPK vs TMRM (data from Wu et al., Biophys J (2021) Fig 4F)
TMRM_norm = 2295.65;
TMRM = [495.65 843.48 1052.17 1426.09 3686.96; 530.43 913.04 1321.74 2469.56 5756.52; 591.30 1060.87 1452.17 2582.61 5921.74; 695.65 1443.48 2295.65 3121.74 5930.43]/TMRM_norm;
TMRM = TMRM(:, [1 2 2 3 4 4 5]);
figure
hold on;
boxplot(TMRM.', 'Whisker', inf);
plot([1,2,3,4], [AMPK(17,1)/AMPK(96,1), AMPK(41,1)/AMPK(96,1), AMPK(65,1)/AMPK(96,1), AMPK(96,1)/AMPK(96,1)]);
hold off;

% Lower alpha (ML7 treatment)
alpha = [2.33e-3, 1e-5]; % Alpha values used corresponding to control and ML7 
stress_alpha = [1.295e-12, 7.198e-13]; % Used stress in optimum cell shape in 1.5 mg/ml collagen (from 3Dcell.mph - more accurate than the interpolated values calculated via matlab link which does not run at exactly 1.5 mg/ml)
x_alpha = zeros(length(stress_alpha),1);
ATP_alpha = zeros(length(stress_alpha),1);
ADP_alpha = zeros(length(stress_alpha),1);
AMP_alpha = zeros(length(stress_alpha),1);
Aratio_alpha = zeros(length(stress_alpha),1);
AMPK_alpha = zeros(length(stress_alpha),1);
for i = 1:length(stress_alpha)
    x_alpha(i) = C * (1 + gamma_prime)^2 / (stress_alpha(i)^n);
    ATP_alpha(i) = (x_alpha(i) + 1 - sqrt(x_alpha(i)^2 + 2 * x_alpha(i)));
    ADP_alpha(i) = ((1 - ATP_alpha(i))/(1 + gamma_prime));
    AMP_alpha(i) = 1 - ATP_alpha(i) - ADP_alpha(i);
    Aratio_alpha(i) = ATP_alpha(i) / ADP_alpha(i);
    AMPK_alpha(i) = AMP_alpha(i) * stress_alpha(i)^m;
end
alpha_data = [[0.7071 0.6158]; [0.8017 0.7974]; [1.1983 0.9876]; [1.5368 1.2273]]; % (data from Zanotelli et al., Mol Bio Cell (2018), Figure 5
alpha_norm = 0.9305; % normalizing to median of control data
figure
hold on;
boxplot(alpha_data);
plot([1,2], [Aratio_alpha(1)/Aratio_alpha(1), Aratio_alpha(2)/Aratio_alpha(1)]);
hold off;

% ATP consumption rate (data from Zanotelli et al., Mol Bio Cell (2018), Figure 3)
alpha_d = 2.33e-3;
figure
hold on;
plot([1,2,3],[Stress(41,1)/Stress(96,1), Stress(65,1)/Stress(96,1), Stress(96,1)/Stress(96,1)])
bar([1,2,3], [0.000322,0.000378,0.000533]/0.000533); % Normalized to measurement at 5 mg/ml
errorbar([1,2,3], [0.000322,0.000378,0.000533]/0.000533, [0,0,0], [0.000378-0.000322,0.000432-0.000378,0.000599-0.000533]/0.000533);
hold off;



