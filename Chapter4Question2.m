%Behrad Bolkhari     Student Number:9919653       Dr Mahdi Imanian
clc
clear
close all
format long
D = 1.66;
P_wellhead = 300; % psia
production_rate = 2000; % STB/day
gas_liquid_ratio = 800; % scf/STB
wc = 0.30;
api = 40; 
gamma_w = 1.05;
gamma_g = 0.7;
bw = 1.2; % bbl/STB
wellhead_temperature = 100; % F
L = 8000; % ft
bottomhole_temperature = 170; % F

gor = gas_liquid_ratio / (1 - wc)
wor = wc / (1 - wc)
qo = production_rate * (1 - wc)
yH2S=0, yCO2=0, yN2=0
gamma_o = 141.5/(api+131.5)
M = 350.17 * (gamma_o + wor * gamma_w) + 0.0765 * gor * gamma_g
Rs = gamma_g * (P_wellhead / 18 * 10^(0.0125 * api) / 10^(0.00091 * wellhead_temperature))^(1.2048)
Bo = 0.9759 + 0.00012 * (Rs * ((gamma_g/gamma_o)^(0.5)) + 1.25 * wellhead_temperature)^(1.2)
Ppc = 678 - 50 * (gamma_g - 0.5) - 206.7 * yN2 + 440 * yCO2 + 606.7 * yH2S
Tpc = 326 + 315.7 * (gamma_g - 0.5) - 240 * yN2 - 83.3 * yCO2 + 133.3 * yH2S
Tpr = (wellhead_temperature + 460) / Tpc;
Ppr = P_wellhead / Ppc;

F = 0.3106 - 0.49 * Tpr + 0.1824 * Tpr^2;
E = 9 * (Tpr -1);
D = 10^F;
C = 0.132 - 0.32 * log10(Tpr);
B = (0.62 - 0.23 * Tpr) * Ppr + (0.066 / (Tpr - 0.86) - 0.037) * Ppr^2 + 0.32 * Ppr^6 / 10^E;
A = 1.39 * (Tpr - 0.92)^(0.5) - 0.36 * Tpr - 0.10;
z_wellead = A + (1 - A) / exp(B) + C * Ppr^D;

V_Wellhead = 5.615 * (Bo + wor * bw) + (gor - Rs) * (14.7 /P_wellhead) * ((wellhead_temperature + 460) / 520) * z_wellead
g_wellhead = M / V_Wellhead;
D_Rov = 1.4737 * 10^(-5) * M * qo * 12 / D
f2F = 10^(1.444 - 2.5 * log10(D_Rov));
k = f2F * qo^2 * M^2 / (7.4137 * 10^10 * D^5)

function PBHP = finding_bhp(P_wellhead, rho_wellhead, k, L, api, bottomhole_temperature, gamma_g, wor, bw, gor, M)
    PBHP = P_wellhead;
    error_h = 1;
    for i = 1:10
        rho_bottomhole = calculating_rho_at_any_point(api, pbh, bottomhole_temperature, gamma_g, wor, bw, gor, M);    
        rho_avg = (rho_wellhead + rho_bottomhole) / 2;
        pbh = P_wellhead + (rho_avg + k / rho_avg) * L / 144;
        error_h = 144 * (PBHP - P_wellhead) / (rho_avg + k / rho_avg) - L;
    end
    fprintf('The bottomhole using Poettmann-Carpenter method is: %.2f psia\n', bhp);
end
