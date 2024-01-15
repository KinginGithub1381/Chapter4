%Behrad Bolkhari     Student Number:9919653       Dr Mahdi Imanian
clc
clear
close all
L = 8000; % ft
d=3.5;%in
D = d / 12; % ft
api=30; % API
mu = 2; % cp
alpha=0%radyan
glr = 500; % scf/bbl
gamma_g = 0.75;
p_head = 1000; % psia
T_head = 120; % F
T_bot = 180; % F
relative_roughness=0.0006;
liquid_rate = 1500; % bbl/day
water_cut = 0.2; % fraction
ift = 30; % dynes/cm
gamma_w = 1.05;
qw=0;
qsc=3000;%MSCF/Day
g_o = 141.5/(api+131.5);
gamma_l = ((liquid_rate - qw) * g_o + qw * gamma_w) / liquid_rate;
qg = glr * liquid_rate;
A = 3.141592653 * (D^2) / 4;
Ppc = 678 - 50 *(gamma_g-0.5);
Tpc = 326 +315.7 *(gamma_g-0.5);
Tpr=T_head/Tpc
Ppr=p_head/Ppc
F = 0.3106 - 0.49 * Tpr + 0.1824 * Tpr^2;
E = 9 * (Tpr -1);
D = 10^F;
C = 0.132 - 0.32 * log10(Tpr);
B = (0.62 - 0.23 * Tpr) * Ppr + (0.066 / (Tpr - 0.86) - 0.037) * Ppr^2 + 0.32 * Ppr^6 / 10^E;
A = 1.39 * (Tpr - 0.92)^(0.5) - 0.36 * Tpr - 0.10;
z_wellead = A + (1 - A) / exp(B) + C * Ppr^D;
f_F = 1 / ( 1.74-2*log10(2*relative_roughness))^2
L=[0 1000 2000 3000 4000 5000 6000 7000 8000];
alpha = deg2rad(alpha);
s =0.0375* L * cos(alpha)*gamma_g/(z_wellead*(460+T_head))
Pwf=sqrt(exp(s)+(6.67*10^-4)*(exp(s)-1)*f_F*(qsc^2)*T_head^2*z_wellead^2)