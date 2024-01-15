
%Behrad Bolkhari     Student Number:9919653       Dr Mahdi Imanian
clc
clear
close all
L = 6000; % ft
D = 1.995 / 12; % ft
api=30; % API
mu = 2; % cp
glr = 500; % scf/bbl
gamma_g = 0.65;
p_head = 100; % psia
T_head = 80; % F
T_bot = 140; % F
liquid_rate = 1500; % bbl/day
water_cut = 0.2; % fraction
ift = 30; % dynes/cm
gamma_w = 1.05;

qw = liquid_rate * water_cut;
g_o = 141.5/(api+131.5);
gamma_l = ((liquid_rate - qw) * g_o + qw * gamma_w) / liquid_rate;
qg = glr * liquid_rate;
A = 3.141592653 * (D^2) / 4;
dz = 100;
u_sl = liquid_rate * 5.615 / 86400 / A;
mu_l = (mu * (liquid_rate - qw) + 0.5 * qw) / liquid_rate;
m_t = gamma_l * 62.4 * liquid_rate * 5.615 + 0.0765 * gamma_g * qg;

function [p, depth] = Hagedorn_Brown_Correlation(L, p_head, T_head, T_bot, dz, gamma_g, gamma_l, ift, D, m_t, qg, A, mu_l, u_sl)
    p = [p_head];
    depth = [0];
    Tpc = calculate_Tpc_ahmed_method(gamma_g, 0, 0, 0);
    Ppc = calculate_Ppc_ahmed_method(gamma_g, 0, 0, 0);
    for i = 1:floor(L / dz)
        eps = 0.0006;
        T = depth(i) * (T_bot - T_head) / L + T_head;
        Tpr = (T + 460) / Tpc;
        Ppr = p(i) / Ppc;
        F = 0.3106 - 0.49 * Tpr + 0.1824 * Tpr^2;
        E = 9 * (Tpr -1);
        D = 10^F;
        C = 0.132 - 0.32 * log10(Tpr);
        B = (0.62 - 0.23 * Tpr) * Ppr + (0.066 / (Tpr - 0.86) - 0.037) * Ppr^2 + 0.32 * Ppr^6 / 10^E;
        A = 1.39 * (Tpr - 0.92)^(0.5) - 0.36 * Tpr - 0.10;
        z = A + (1 - A) / exp(B) + C * Ppr^D;
        z = calculate_zfactor_brill_and_beggs_method(Tpr, Ppr);
        mu_g = calculate_mu_carr_et_al_method(gamma_g, 0, 0, 0, T, Ppr, Tpr);
        u_sg = 1 / A * qg * z * (460 + T) / (460+60) * (14.7 / p(i)) /86400;
        N_vl = 1.938 * u_sl * (62.4 * gamma_l / ift)^0.25;
        N_vg = 1.938 * u_sg * (62.4 * gamma_l / ift)^0.25;
        N_d = 120.872 * D * sqrt(62.4 * gamma_l / ift);
        N_l = 0.15726 * mu_l * (1 / (62.4 * gamma_l * ift)^3)^0.25;
        x1 = log10(N_l) + 3;
        Y = - 2.69851 + 0.15840954 * x1 - 0.55099756 * x1^2 + 0.54784917 * x1^3 - 0.12194578 * x1^4;
        Cn_l = 10 ^ Y;
        x2 = N_vl / N_vg^0.575 * (p(i) / 14.7)^0.1 * Cn_l / N_d;
        y_L_psi = - 0.10306578 + 0.617774 * (log10(x2) + 6) - 0.632946 * (log10(x2) + 6)^2 + 0.29598 * (log10(x2) + 6)^3 - 0.0401 * (log10(x2) + 6)^4;
        x3 = 0.012;
        % x3 = N_vg * N_l^0.38 / N_d^2.14
        psi = 0.91162574 - 4.82175636 * x3 + 1232.25036621 * x3^2 - 22253.57617 * x3^3 + 116174.28125 * x3^4;
        y_L = psi * y_L_psi;
        Re = 0.022 * m_t / ((D * 12) * mu_l * y_L * mu_g * (1 - y_L));
        f = 1 / (- 4 * log10(eps / 3.7065 - 5.0452 / Re * log10(eps - 1.1098 / 2.8257 + (7.149 / Re) - 0.8981)))^2;
        rho_g = 28.97 * gamma_g * p(i) / z / 10.73 / (460 + T);
        rho_avg = y_L * gamma_l * 62.4 + (1 - y_L) * rho_g;
        p_gradient = 1 / 144 * (rho_avg + f * m_t^2 / 7.413 / 10000000000 / D^5 / rho_avg);
        p = [p, p(i) + p_gradient * dz]
        depth = [depth, depth(i) + dz]
    end
    plot(p_depth(1,:), p_depth(2,:))
set(gca,'YDir','reverse')
xlabel('Pressure (psia)')
ylabel('Depth (ft)')
grid on
ylim([max(p_depth(2,:)), 0])
end

