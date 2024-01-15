
%Behrad Bolkhari     Student Number:9919653       Dr Mahdi Imanian
clc
clear
close all
format long
q = 1000;  %bbl/day
api = 16; %degree API
mu = 5; %cp
alpha = 3; %degree
L = 1000; %ft
d = 2.259; %in
D = d / 12; %ft
g=32.17;
Delta_u=0
relative_roughness = 0.001
g_o = 141.5/(api+131.5)
Density=62.4*g_o

alpha = deg2rad(alpha);
dz = L * cos(alpha)
u = 4 * q * 5.615 / (86400 * pi * D^2)
Re = 1.48 * q * Density / (mu * D)
f_F = 1 / (16 * log10(relative_roughness/ 3.7065 - 5.0452 / Re * log10(relative_roughness^1.1098 / 2.8257 + (7.149 / Re)^0.8981)))^2
delta_p =(32.17/32.17)*(Density*dz)+(Density*Delta_u^2)/(2*37.14)+(2*f_F*Density*u^2*L)/(32.17*D)
