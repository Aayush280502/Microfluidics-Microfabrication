
clear; clc; close all;

rho_c = 1000;      % Density of continuous phase [kg/m^3]
rho_d = 900;       % Density of dispersed phase [kg/m^3]
mu_c  = 1e-3;     
mu_d  = 2e-3;      
sigma = 0.03;      

W = 100e-6;        
H = 100e-6;       
W_orifice = 50e-6; 
% Flow rate ratios to study
QcQd = linspace(1, 20, 10);  
Qd = 1e-9;                    


Droplet_D = zeros(size(QcQd));
Droplet_S = zeros(size(QcQd));
Droplet_f = zeros(size(QcQd));


for i = 1:length(QcQd)
    Qc = QcQd(i) * Qd;   
    Uc = Qc / (W * H);
    Ud = Qd / (W_orifice * H);
    Ca = mu_c * Uc / sigma;
    Droplet_D(i) = W * (0.5 + 1.5 * (Qd/Qc)^0.3) / (1 + 10*Ca);
    Droplet_S(i) = 3 * Droplet_D(i);
    Droplet_f(i) = Uc / (Droplet_D(i) + Droplet_S(i));
end


Droplet_D_um = Droplet_D * 1e6;
Droplet_S_um = Droplet_S * 1e6;


fprintf('\n----- DROPLET FORMATION RESULTS -----\n');
fprintf(' W_orifice = %.0f µm\n', W_orifice*1e6);
fprintf(' Qd = %.2e m^3/s\n\n', Qd);
fprintf(' Qc/Qd\tDiameter(µm)\tSpacing(µm)\tFrequency(Hz)\n');
for i = 1:length(QcQd)
    fprintf(' %.1f\t\t%.1f\t\t%.1f\t\t%.1f\n', ...
        QcQd(i), Droplet_D_um(i), Droplet_S_um(i), Droplet_f(i));
end


figure('Color','w');

subplot(3,1,1)
plot(QcQd, Droplet_D_um, 'o-b', 'LineWidth', 1.5)
xlabel('Qc/Qd'); ylabel('Droplet Diameter [µm]');
title('Droplet Diameter vs Flow Rate Ratio'); grid on;

subplot(3,1,2)
plot(QcQd, Droplet_S_um, 's-r', 'LineWidth', 1.5)
xlabel('Qc/Qd'); ylabel('Droplet Spacing [µm]');
title('Droplet Spacing vs Flow Rate Ratio'); grid on;

subplot(3,1,3)
plot(QcQd, Droplet_f, 'd-g', 'LineWidth', 1.5)
xlabel('Qc/Qd'); ylabel('Formation Frequency [Hz]');
title('Droplet Formation Frequency'); grid on;
