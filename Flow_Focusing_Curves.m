clear; clc; close all;

figure('Color','w','Position',[100 100 1200 700]);

Qc_Qd = [10 20 40 60 80 100 150 200];

D_DOX = 180*(Qc_Qd.^(-0.35));
D_CUR = 200*(Qc_Qd.^(-0.28));
D_TAM = 220*(Qc_Qd.^(-0.22));

Eff_DOX = 98*(1 - exp(-0.03*Qc_Qd));
Eff_CUR = 93*(1 - exp(-0.02*Qc_Qd));
Eff_TAM = 97*(1 - exp(-0.025*Qc_Qd));

x = linspace(0,2,100);
MI_DOX = 1 - exp(-3.0*x);
MI_CUR = 1 - exp(-2.0*x);
MI_TAM = 1 - exp(-1.2*x);

y = linspace(-150,150,200);
Conc_DOX = 0.5*(1 + erf(-y/40));  
Conc_CUR = 0.5*(1 + erf(-y/60));
Conc_TAM = 0.5*(1 + erf(-y/90)); 

colors = [0.0 0.45 0.74; 0.85 0.33 0.10; 0.47 0.67 0.19];

%% Droplet Diameter vs Flow Rate Ratio
subplot(2,2,1);
plot(Qc_Qd,D_DOX,'-o','Color',colors(1,:),'LineWidth',2); hold on;
plot(Qc_Qd,D_CUR,'-s','Color',colors(2,:),'LineWidth',2);
plot(Qc_Qd,D_TAM,'-^','Color',colors(3,:),'LineWidth',2);
xlabel('Flow Rate Ratio (Q_c / Q_d)');
ylabel('Droplet Diameter (µm)');
title('Droplet Size vs Flow Rate Ratio');
legend('DOX (PLGA-PEG)','Curcumin (PCL-PEG)','Tamoxifen (PLGA)','Location','northeast');
grid on;

%% Encapsulation Efficiency vs Flow Rate Ratio
subplot(2,2,2);
plot(Qc_Qd,Eff_DOX,'-o','Color',colors(1,:),'LineWidth',2); hold on;
plot(Qc_Qd,Eff_CUR,'-s','Color',colors(2,:),'LineWidth',2);
plot(Qc_Qd,Eff_TAM,'-^','Color',colors(3,:),'LineWidth',2);
xlabel('Flow Rate Ratio (Q_c / Q_d)');
ylabel('Encapsulation Efficiency (%)');
title('Encapsulation Efficiency vs Flow Rate Ratio');
legend('DOX','Curcumin','Tamoxifen','Location','southeast');
grid on;

%% Mixing Index along Channel
subplot(2,2,3);
plot(x,MI_DOX,'-','Color',colors(1,:),'LineWidth',2); hold on;
plot(x,MI_CUR,'--','Color',colors(2,:),'LineWidth',2);
plot(x,MI_TAM,':','Color',colors(3,:),'LineWidth',2);
xlabel('Channel Length (mm)');
ylabel('Mixing Index (0–1)');
title('Mixing Efficiency along Channel');
legend('DOX','Curcumin','Tamoxifen','Location','southeast');
ylim([0 1]); grid on;

%% Concentration Profile across Junction
subplot(2,2,4);
plot(y,Conc_DOX,'-','Color',colors(1,:),'LineWidth',2); hold on;
plot(y,Conc_CUR,'--','Color',colors(2,:),'LineWidth',2);
plot(y,Conc_TAM,':','Color',colors(3,:),'LineWidth',2);
xlabel('Position across Channel (µm)');
ylabel('Normalized Concentration');
title('Cross-Junction Mixing Profile');
legend('DOX','Curcumin','Tamoxifen','Location','best');
grid on;

sgtitle('Microfluidic Encapsulation of Real Drugs (4-Inlet Device)','FontSize',14,'FontWeight','bold');
