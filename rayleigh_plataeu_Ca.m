
clear; clc; close all;

% Physical properties
rho_c = 1000;    % continuous density (kg/m^3)
rho_d = 900;     % dispersed density (kg/m^3)
mu_c  = 1e-3;    % continuous viscosity (Pa.s)
mu_d  = 2e-3;    % dispersed viscosity (Pa.s)
sigma = 0.03;    % interfacial tension (N/m)

% Geometry: all meters
W_main = 120e-6;      % main channel height (y-direction)
W_orifice = 40e-6;    % orifice width (narrow constriction)
depth = 50e-6;        % channel depth (z-direction) used for volume estimate
L_display = 1.6e-3;   % downstream display length (x-direction for animation)

% Flow rates and ratios
Qd_each = 5e-10;      % flow rate of dispersed phase from each side inlet (top and bottom) [m^3/s]
% therefore total dispersed Qd_total = 2*Qd_each
Qc_over_Qd_total = 40;  % set Qc/Qd_total ratio (continuous flow relative to total dispersed)
Qd_total = 2 * Qd_each;
Qc = Qc_over_Qd_total * Qd_total;

% Derived velocities (mean using cross-sectional area = width * depth)
U_c = Qc / (W_main * depth);          % continuous mean velocity (m/s)
U_d_each = Qd_each / (W_orifice * depth); % each jet mean velocity (m/s)

% Capillary number (based on continuous phase)
Ca = mu_c * U_c / sigma;

% regime tuning constants (empirical / adjustable)
A = 1.0; B = 0.9; C = 0.35;

% which Qc/Qd points to sweep for plot (optional)
QcQd_sweep = logspace(log10(1), log10(200), 50);


Ca_squeeze_drip = 0.01;
Ca_drip_jet = 0.1;
QcQd_jet_threshold = 80;


Qratio_total = Qc / Qd_total;
qratio_inv = 1 / Qratio_total; % Qd_total / Qc

if (Ca < Ca_squeeze_drip) && (Qratio_total < QcQd_jet_threshold)
    regime = 'Squeezing';
    D = W_main * (1 + A * qratio_inv);    % geometry-dominated estimate
    k = 4.0;
elseif (Ca >= Ca_drip_jet) || (Qratio_total > QcQd_jet_threshold)
    regime = 'Jetting';
    D = W_main * C * (Ca)^(-0.10) * (qratio_inv)^0.05;
    k = 1.2;
else
    regime = 'Dripping';
    D = W_main * B * (qratio_inv)^0.30 * (1 + 0.5*Ca^0.05);
    k = 2.8;
end

% Enforce physical bounds
D = max(0.05*W_main, min(3*W_main, D));
S = k * D;

% Estimate droplet frequency from mass conservation:
% Use volume approx: droplet volume Vdrop = (pi/6)*D^3 (assume spherical)
Vdrop = (pi/6) * D^3;
f = Qd_total / Vdrop;   % Hz

% Print computed values
fprintf('\nOperating point (cross-constriction):\n');
fprintf(' Regime: %s\n', regime);
fprintf(' Qc/Qd_total = %.1f\n', Qratio_total);
fprintf(' Ca = %.3e\n', Ca);
fprintf(' Predicted droplet diameter D = %.2f µm\n', D*1e6);
fprintf(' Predicted spacing S = %.2f µm\n', S*1e6);
fprintf(' Predicted formation frequency f = %.2f Hz\n\n', f);


% Animation parameters
num_drops_display = 12;
% initial positions upstream so they appear to form as animation runs
drop_positions = - (0:num_drops_display-1) * S * 1.1;   % in meters
drop_radius = D / 2;
dt_anim = min(1/(20*max(1,f)), 1e-3);   % animation time step
frames = 300;

% Setup figure
figure('Color','w','Name','Cross-constriction droplet animation');
for frame = 1:frames
    clf; hold on; axis equal;
    % plotting window
    xlim([-0.2e-3, L_display*1e0]*1e6);   % µm scale for display
    ylim([-1.3*W_main, 1.3*W_main]*1e6);
    xlabel('x (\mum)'); ylabel('y (\mum)');
    title(sprintf('Cross-constriction animation — %s regime', regime));
    
    % Draw main channel boundary (centered at y=0)
    rectangle('Position',[-0.2e-3*1e6 -W_main*1e6 (L_display+0.2e-3)*1e6 2*W_main*1e6], 'EdgeColor','k','LineWidth',1.5);
    % Draw top and bottom side inlet rectangles (left side)
    inlet_w_display = 0.22e-3; % display upstream inlet length (m)
    rectangle('Position',[-0.2e-3*1e6 W_main*0.5*1e6 inlet_w_display*1e6 (W_main*0.5)*1e6], 'FaceColor',[0.95 0.95 0.95],'EdgeColor','none');
    rectangle('Position',[-0.2e-3*1e6 -W_main*1e6 inlet_w_display*1e6 (W_main*0.5)*1e6], 'FaceColor',[0.95 0.95 0.95],'EdgeColor','none');
    
    % Draw orifice in center (narrow constriction)
    or_x = 0; % orifice location x=0
    rectangle('Position',[or_x*1e6 -W_orifice*0.5*1e6 W_orifice*1e6 W_orifice*1e6], 'EdgeColor','r','LineWidth',1.8);
    
    % Advance droplet positions (advect with main flow U_c)
    drop_positions = drop_positions + U_c * dt_anim;
    
    % When first (leading) droplet has moved sufficiently, spawn a new droplet at x ~ 0
    if drop_positions(1) > (S*0.9)
        drop_positions = [0, drop_positions(1:end-1)];
    end
    
    % Draw droplets as filled ellipses centered at mid-channel (y=0)
    for k = 1:num_drops_display
        xk = drop_positions(k);
        if xk > -0.5e-3 && xk < L_display
            % to keep droplet visible in 2D, draw ellipse with height = D*0.9, width = D*1.1
            hx = 1.1 * drop_radius; hy = 0.9 * drop_radius;
            rectangle('Position',[(xk-hx)*1e6 -hy*1e6 2*hx*1e6 2*hy*1e6], 'Curvature',[1 1], ...
                'FaceColor',[0.2 0.6 0.9],'EdgeColor','none');
        end
    end
    
    % Display info text
    text(L_display*0.6*1e6, 1.05*W_main*1e6, sprintf('D = %.1f µm', D*1e6), 'FontSize', 11);
    text(L_display*0.6*1e6, 0.85*W_main*1e6, sprintf('S = %.1f µm', S*1e6), 'FontSize', 11);
    text(L_display*0.6*1e6, 0.65*W_main*1e6, sprintf('f = %.1f Hz', f), 'FontSize', 11);
    text(L_display*0.6*1e6, 0.45*W_main*1e6, sprintf('Qc/Qd_total = %.1f', Qratio_total), 'FontSize', 11);
    text(L_display*0.6*1e6, 0.25*W_main*1e6, sprintf('Ca = %.2e', Ca), 'FontSize', 11);
    
    drawnow;
end

QcQd_vals = logspace(log10(1), log10(200), 40);
Ca_vals = [1e-4, 5e-4, 1e-3, 5e-3, 2e-2, 1e-1]; % sample Ca values
Dmap = zeros(length(Ca_vals), length(QcQd_vals));
Smap = Dmap;
Fmap = Dmap;

for ii = 1:length(Ca_vals)
    Ca_here = Ca_vals(ii);
    for jj = 1:length(QcQd_vals)
        Qratio = QcQd_vals(jj);
        qinv = 1 / Qratio;
        % regime decision (same logic)
        if (Ca_here < Ca_squeeze_drip) && (Qratio < QcQd_jet_threshold)
            Dtmp = W_main * (1 + A * qinv);
            ktmp = 4.0;
        elseif (Ca_here >= Ca_drip_jet) || (Qratio > QcQd_jet_threshold)
            Dtmp = W_main * C * (Ca_here)^(-0.10) * (qinv)^0.05;
            ktmp = 1.2;
        else
            Dtmp = W_main * B * (qinv)^0.30 * (1 + 0.5*Ca_here^0.05);
            ktmp = 2.8;
        end
        Dtmp = max(0.05*W_main, min(3*W_main, Dtmp));
        Dmap(ii,jj) = Dtmp;
        Smap(ii,jj) = ktmp * Dtmp;
        Vtmp = (pi/6) * Dtmp^3;
        Fmap(ii,jj) = (2*Qd_each) / Vtmp; % use total dispersed flow
    end
end

% Plot sweep results
figure('Color','w','Name','Parametric sweep');
subplot(3,1,1)
semilogx(QcQd_vals, Dmap'*1e6);
xlabel('Qc/Qd_{total}'); ylabel('D (µm)');
title('Droplet diameter vs Qc/Qd_{total} for different Ca (curves = different Ca)');
legend(arrayfun(@(c) sprintf('Ca=%.0e',Ca_vals(c)),1:length(Ca_vals),'UniformOutput',false),'Location','northeast');

subplot(3,1,2)
semilogx(QcQd_vals, Smap'*1e6);
xlabel('Qc/Qd_{total}'); ylabel('Spacing S (µm)');
title('Spacing vs Qc/Qd_{total}');

subplot(3,1,3)
loglog(QcQd_vals, Fmap');
xlabel('Qc/Qd_{total}'); ylabel('Formation frequency f (Hz)');
title('Frequency vs Qc/Qd_{total}');
legend(arrayfun(@(c) sprintf('Ca=%.0e',Ca_vals(c)),1:length(Ca_vals),'UniformOutput',false),'Location','northeast');

% End of script
