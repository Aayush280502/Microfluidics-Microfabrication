clear; clc; close all;

% --- Channel geometry (micrometers)
L = 500;     % channel length (µm)
H = 150;     % channel height (µm)

% --- Grid
Nx = 300; Ny = 100;
x = linspace(0,L,Nx);
y = linspace(0,H,Ny);
[X,Y] = meshgrid(x,y);

% --- Flow parameters
U_center = 0.02;   % central inlet velocity (m/s)
U_side   = 0.08;   % side inlet velocity (m/s)
mu = 1e-3; rho = 1000;  % physical properties (unused here)

% --- Define velocity profile manually
u = zeros(Ny,Nx);
v = zeros(Ny,Nx);

for j = 1:Nx
    if X(1,j) < L/5
        % Left inlet region
        for i = 1:Ny
            if Y(i,j) < H/3
                u(i,j) = U_side * (1 - (Y(i,j)/(H/3) - 0.5)^2); % bottom sheath
            elseif Y(i,j) > 2*H/3
                u(i,j) = U_side * (1 - ((Y(i,j)-2*H/3)/(H/3) - 0.5)^2); % top sheath
            else
                u(i,j) = U_center * (1 - ((Y(i,j)-H/2)/(H/6))^2); % center stream
            end
        end
    else
        % After merging, smooth focusing profile
        width_focus = 0.2*H*(1 - X(1,j)/L) + 0.05*H; % narrowing width
        for i = 1:Ny
            y_rel = (Y(i,j)-H/2)/width_focus;
            u(i,j) = U_side * exp(-y_rel.^2);  % Gaussian-shaped velocity
        end
    end
end

% Smooth flow downstream
u = smoothdata(u,2,'gaussian',15);
v(:,2:end) = (u(:,2:end)-u(:,1:end-1))/10;  % weak downward motion

% --- Velocity magnitude
Umag = sqrt(u.^2 + v.^2);

% -------------------- PLOTS --------------------
figure('Color','w','Position',[100 100 1000 600]);

% Velocity magnitude
subplot(2,2,1);
contourf(X,Y,Umag,20,'LineStyle','none');
colorbar;
xlabel('x (µm)'); ylabel('y (µm)');
title('Velocity Magnitude (m/s)');
set(gca,'YDir','normal');

% Quiver (velocity vectors)
subplot(2,2,2);
stepx = 5; stepy = 3;
quiver(X(1:stepy:end,1:stepx:end), Y(1:stepy:end,1:stepx:end), ...
       u(1:stepy:end,1:stepx:end), v(1:stepy:end,1:stepx:end), 2);
xlabel('x (µm)'); ylabel('y (µm)');
title('Velocity Field (Quiver)');
set(gca,'YDir','normal');

% Streamlines
subplot(2,2,[3 4]);
contourf(X,Y,Umag,20,'LineStyle','none'); hold on;
starty = linspace(0,H,12);
streamline(X,Y,u,v,zeros(size(starty)),starty);
xlabel('x (µm)'); ylabel('y (µm)');
title('Streamlines - Flow Focusing Visualization');
set(gca,'YDir','normal');
colorbar;

sgtitle('Flow Focusing Microfluidic Device - Simplified Model');