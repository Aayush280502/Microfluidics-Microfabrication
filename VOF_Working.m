
clear; clc; close all;

Lx = 400e-6;    % domain length (m)
Ly = 150e-6;    % domain height (m)
Nx = 200; Ny = 75;
dx = Lx/Nx; dy = Ly/Ny;
[x, y] = meshgrid(linspace(0, Lx, Nx), linspace(0, Ly, Ny));

% ------------------ Physical properties ------------------
rho1 = 1000; rho2 = 800;
mu1 = 1e-3; mu2 = 5e-3;
sigma = 0.03;        % surface tension (N/m)
U_in = 0.02;         % inlet velocity (m/s)

% ------------------ Time parameters ------------------
dt = 1e-6; Nt = 30000;  % smaller dt for stability

% ------------------ Initialization ------------------
phi = ones(Ny, Nx);             % level set (1 = continuous, -1 = dispersed)
phi(:, 1:round(Nx/6)) = -1;     % left dispersed phase
phi = smoothdata(phi, 2, 'gaussian', 5);

u = zeros(Ny, Nx); v = zeros(Ny, Nx);

% Inlet velocity profile
u(:,1) = U_in; u(:,end) = U_in;
v(:,1) = 0; v(:,end) = 0;

% ------------------ Simulation loop ------------------
figure('Color', 'w', 'Position', [100 100 1000 400]);

for n = 1:Nt
    % --- Compute gradients (finite differences)
    [phi_x, phi_y] = gradient(phi, dx, dy);
    mag_grad = sqrt(phi_x.^2 + phi_y.^2 + 1e-12);
    nx = phi_x ./ mag_grad;
    ny = phi_y ./ mag_grad;

    % --- Curvature (∇·n)
    [nx_x, ~] = gradient(nx, dx, dy);
    [~, ny_y] = gradient(ny, dx, dy);
    curvature = nx_x + ny_y;
    curvature(~isfinite(curvature)) = 0;

    % --- Surface tension force
    Fst = sigma * curvature .* nx;

    % --- Velocity update (advection only)
    u(:,2:end-1) = U_in .* exp(-((y(:,2:end-1)-Ly/2).^2)/(2*(Ly/6)^2));
    v = 0.0001 * sin(2*pi*y/Ly);  % weak perturbation

    % --- Advect phi (upwind scheme)
    phi_xadv = [diff(phi,1,2), zeros(Ny,1)] / dx;
    phi_yadv = [diff(phi,1,1); zeros(1,Nx)] / dy;
    phi_new = phi - dt*(u.*phi_xadv + v.*phi_yadv);

    % --- Clamp and smooth phi
    phi = min(max(phi_new, -1), 1);
    phi = smoothdata(phi, 2, 'gaussian', 3);

    % --- Visualization ---
    if mod(n,40)==0
        contourf(x*1e6, y*1e6, phi, 20, 'LineStyle', 'none');
        xlabel('x (µm)'); ylabel('y (µm)');
        title(['Flow Focusing VOF Approximation | t = ' num2str(n*dt*1e3, '%.2f') ' ms']);
        colorbar; clim([-1 1]);
        set(gca, 'YDir', 'normal');
        drawnow;
    end
end