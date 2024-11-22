%%%%%%%%%%%%%%%%    Kundu
    %%%%%%%%%%%%     En esta primera version se implementa el metodo
    %%%%%%%%%%%%     MacCormak utilizando SOR. se representa la direccion y
    %%%%%%%%%%%%     modulo de la velocidad para cada timestep

% Driven Cavity by the MAC Method
clear; clc;

% Dimensiones físicas
Nx = 32; 
Ny = 32; 
Lx = 1; 
Ly = 1; 
dx = Lx / Nx; 
dy = Ly / Ny; 

Visc = 0.1; 
rho = 1.0;

% Parámetros para SOR
MaxIt = 100; 
Beta = 1.5; 
MaxErr = 0.001;

% Velocidades de la cavidad
un = 10; 
us = 0; 
ve = 0; 
vw = 0;

% Tiempo
MaxStep = 500; 
time = 0.0; 
dt = 0.002;

% Inicialización de arrays
u = zeros(Nx+1, Ny+2); v = zeros(Nx+2, Ny+1);
p = zeros(Nx+2, Ny+2); ut = zeros(Nx+1, Ny+2);
vt = zeros(Nx+2, Ny+1); pold = zeros(Nx+2, Ny+2);
c = ones(Nx+1, Ny+2) / (2/dx^2 + 2/dy^2);

% Ajuste de coeficientes en las fronteras para presión
c(2, 3:Ny) = 1 / (1/dx^2 + 2/dy^2);
c(Nx+1, 3:Ny) = 1 / (1/dx^2 + 2/dy^2);
c(3:Nx, 2) = 1 / (1/dx^2 + 2/dy^2);
c(3:Nx, Ny+1) = 1 / (1/dx^2 + 2/dy^2);
c(2, 2) = 1 / (1/dx^2 + 1/dy^2); 
c(2, Ny+1) = 1 / (1/dx^2 + 1/dy^2);
c(Nx+1, 2) = 1 / (1/dx^2 + 1/dy^2); 
c(Nx+1, Ny+1) = 1 / (1/dx^2 + 1/dy^2);

% Puntos de la malla
[x, y] = meshgrid(linspace(0, Lx, Nx+1), linspace(0, Ly, Ny+1));

% Bucle de tiempo
for is = 1:MaxStep
    % Condiciones de frontera
    u(:, 1) = 2 * us - u(:, 2);
    u(:, Ny+2) = 2 * un - u(:, Ny+1);
    v(1, :) = 2 * vw - v(2, :);
    v(Nx+2, :) = 2 * ve - v(Nx+1, :);

    % Calcular velocidades temporales
    for i = 2:Nx
        for j = 2:Ny+1
            ut(i, j) = u(i, j) + dt * ( ...
                0.25 * ( ...
                    ((u(i+1, j) + u(i, j))^2 - (u(i, j) + u(i-1, j))^2) / dx + ...
                    ((u(i, j+1) + u(i, j)) * (v(i+1, j) + v(i, j)) - ...
                     (u(i, j) + u(i, j-1)) * (v(i+1, j-1) + v(i, j-1))) / dy ...
                ) + ...
                Visc * ((u(i+1, j) - 2*u(i, j) + u(i-1, j)) / dx^2 + ...
                        (u(i, j+1) - 2*u(i, j) + u(i, j-1)) / dy^2) ...
            );
        end
    end

    for i = 2:Nx+1
        for j = 2:Ny
            vt(i, j) = v(i, j) + dt * ( ...
                0.25 * ( ...
                    ((u(i, j+1) + u(i, j)) * (v(i+1, j) + v(i, j)) - ...
                     (u(i-1, j+1) + u(i-1, j)) * (v(i, j) + v(i-1, j))) / dx + ...
                    ((v(i, j+1) + v(i, j))^2 - (v(i, j) + v(i, j-1))^2) / dy ...
                ) + ...
                Visc * ((v(i+1, j) - 2*v(i, j) + v(i-1, j)) / dx^2 + ...
                        (v(i, j+1) - 2*v(i, j) + v(i, j-1)) / dy^2) ...
            );
        end
    end

    % Resolver para presión usando SOR
    for it = 1:MaxIt
        pold = p;
        for i = 2:Nx+1
            for j = 2:Ny+1
                p(i, j) = Beta * c(i, j) * ( ...
                    ((p(i+1, j) + p(i-1, j)) / dx^2 + ...
                     (p(i, j+1) + p(i, j-1)) / dy^2) - ...
                    (rho / dt) * ((ut(i, j) - ut(i-1, j)) / dx + ...
                                  (vt(i, j) - vt(i, j-1)) / dy) ...
                ) + (1 - Beta) * p(i, j);
            end
        end

        % Comprobar convergencia
        Err = sum(sum(abs(p - pold)));
        if Err <= MaxErr, break; end
    end

    % Corregir velocidades
    u(2:Nx, 2:Ny+1) = ut(2:Nx, 2:Ny+1) - dt/dx * (p(3:Nx+1, 2:Ny+1) - p(2:Nx, 2:Ny+1));
    v(2:Nx+1, 2:Ny) = vt(2:Nx+1, 2:Ny) - dt/dy * (p(2:Nx+1, 3:Ny+1) - p(2:Nx+1, 2:Ny));


    time = time + dt;

    % módulo de velocidad total
    uu = 0.5 * (u(1:Nx+1, 2:Ny+2) + u(1:Nx+1, 1:Ny+1));
    vv = 0.5 * (v(2:Nx+2, 1:Ny+1) + v(1:Nx+1, 1:Ny+1));
    vel_mod = sqrt(uu.^2 + vv.^2);

    % Graficar resultados
    subplot(1, 2, 1);
    quiver(x, y, uu', vv', 'linewidth', 1);
    title('Campo de velocidades');
    axis equal; axis([0, Lx, 0, Ly]);
    subplot(1, 2, 2);
    contourf(x', y', vel_mod, 20, 'LineColor', 'none');
    title('Módulo de la velocidad');
    colorbar;
    axis equal; axis([0, Lx, 0, Ly]);

    pause(0.01);
end
