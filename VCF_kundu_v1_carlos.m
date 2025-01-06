%%%%%%%%%%%%%%%%    Kundu
    %%%%%%%%%%%%     En esta primera version se implementa el metodo
    %%%%%%%%%%%%     MacCormak utilizando SOR. se representa la direccion y
    %%%%%%%%%%%%     modulo de la velocidad para cada timestep

%% Driven Cavity by the MAC Method
clear; clc; close all;

%% Dimensiones físicas
% número de Vcontrol
Nx_values = [8,12,16,24,32]; % Diferentes tamaños de malla
% Tamaño rectángulo
Lx = 1;Ly = 1; 
Visc = 0.1; 
rho = 1.0;      
 
% Parámetros para resolución de la presión mediante SOR
MaxIt = 100; 
% beta_values = 1.2:0.1:1.5;
MaxErr = 0.001;

% Velocidades de la cavidad
un = 10;    % u tapa
us = 0;     % u pared
ve = 0; 
vw = 0;

Re= un / (Visc/rho);        % tomamos L=1m 

% Tiempo
MaxStep = 500; 
time = 0.0; 
dt = 0.002;


% Inicialización de contadores de tiempo
time_u = 0;
time_v = 0;
time_p = 0;
time_plot = 0;

% Inicializar variables para almacenar el flujo neto
time_steps = dt * (1:MaxStep); % Eje temporal
% flux_net = zeros(length(Nx_values), MaxStep);
% time_p_values = zeros(size(beta_values)); % Almacenar tiempos


%% Bucle para diferentes tamaños de malla
for n = 1:length(Nx_values)
   
    Nx=Nx_values(n); % Malla cuadrada
    Ny= Nx;
    Beta =  1.2; % Asignar Beta actual

    dx = Lx / Nx;
    dy = Ly / Ny;
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
    [Xp, Yp] = meshgrid(linspace(dx/2, Lx-dx/2, Nx), linspace(dy/2, Ly-dy/2, Ny));
    
    
   % time_p = 0; % Reiniciar tiempo de presión para esta beta


    %% Bucle tiempo
    for is = 1:MaxStep

        fprintf('Iteración %d\n', is);
    
        % Condiciones de frontera
        u(:, 1) = 2 * us - u(:, 2);
        u(:, Ny+2) = 2 * un - u(:, Ny+1);
        v(1, :) = 2 * vw - v(2, :);
        v(Nx+2, :) = 2 * ve - v(Nx+1, :);

        % Calcular velocidades temporales sin tener en cuenta la presión
        % componente horizontal
        tic;
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
        time_u = time_u + toc; % Almacena el tiempo transcurrido en calcular u

        % componente vertical
        tic;
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
        time_v = time_v + toc; % Almacena el tiempo transcurrido en calcular v

        % Resolver para presión usando SOR
        tic;
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

            % Comprobar convergencia metodo SOR
            Err = sum(sum(abs(p - pold)));
           if Err <= MaxErr, break; end
        end
        time_p = time_p + toc; % Almacena el tiempo transcurrido en calcular la presion
    

    % Corregir velocidades teniendo en cuenta la presión
    u(2:Nx, 2:Ny+1) = ut(2:Nx, 2:Ny+1) - dt/dx * (p(3:Nx+1, 2:Ny+1) - p(2:Nx, 2:Ny+1));
    v(2:Nx+1, 2:Ny) = vt(2:Nx+1, 2:Ny) - dt/dy * (p(2:Nx+1, 3:Ny+1) - p(2:Nx+1, 2:Ny));


    time = time + dt;

    % módulo de velocidad total
    uu = 0.5 * (u(1:Nx+1, 2:Ny+2) + u(1:Nx+1, 1:Ny+1));
    vv = 0.5 * (v(2:Nx+2, 1:Ny+1) + v(1:Nx+1, 1:Ny+1));
    vel_mod = sqrt(uu.^2 + vv.^2);

    % Graficar resultados
    tic;
   %  subplot(1, 3, 1);
   %  quiver(x, y, uu', vv', 'linewidth', 1);
   %  title('Campo de velocidades');
   %  axis equal; axis([0, Lx, 0, Ly]);
   %  subplot(1, 3, 2);
   %  contourf(x', y', vel_mod, 20, 'LineColor', 'none');
   % title('Módulo de la velocidad');
   %  colorbar;
   %  axis equal; axis([0, Lx, 0, Ly]);
   %  subplot(1,3,3);
   %  p_plot = p(2:Nx+1, 2:Ny+1);
   %  contourf(Xp', Yp', p_plot, 20, 'LineColor', 'none');
   %  title('Distribución de Presión');
   %  colorbar; 
   %  axis equal; axis([0, Lx, 0, Ly]);

    % pause(0.001);


    end
    % Extraer perfil de velocidades en la línea central
    center_x = round(Nx / 2); % Índice aproximado de x = 0.5
    u_centerline = uu(center_x, :)'; % Extraer las velocidades u en el centro
    y_coords = linspace(0, Ly, Ny+1); % Coordenadas y

    % Almacenar resultados para diferentes mallas
    u_profiles{n} = u_centerline;
    y_profiles{n} = y_coords;
    % 
    % Almacenar tiempo de presión para este Beta
    time_t_values(n) = time_u + time_v + time_p;
end
hold off;
figure;
b = bar(Nx_values, time_t_values); % Crear el gráfico de barras
b.FaceColor = 'flat'; % Habilitar FaceColor para colores personalizados
b.CData = jet(length(Nx_values)); % Aplicar una escala de colores (ej. 'jet')

% Personalizar etiquetas del eje X
xticks(Nx_values); % Asegura que los ticks coincidan con los valores de Nx_values
xticklabels(arrayfun(@(n) sprintf('%dx%d', n, n), Nx_values, 'UniformOutput', false));

xlabel('Mallado');
ylabel('Tiempo de cómputo total (s)');
% title('Tiempo de Cómputo para Diferentes Mallados');
grid on;


% Gráfico del campo de velocidades usando streamslice
figure;
hold on;

[Xgrid, Ygrid] = meshgrid(linspace(0, Lx, Nx+1), linspace(0, Ly, Ny+1));
% Representar las líneas de corriente usando streamslice
streamslice(Xgrid, Ygrid, uu', vv', 'linear', 'noarrows');
axis equal;
axis([0, Lx, 0, Ly]);
grid on;

hold off;


fprintf('Tiempo total en calcular u*: %.2f segundos\n', time_u);
fprintf('Tiempo total en calcular v*: %.2f segundos\n', time_v);
fprintf('Tiempo total en calcular presión (SOR): %.2f segundos\n', time_p);
fprintf('Tiempo total en graficar: %.2f segundos\n', time_plot);
fprintf('Tiempo total de la simulación: %.2f segundos\n', time_u + time_v + time_p+time_plot);
  

% Datos para el gráfico
labels = {'u*', 'v*', 'Presión (SOR)'};
times = [time_u, time_v, time_p];

% Crear el diagrama circular
figure;
pie(times, labels);

% Añadir leyenda
% legend(labels, 'Location', 'bestoutside');

% Añadir rejilla (opcional, no se aplica en gráficos pie)
% grid on;


figure;
hold on;
plot(u_profiles{1}, y_profiles{1}, 'k-', 'LineWidth', 1.5, 'DisplayName', '8x8');
if length(u_profiles) >= 2
    plot(u_profiles{2}, y_profiles{2}, 'r-', 'LineWidth', 1.5, 'DisplayName', '16x16');
end
if length(u_profiles) >= 3
    plot(u_profiles{3}, y_profiles{3}, 'b-', 'LineWidth', 1.5, 'DisplayName', '34x34');
end

xlabel('Velocidad U (m/s)');
ylabel('Coordenada Y (m)');
legend('Location', 'best');
% title('Perfil de Velocidades en la Línea Central');
grid on;
hold off;