clear; clc; close all;

% ========================================================================
% TEMPLATE - DIFERENCIAS FINITAS 2D
% ========================================================================
% Ecuacion general estacionaria:
%   -k*(Txx + Tyy) + c*T = G
%
% model.ts = 2 -> estacionario
% model.ts = 1 -> Euler implicito
% model.ts = 0 -> Euler explicito
%
% Matrices de condiciones:
%   DIR = [nodo, valor]
%   NEU = [nodo, q, direccion]
%   ROB = [nodo, hconv, Tinf, direccion]
%
% Direcciones: 1 Sur, 2 Este, 3 Norte, 4 Oeste
%
% ========================================================================
% 1. PARAMETROS DEL EJERCICIO: REEMPLAZAR ESTA PARTE
% ========================================================================
Lx = 1;
Ly = 1;
nx = 11;                              % nodos en x
ny = 11;                              % nodos en y

model.k = 1;
model.c = 0;
model.G = 0;
model.ts = 2;

% Parametros transitorios; no se usan si model.ts = 2
model.rho = 1;
model.cp = 1;
lambda = 0.1;                         % lambda = k*dt/(rho*cp*h^2)
model.dt = 0.001;                     % se reemplaza abajo si ts = 1
model.maxit = 1000;
model.tol = 1e-6;                     % cambio relativo menor que tol

% ========================================================================
% 2. MALLA RECTANGULAR Y CONECTIVIDAD
% ========================================================================
dx = Lx/(nx-1);                       % dx = (Lx - 0)/(nx - 1)
dy = Ly/(ny-1);                       % dy = (Ly - 0)/(ny - 1)
% Para una malla cuadrada uniforme: h = dx = dy = L/(N - 1).
[X, Y] = meshgrid(0:dx:Lx, 0:dy:Ly);
xnode = [X(:), Y(:)];
model.nnodes = size(xnode,1);

icone = zeros((nx-1)*(ny-1),4);
elemento = 0;
for j = 1:ny-1
    for i = 1:nx-1
        elemento = elemento + 1;
        n1 = (j-1)*nx + i;
        icone(elemento,:) = [n1, n1+1, n1+nx+1, n1+nx];
    endfor
endfor

% k, c y G pueden ser constantes o vectores de nnodes valores.
if isscalar(model.k)
    model.k = model.k * ones(model.nnodes,1);
endif
if isscalar(model.c)
    model.c = model.c * ones(model.nnodes,1);
endif
if isscalar(model.G)
    model.G = model.G * ones(model.nnodes,1);
endif

% Para Euler implicito se puede elegir dt mediante lambda.
% En una malla no cuadrada se usa el menor paso para ser conservador.
if model.ts == 1
    h_temporal = min(dx, dy);
    model.dt = lambda * model.rho * model.cp * h_temporal^2 / min(model.k);
endif

% ========================================================================
% 3. CONDICIONES DE BORDE: REEMPLAZAR ESTA PARTE
% ========================================================================
sur = 1:nx;
norte = (ny-1)*nx + (1:nx);
oeste = 1:nx:(ny-1)*nx+1;
este = nx:nx:ny*nx;
dirichlet_nodes = unique([sur, norte, oeste, este]);

% Ejemplo: Dirichlet T=0 en todo el contorno.
DIR = [dirichlet_nodes(:), zeros(length(dirichlet_nodes),1)];

% Quitar de DIR los nodos que se quieran usar como Neumann o Robin.
NEU = [];                              % [nodo, q, direccion]
ROB = [];                              % [nodo, hconv, Tinf, direccion]

% ========================================================================
% 4. CONDICION INICIAL: REEMPLAZAR SI ES TRANSITORIO
% ========================================================================
T_inicial = @(x,y) 0*x + 0*y;
model.PHI_n = T_inicial(xnode(:,1), xnode(:,2));

% ========================================================================
% 5. SOLUCION EXACTA: OPCIONAL, REEMPLAZAR O DEJAR []
% ========================================================================
% Debe ser una funcion de dos variables. Ejemplo:
% exacta = @(x,y) sin(pi*x).*sin(pi*y);
exacta = [];

% ========================================================================
% 6. RESOLUCION
% ========================================================================
[PHI, Q] = fdm2d(xnode, icone, DIR, NEU, ROB, model);
PHIfinal = PHI(:,end);

if model.ts ~= 2
    cambio_final = max(abs(PHI(:,end) - PHI(:,end-1)));
    fprintf('Iteraciones calculadas: %d\n', columns(PHI)-1);
    fprintf('Tiempo alcanzado: %g\n', (columns(PHI)-1)*model.dt);
    fprintf('Cambio final: %e (tolerancia: %e)\n', cambio_final, model.tol);
endif

% ========================================================================
% 7. COMPARACION CON LA SOLUCION EXACTA
% ========================================================================
if !isempty(exacta)
    Texacta = exacta(xnode(:,1), xnode(:,2));
    Texacta = Texacta(:);
    error_abs = abs(PHIfinal - Texacta);
    error_rel = error_abs ./ max(abs(Texacta), eps);

    % Errores globales de la solucion completa.
    error_rel_inf = norm(PHIfinal - Texacta, inf) / norm(Texacta, inf);
    error_rel_2 = norm(PHIfinal - Texacta, 2) / norm(Texacta, 2);

    resultados = [xnode, PHIfinal, Texacta, error_abs, error_rel];
    disp('   x          y          T_num       T_exacta       error_abs       error_rel');
    disp(resultados);
    fprintf('Error relativo global (norma infinito): %e\n', error_rel_inf);
    fprintf('Error relativo global (norma 2):        %e\n', error_rel_2);
endif

% ========================================================================
% 8. GRAFICO
% ========================================================================
figure;
trisurf([icone(:,1:3); icone(:,[1,3,4])], ...
        xnode(:,1), xnode(:,2), PHIfinal);
shading interp;
colorbar;
view(2);
axis equal tight;
xlabel('x');
ylabel('y');
title('Solucion 2D');
