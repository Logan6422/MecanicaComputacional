clear; clc; close all;

% ========================================================================
% TEMPLATE - DIFERENCIAS FINITAS 1D
% ========================================================================
% Ecuacion general:
%   rhoCp*Tt = k*Txx - c*T + G(x,t)
%
% model.transitorio = 0 -> estacionario
% model.transitorio = 1 -> transitorio
% model.tita = 0        -> Euler explicito
% model.tita = 1        -> Euler implicito
% model.tita = 1/2      -> Crank-Nicolson
%
% Condiciones de borde: [tipo, dato2, dato3]
%   tipo 1: Dirichlet -> [1, T, -1]
%   tipo 2: Neumann   -> [2, q, -1]
%   tipo 3: Robin     -> [3, hconv, Tinf]
%
% ========================================================================
% 1. PARAMETROS DEL EJERCICIO: REEMPLAZAR ESTA PARTE
% ========================================================================

model.k = 1;
model.rhoCp = 1;
model.c = 0;
model.G = 100;                       % escalar o @(x) fuente variable

mesh.xi = 0;
mesh.xf = 1;
Nnodos = 41;
mesh.h = (mesh.xf - mesh.xi)/(Nnodos - 1); % h = (xf - xi)/(N - 1)

model.transitorio = 1;

% La funcion existente interpreta Robin como:
%   -k*T_x = hconv*(T - Tinf)
% Para T_x + T = 0, con k=1, usar hconv=1 y Tinf=0.
cb = [1, 0, -1;
      3, 1,  0];

% Parametros usados solamente en el problema transitorio
model.tita = 1;
lambda = 40;                           % lambda = k*dt/(rhoCp*h^2)
model.dt = lambda * model.rhoCp * mesh.h^2 / model.k;
model.tol = 1e-4;                      % detiene si max(abs(Tnew - Told)) < tol
model.tf = 10;                         % tiempo maximo de seguridad

% ========================================================================
% 2. CONDICION INICIAL: REEMPLAZAR SI ES TRANSITORIO
% ========================================================================
xnode = mesh.xi:mesh.h:mesh.xf;
T_inicial = @(x) 10*x.*(1-x).^2;
model.T0 = T_inicial(xnode(:));

% ========================================================================
% 3. SOLUCION EXACTA: OPCIONAL, REEMPLAZAR O DEJAR []
% ========================================================================
% Debe ser una funcion de x. Ejemplo:
% exacta = @(x) 25*x.*(3-2*x);
exacta = @(x) 25*x.*(3-2*x);

% ========================================================================
% 4. RESOLUCION
% ========================================================================
T = difFinitas(mesh, model, cb);

if model.transitorio == 0
    Tfinal = T;
else
    Tfinal = T(:,end);
    cambio_final = max(abs(T(:,end) - T(:,end-1)));
    fprintf('Pasos calculados: %d\n', columns(T)-1);
    fprintf('Tiempo alcanzado: %g\n', (columns(T)-1)*model.dt);
    fprintf('Cambio final: %e (tolerancia: %e)\n', cambio_final, model.tol);

    % Primer instante en que T(x_objetivo,t) supera la temperatura buscada.
    x_objetivo = 0.5;
    T_objetivo = 23;
    [~, indice_x] = min(abs(xnode - x_objetivo));
    indice_t = find(T(indice_x,:) > T_objetivo, 1, 'first');

    if isempty(indice_t)
        fprintf('No se alcanzo T(%g,t) > %g.\n', x_objetivo, T_objetivo);
    else
        tiempo_objetivo = (indice_t - 1) * model.dt;
        fprintf(['Primer instante con T(%g,t) > %g: t = %g, ', ...
                 'T = %g.\n'], xnode(indice_x), T_objetivo, ...
                 tiempo_objetivo, T(indice_x,indice_t));
    endif
endif

% ========================================================================
% 5. COMPARACION CON LA SOLUCION EXACTA
% ========================================================================
if !isempty(exacta)
    Texacta = exacta(xnode(:));
    error_abs = abs(Tfinal - Texacta);
    error_rel = error_abs ./ max(abs(Texacta), eps);

    % Errores globales de la solucion completa.
    error_rel_inf = norm(Tfinal - Texacta, inf) / norm(Texacta, inf);
    error_rel_2 = norm(Tfinal - Texacta, 2) / norm(Texacta, 2);

    resultados = [xnode(:), Tfinal, Texacta, error_abs, error_rel];
    disp('   x          T_num       T_exacta       error_abs       error_rel');
    disp(resultados);
    fprintf('Error relativo global (norma infinito): %e\n', error_rel_inf);
    fprintf('Error relativo global (norma 2):        %e\n', error_rel_2);
endif

% ========================================================================
% 6. GRAFICO
% ========================================================================
figure;
plot(xnode, Tfinal, 'o-', 'linewidth', 1.2);
hold on;
if !isempty(exacta)
    xgraf = linspace(mesh.xi, mesh.xf, 300);
    plot(xgraf, exacta(xgraf), '--', 'linewidth', 1.2);
    legend('Diferencias finitas', 'Solucion exacta', 'location', 'best');
else
    legend('Diferencias finitas', 'location', 'best');
endif
grid on;
xlabel('x');
ylabel('T');
title('Solucion 1D');
