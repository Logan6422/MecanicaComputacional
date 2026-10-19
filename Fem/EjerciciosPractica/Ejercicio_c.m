clear;
close all;
clc;

% ================================================================
% TEMPLATE PARA EJERCICIOS FEM 1D
% Reemplaza solamente los datos marcados como "CAMBIAR".
% ================================================================

% Agrega las funciones FEM 1D.
scriptFolder = fileparts(mfilename('fullpath'));
addpath(fullfile(scriptFolder,'..','fem1D'));

% 1) Dominio y malla.
L1 = 1;                 % CAMBIAR: limite izquierdo
L2 = 5;                 % CAMBIAR: limite derecho
N = 10;                % CAMBIAR: cantidad de nodos
xnode = linspace(L1,L2,N)';

% 2) Constantes de la ecuacion:
%    rho*cp*dT/dt = k*d2T/dx2 - c*T + G.
model.k = 1;            % CAMBIAR: conductividad k
model.c = 0;            % CAMBIAR: coeficiente de -c*T
model.rho = 0;          % CAMBIAR: densidad rho
model.cp = 0;           % CAMBIAR: calor especifico cp
model.G = @(x) 100*(x-3).^2; % Fuente variable G(x)
% Cantidad de puntos de Gauss-Legendre para integrar la fuente G(x).
model.quadraturePoints = 2; % CAMBIAR: cantidad de puntos (1 a 4)

% 3) Condiciones de borde.
%    Fila 1: extremo izquierdo x=L1.
%    Fila 2: extremo derecho x=L2.
%
%    Dirichlet: [1 temperatura -1]
%    Neumann:   [2 flujo -1]
%    Robin:     [3 h temperatura_externa]
cb = [2 2 -1;             % CAMBIAR: borde izquierdo
      1 0 -1];            % CAMBIAR: borde derecho

% 4) Esquema temporal.
%    et=0: Euler explicito.
%    et=1: Euler implicito.
%    et=2: estacionario.
%    et=3: Crank-Nicolson.
et = 2;                    % CAMBIAR

% 5) Datos temporales. Se usan si et=0, et=1 o et=3.
if et == 0 || et == 1 || et == 3
    model.dt = 0.001;     % CAMBIAR: paso temporal
    model.maxit = 1000;   % CAMBIAR: maximo de iteraciones
    model.tol = 1e-8;     % CAMBIAR: tolerancia
    model.T0 = zeros(size(xnode));
end

% 6) Mostrar componentes para estudiar matrices y vectores.
%    Se usa una malla chica para que la consola sea legible.
xnode_componentes = linspace(L1,L2,5)';
componentes = mostrar_componentes(xnode_componentes,model,cb);

% 7) Resolver el problema.
T = elemFinitos(xnode,model,cb,et);
Tfinal = T(:,end);

% 8) Graficar el estado final.
figure;
plot(xnode,Tfinal,'o-','LineWidth',1.5);
xlabel('x');
ylabel('Temperatura');
title('Solucion FEM 1D');
grid on;

% 9) Mostrar evolucion si el problema es temporal.
if et == 0 || et == 1 || et == 3
    graficaEvolucion(xnode,T,model.dt);
end

% Si el enunciado da una solucion analitica, escribela aqui y compara:
T_analitica = @(x) (-25.*x.^4 + 300.*x.^3 -1350.*x.^2 +1906.*x + 2345)./3;
T_analitica_nodos = T_analitica(xnode);
errorMaximo = max(abs(Tfinal - T_analitica_nodos));
fprintf('Error maximo: %e\n',errorMaximo);

figure;
plot(xnode,Tfinal,'o-',xnode,T_analitica_nodos,'--','LineWidth',1.5);
legend('FEM','Solucion analitica');
xlabel('x');
ylabel('Temperatura');
title('Comparacion FEM y solucion analitica');
grid on;
