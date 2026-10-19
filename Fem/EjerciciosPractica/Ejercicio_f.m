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
L1 = 0;                 % CAMBIAR: limite izquierdo
L2 = 1;                 % CAMBIAR: limite derecho
N = 10;                % CAMBIAR: cantidad de nodos
xnode = linspace(L1,L2,N)';

% 2) Constantes de la ecuacion:
%    rho*cp*dT/dt = k*d2T/dx2 - c*T + G.
model.k = 2;            % CAMBIAR: conductividad k
model.c = 2;            % CAMBIAR: coeficiente de -c*T
model.rho = 2;          % CAMBIAR: densidad rho
model.cp = 1;           % CAMBIAR: calor especifico cp
model.G = 75;            % CAMBIAR: fuente constante G
% Cantidad de puntos de Gauss-Legendre para integrar la fuente G.
model.quadraturePoints = 2; % CAMBIAR: cantidad de puntos (1 a 4)
% Para una fuente variable usar, por ejemplo:
% model.G = @(x) 100*(x-3).^2;

% 3) Condiciones de borde.
%    Fila 1: extremo izquierdo x=L1.
%    Fila 2: extremo derecho x=L2.
%
%    Dirichlet: [1 temperatura -1]
%    Neumann:   [2 flujo -1]
%    Robin:     [3 h temperatura_externa]
cb = [1 0 -1;             % CAMBIAR: borde izquierdo
      3 2 10];            % CAMBIAR: borde derecho

% 4) Elegi el tipo de problema.
%    et=0: Euler explicito.
%    et=1: Euler implicito.
%    et=2: estacionario.
%    et=3: Crank-Nicolson.
%
%    Alternativa A: problema estacionario.
%    et = 2;
%
%    Alternativa B: problema transitorio con Euler explicito.
%    et = 0;
%
%    Alternativa C: problema transitorio con Euler implicito.
%    et = 1;
%
%    Alternativa D: problema transitorio con Crank-Nicolson.
%    et = 3;
et = 0;                    % CAMBIAR

% 5) Datos temporales. Se usan si et=0, et=1 o et=3.
if et == 0 || et == 1 || et == 3
    % Para Euler explicito, elegir dt pequeno por estabilidad.
    % Para Euler implicito, dt puede ser mas grande, aunque afecta
    % la precision temporal.
    model.dt = 0.001;     % CAMBIAR: paso temporal
    model.maxit = 1000;   % CAMBIAR: maximo de iteraciones
    model.tol = 1e-8;     % CAMBIAR: tolerancia de convergencia

    % Opcion 1: la barra comienza completa a temperatura cero.
    model.T0 = zeros(size(xnode));

    % Opcion 2: si el enunciado da una temperatura inicial T0(x),
    % reemplazar la linea anterior por una de estas formas:
    % model.T0 = 20*ones(size(xnode));
    % model.T0 = 20 + 10*xnode;
    % model.T0 = @(x) 20 + 10*x;  % requiere evaluarla en los nodos
    % model.T0 = model.T0(xnode);
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
if et == 0 || et == 1
    graficaEvolucion(xnode,T,model.dt);
end

% 10) Comparacion opcional con solucion analitica.
% Si el enunciado proporciona T(x), escribirla como funcion anonima:
% T_analitica = @(x) ...; % CAMBIAR: formula analitica
% T_analitica_nodos = T_analitica(xnode);
% errorMaximo = max(abs(Tfinal - T_analitica_nodos));
% fprintf('Error maximo = %e\n',errorMaximo);
% figure;
% plot(xnode,Tfinal,'o-',xnode,T_analitica_nodos,'--');
% legend('FEM','Solucion analitica');
% grid on;
