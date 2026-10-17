clear;
close all;
clc;

% Agrega la carpeta donde esta el solver FEM 1D.
scriptFolder = fileparts(mfilename('fullpath'));
addpath(fullfile(scriptFolder,'..','fem1D'));

% Dominio de la barra y cantidad de nodos.
L1 = 0;
L2 = 1;
N = 100;

xnode = linspace(L1,L2,N)';

% Constantes de la ecuacion:
% rho*cp*dT/dt = k*d2T/dx2 - c*T + G.
model.k = 2;
model.c = 0;
model.rho = 1;
model.cp = 1;
model.G = 100;

% Cada fila corresponde a un extremo: [tipo valor dato3].
% Tipo 1: Dirichlet [1 temperatura -1].
% Tipo 2: Neumann   [2 flujo -1].
% Tipo 3: Robin     [3 h temperatura_externa].
cb = [1 10 -1;
      1 50 -1];

% Para estudiar las matrices sin llenar la consola, usamos una malla chica.
xnode_componentes = linspace(L1,L2,6)';
componentes = mostrar_componentes(xnode_componentes,model,cb);

% et = 0: explicito, et = 1: implicito, otro valor: estacionario.
et = 2;

% Resuelve el problema completo.
T = elemFinitos(xnode,model,cb,et);

% Solucion analitica opcional para comparar.
% Primero se define como funcion y luego se evalua en los nodos.
T_analitica = @(x) -25*x.^2 + 65*x + 10;
T_analitica_nodos = T_analitica(xnode);
errorMaximo = max(abs(T(:,end) - T_analitica_nodos));

fprintf('Error maximo: %e\n',errorMaximo);
fprintf('FEM:       T(L1) = %.6f, T(L2) = %.6f\n',T(1),T(end));
fprintf('Analitica: T(L1) = %.6f, T(L2) = %.6f\n', ...
      T_analitica_nodos(1),T_analitica_nodos(end));

plot(xnode,T(:,end),'o-', ...
     xnode,T_analitica_nodos,'--','LineWidth',1.5);

legend('FEM','Solucion analitica');
xlabel('x');
ylabel('Temperatura');
grid on;