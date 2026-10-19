clear;
close all;
clc;

% Agrega la carpeta donde esta el solver FEM 1D.
scriptFolder = fileparts(mfilename('fullpath'));
addpath(fullfile(scriptFolder,'..','fem1D'));

% 1) Copia aqui los limites del dominio del enunciado.
L1 = 0;
L2 = 2;

% 2) Elegi la cantidad de nodos. Mas nodos suele dar mayor precision.
N = 20;
xnode = linspace(L1,L2,N)';

% 3) Carga las constantes de la ecuacion:
%    rho*cp*dT/dt = k*d2T/dx2 - c*T + G.
model.k = 1;
model.c = 1;
model.rho = 0;
model.cp = 0;
model.G = 0;
% Cantidad de puntos de Gauss-Legendre para integrar la fuente G.
model.quadraturePoints = 2;

% 4) Completa las condiciones de borde.
%    Fila 1: extremo izquierdo x=L1.
%    Fila 2: extremo derecho x=L2.
%
%    Dirichlet: [1 temperatura -1]
%    Neumann:   [2 flujo -1]
%    Robin:     [3 h temperatura_externa]
cb = [1 100 -1;
      2 0 -1];

% 4b) Muestra las matrices del FEM para estudiarlas.
%     Usamos pocos nodos para que la salida sea legible.
xnode_componentes = linspace(L1,L2,5)';
componentes = mostrar_componentes(xnode_componentes,model,cb);

% 5) Elegi el esquema.
%    et=0: Euler explicito.
%    et=1: Euler implicito.
%    et=2: estacionario.
%    et=3: Crank-Nicolson.
et = 2;

% 6) Si et es 0, 1 o 3, agrega estos datos temporales.
%    T0 debe tener una temperatura inicial por nodo.
if et == 0 || et == 1 || et == 3
    model.dt = 0.001;
    model.maxit = 1000;
    model.tol = 1e-8;
    model.T0 = zeros(size(xnode));
end

% 7) Resuelve todo el problema.
T = elemFinitos(xnode,model,cb,et);

% 8) La solucion final es la ultima columna.
Tfinal = T(:,end);

% 9) Grafica la solucion final.
figure;
plot(xnode,Tfinal,'o-','LineWidth',1.5);
xlabel('x');
ylabel('Temperatura');
title('Solucion FEM 1D');
grid on;

% 10) En un problema temporal, muestra toda la evolucion.
if et == 0 || et == 1 || et == 3
    graficaEvolucion(xnode,T,model.dt);
end

% Si el enunciado da una solucion analitica, escribela aqui y compara:
T_analitica = @(x) (100.*exp(-x).*(exp(2.*x) + exp(4)))./(1+exp(4));
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
