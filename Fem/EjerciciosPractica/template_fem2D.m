clear;
close all;
clc;

% ========================================================================
% TEMPLATE FEM 2D - ECUACION DEL CALOR
% ========================================================================
% rho*cp*dT/dt - div(k*grad(T)) + c*T = G
%
% model.ts = 0: Euler explicito
% model.ts = 1: Euler implicito
% model.ts = 2: estacionario
%
% DIR: [nodo, temperatura]
% NEU: [nodo1, nodo2, flujo q]
% ROB: [nodo1, nodo2, hconv, Tinf]
% PUN: [x, y, fuente puntual]
%
% Las filas de ROB y NEU representan lados de elementos, no nodos aislados.
% ========================================================================

scriptFolder = fileparts(mfilename('fullpath'));
addpath(fullfile(scriptFolder,'..','fem2D'));

% 1) DOMINIO, MALLA Y NUMERACION
Lx = 1;
Ly = 1;
nx = 7;                              % cantidad de nodos en x
ny = 7;                              % cantidad de nodos en y

dx = Lx/(nx-1);                     % dx = Lx/(nx-1)
dy = Ly/(ny-1);                     % dy = Ly/(ny-1)
[X,Y] = meshgrid(0:dx:Lx,0:dy:Ly);
xnode = [X(:),Y(:)];
model.nnodes = size(xnode,1);

icone = zeros((nx-1)*(ny-1),4);
element = 0;
for j = 1:ny-1
    for i = 1:nx-1
        element = element + 1;
        n1 = (j-1)*nx + i;
        icone(element,:) = [n1,n1+1,n1+nx+1,n1+nx];
    endfor
endfor
model.nelem = size(icone,1);

% 2) PROPIEDADES FISICAS
model.kx = 1;
model.ky = 1;
model.c = 0;
model.rho = 1;
model.cp = 1;
model.G = zeros(model.nelem,1);      % escalar o un valor por elemento

% 3) CONDICIONES DE BORDE
% Dirichlet: seleccionar nodos y asignar temperatura.
DIR = [];
% Neumann: [nodo1, nodo2, q], con q = -k*dT/dn segun la convencion.
NEU = [];
% Robin: [nodo1, nodo2, hconv, Tinf].
ROB = [];
% Fuentes puntuales: [x, y, valor].
PUN = [];

% 4) ESQUEMA Y DATOS TEMPORALES
model.ts = 0;                         % CAMBIAR: 0 explicito, 1 implicito, 2 estacionario
model.maxit = 100000;
model.tol = 1e-5;
model.PHI_n = zeros(model.nnodes,1);  % CAMBIAR: T inicial

% Para ts = 1 se usa model.dt directamente.
model.dt = 0.001;

% Para ts = 0 el solver calcula automaticamente un dt estable conservador.
if model.ts == 0
    dt = fem2d_heat_explicit_delta_t(xnode,icone,model);
else
    dt = model.dt;
endif

% 5) RESOLUCION
[PHI,Q] = fem2d_heat(xnode,icone,DIR,NEU,ROB,PUN,model);
PHIfinal = PHI(:,end);

% 6) TIEMPO Y CONVERGENCIA
if model.ts ~= 2
    numberOfSteps = columns(PHI)-1;
    finalTime = numberOfSteps*dt;
    finalDifference = PHI(:,end)-PHI(:,end-1);
    finalDenominator = norm(PHI(:,end),2);
    if finalDenominator == 0
        finalChange = norm(finalDifference,2);
    else
        finalChange = norm(finalDifference,2)/finalDenominator;
    endif
    fprintf('Pasos calculados: %d\n',numberOfSteps);
    fprintf('dt utilizado: %e\n',dt);
    fprintf('Tiempo alcanzado: %g\n',finalTime);
    fprintf('Cambio relativo final: %e\n',finalChange);
    fprintf('Tolerancia: %e\n',model.tol);
endif

% 7) TEMPERATURAS EN PUNTOS PEDIDOS
% Para un punto que no es nodo, buscar el elemento y usar interpolacion.
puntos = [0.25,0.25; 1.0,0.0; 0.5,0.5];
for p = 1:size(puntos,1)
    x = puntos(p,1);
    y = puntos(p,2);
    elementoPunto = fem2d_findElement(xnode,icone,x,y);
    if elementoPunto == 0
        error('El punto pedido no pertenece a la malla.');
    endif
    nodosPunto = icone(elementoPunto,:);
    N = fem2d_heat_blerp(xnode(nodosPunto,:),x,y);
    temperatura = N'*PHIfinal(nodosPunto);
    fprintf('T(%g,%g) = %g\n',x,y,temperatura);
endfor

% 8) FLUJO EN UN PUNTO INTERNO
puntoFlujo = [0.5,0.5];
elementoFlujo = fem2d_findElement(xnode,icone,puntoFlujo(1),puntoFlujo(2));
Qfinal = Q(:,end-1:end);
nodosFlujo = icone(elementoFlujo,:);
% En un nodo del elemento, Qfinal contiene el flujo promedio de elementos vecinos.
indiceCentro = find(abs(xnode(:,1)-puntoFlujo(1)) < eps & ...
                    abs(xnode(:,2)-puntoFlujo(2)) < eps,1);
flujo = Qfinal(indiceCentro,:);
fprintf('Flujo en P(%g,%g): qx = %g, qy = %g, |q| = %g\n', ...
        puntoFlujo(1),puntoFlujo(2),flujo(1),flujo(2),norm(flujo,2));

% 9) SOLUCION EXACTA OPCIONAL
% exacta = @(x,y) ...;
exacta = [];
if !isempty(exacta)
    Texacta = exacta(xnode(:,1),xnode(:,2));
    Texacta = Texacta(:);
    errorAbs = abs(PHIfinal-Texacta);
    errorRelNodal = errorAbs./max(abs(Texacta),eps);
    errorRelInf = norm(PHIfinal-Texacta,inf)/norm(Texacta,inf);
    disp([xnode,PHIfinal,Texacta,errorAbs,errorRelNodal]);
    fprintf('Error relativo global (norma infinito): %e\n',errorRelInf);
endif

% 10) GRAFICO FINAL
figure;
trisurf([icone(:,1:3);icone(:,[1,3,4])],xnode(:,1),xnode(:,2),PHIfinal);
shading interp;
colorbar;
view(2);
xlim([0,Lx]);
ylim([0,Ly]);
title('Solucion FEM 2D');
xlabel('x');
ylabel('y');

