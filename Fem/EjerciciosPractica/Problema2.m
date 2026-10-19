clear;
close all;
clc;

% ========================================================================
% PROBLEMA 2 - FEM 2D
% Placa 1x1, k=10, G=500, rho*cp=1.
% Euler explicito, T0=125 y tolerancia 1e-5.
% ========================================================================

scriptFolder = fileparts(mfilename('fullpath'));
addpath(fullfile(scriptFolder,'..','fem2D'));

% 1) Dominio: dx = dy = 1/6 -> 7 nodos por direccion.
Lx = 1;
Ly = 1;
nx = 7;
ny = 7;
dx = Lx/(nx-1);
dy = Ly/(ny-1);
[X,Y] = meshgrid(0:dx:Lx,0:dy:Ly);
xnode = [X(:),Y(:)];

model.nnodes = size(xnode,1);
model.nelem = (nx-1)*(ny-1);

% Conectividad antihoraria: izquierda-abajo, derecha-abajo,
% derecha-arriba, izquierda-arriba.
icone = zeros(model.nelem,4);
element = 0;
for j = 1:ny-1
    for i = 1:nx-1
        element = element + 1;
        n1 = (j-1)*nx+i;
        icone(element,:) = [n1,n1+1,n1+nx+1,n1+nx];
    endfor
endfor

% 2) Propiedades fisicas.
model.kx = 10;
model.ky = 10;
model.c = 0;
model.rho = 1;
model.cp = 1;
model.G = 500*ones(model.nelem,1);
model.ts = 0;                         % Euler explicito
model.maxit = 100000;
model.tol = 1e-5;
model.PHI_n = 125*ones(model.nnodes,1);

% 3) Condiciones de borde.
% Superior: T=150. Se prioriza esta temperatura en la esquina superior
% izquierda, donde tambien pertenece al borde izquierdo.
top = (ny-1)*nx+(1:nx);
leftInterior = 1+nx: nx:(ny-2)*nx+1;
DIR = [top(:),150*ones(length(top),1);
       leftInterior(:),100*ones(length(leftInterior),1)];

% Inferior: -k*dT/dn = h*(T-Tinf), h=50, Tinf=20.
bottomEdges = zeros(nx-1,4);
for i = 1:nx-1
    bottomEdges(i,:) = [i,i+1,50,20];
endfor
ROB = bottomEdges;

% Derecho aislado: q=0. Cada fila representa un lado de elemento.
rightEdges = zeros(ny-1,3);
for j = 1:ny-1
    n1 = j*nx;
    rightEdges(j,:) = [n1,n1+nx,0];
endfor
NEU = rightEdges;
PUN = [];

% 4) Paso explicito estable calculado por la implementacion.
dt = fem2d_heat_explicit_delta_t(xnode,icone,model);
[PHI,Q] = fem2d_heat(xnode,icone,DIR,NEU,ROB,PUN,model);
PHIfinal = PHI(:,end);
Qfinal = Q(:,end-1:end);

% 5) Tiempo de convergencia.
numberOfSteps = columns(PHI)-1;
finalTime = numberOfSteps*dt;
finalDifference = PHI(:,end)-PHI(:,end-1);
finalChange = norm(finalDifference,2)/max(norm(PHI(:,end),2),eps);
fprintf('Nodos: %d | Elementos: %d\n',model.nnodes,model.nelem);
fprintf('dx = dy = %g\n',dx);
fprintf('dt = %e\n',dt);
fprintf('Pasos calculados: %d\n',numberOfSteps);
fprintf('Tiempo de estado estacionario: %g\n',finalTime);
fprintf('Cambio relativo final: %e\n',finalChange);

% 6) Temperaturas pedidas. (0.25,0.25) se interpola dentro de un elemento.
puntos = [0.25,0.25;1.0,0.0;0.5,0.5];
for p = 1:size(puntos,1)
    xp = puntos(p,1);
    yp = puntos(p,2);
    elementoPunto = fem2d_findElement(xnode,icone,xp,yp);
    nodosPunto = icone(elementoPunto,:);
    N = fem2d_heat_blerp(xnode(nodosPunto,:),xp,yp);
    temperatura = N'*PHIfinal(nodosPunto);
    fprintf('T(%g,%g) = %g\n',xp,yp,temperatura);
endfor

% 7) Flujo de calor en P=(0.5,0.5).
% La funcion entrega Q=[qx,qy] promedio de los elementos que comparten el nodo.
centro = find(abs(xnode(:,1)-0.5)<eps & abs(xnode(:,2)-0.5)<eps,1);
flujoCentro = Qfinal(centro,:);
fprintf('Flujo en P(0.5,0.5): qx = %g, qy = %g, |q| = %g\n', ...
        flujoCentro(1),flujoCentro(2),norm(flujoCentro,2));

% 8) Grafico de temperatura final.
figure;
trisurf([icone(:,1:3);icone(:,[1,3,4])],xnode(:,1),xnode(:,2),PHIfinal);
shading interp;
colorbar;
view(2);
xlim([0,Lx]);
ylim([0,Ly]);
title('Problema 2 - Temperatura estacionaria FEM 2D');
xlabel('x');
ylabel('y');

