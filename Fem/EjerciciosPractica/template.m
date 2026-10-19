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
N = 10;                 % CAMBIAR: cantidad de nodos
xnode = linspace(L1,L2,N)';
mesh_h = (L2 - L1)/(N - 1); % h = (L2 - L1)/(N - 1)

% 2) Constantes de la ecuacion:
%    rho*cp*dT/dt = k*d2T/dx2 - c*T + G.
model.k = 1;            % CAMBIAR: conductividad k
model.c = 0;            % CAMBIAR: coeficiente de -c*T
model.rho = 1;          % CAMBIAR: densidad rho
model.cp = 1;           % CAMBIAR: calor especifico cp
model.G = 0;            % CAMBIAR: fuente constante G
% Para una fuente variable usar, por ejemplo:
% model.G = @(x) 100*(x-3).^2;
% Cuadratura de Gauss-Legendre para integrar G en cada elemento:
% 1 punto: exacta hasta grado 1; 2: hasta grado 3; 3: hasta grado 5;
% 4: hasta grado 7. Para la mayoria de las fuentes, 2 o 3 puntos bastan.
model.quadraturePoints = 2; % CAMBIAR: cantidad de puntos (1 a 4)

% 3) Condiciones de borde.
%    Fila 1: extremo izquierdo x=L1.
%    Fila 2: extremo derecho x=L2.
%
%    Dirichlet: [1 temperatura -1]
%    Neumann:   [2 flujo -1]
%    Robin:     [3 h temperatura_externa]
% En FEM, h es el coeficiente convectivo de la condicion
% -k*dT/dn = h*(T - Tinf), no h/k.
cb = [1 0 -1;             % CAMBIAR: borde izquierdo
      1 0 -1];            % CAMBIAR: borde derecho

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
et = 2;                    % CAMBIAR

% 5) Datos temporales. Se usan si et=0, et=1 o et=3.
if et == 0 || et == 1 || et == 3
    % Para Euler explicito, elegir dt pequeno por estabilidad.
    % Para Euler implicito, dt puede ser mas grande, aunque afecta
    % la precision temporal.
    usarLambda = false;    % true si el enunciado da lambda
    lambda = 0.1;          % lambda = k*dt/(rho*cp*h^2)
    model.dt = 0.001;      % CAMBIAR: dt si usarLambda = false
    if usarLambda
        model.dt = lambda*model.rho*model.cp*mesh_h^2/model.k;
    endif
    model.maxit = 1000;   % CAMBIAR: maximo de iteraciones
    model.tol = 1e-8;     % CAMBIAR: tolerancia de convergencia

    % Temperatura inicial: CAMBIAR ESTE BLOQUE segun el enunciado.
    % Opcion 1: temperatura inicial uniforme.
    temperaturaInicial = 0;
    model.T0 = temperaturaInicial*ones(size(xnode));

    % Opcion 2: temperatura inicial variable T0(x).
    % T0_funcion = @(x) 20 + 10*x;
    % model.T0 = T0_funcion(xnode);
end

% 6) Mostrar componentes para estudiar matrices y vectores.
%    Se usa una malla chica para que la consola sea legible.
xnode_componentes = linspace(L1,L2,5)';
componentes = mostrar_componentes(xnode_componentes,model,cb);

% 7) Resolver el problema.
T = elemFinitos(xnode,model,cb,et);
Tfinal = T(:,end);

if et == 0 || et == 1 || et == 3
    numberOfSteps = columns(T) - 1;
    finalDifference = T(:,end) - T(:,end-1);
    finalDenominator = norm(T(:,end),2);
    if finalDenominator == 0
        finalChange = norm(finalDifference,2);
    else
        finalChange = norm(finalDifference,2)/finalDenominator;
    endif
    finalTime = numberOfSteps*model.dt;
    fprintf('Pasos calculados: %d\n',numberOfSteps);
    fprintf('Tiempo alcanzado: %g\n',finalTime);
    fprintf('Cambio relativo final (norma 2): %e\n',finalChange);
    fprintf('Tolerancia utilizada: %e\n',model.tol);
    if finalChange < model.tol
        fprintf('Estado estacionario alcanzado por tolerancia.\n');
    else
        fprintf('No se alcanzo la tolerancia antes de maxit.\n');
    endif
endif

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

% 10) Solucion exacta y error: reemplazar [] por la funcion del enunciado.
% Ejemplo: T_analitica = @(x) 25*x.*(3-2*x);
T_analitica = [];
if !isempty(T_analitica)
    T_analitica_nodos = T_analitica(xnode);
    errorAbs = abs(Tfinal - T_analitica_nodos);
    errorRelNodal = errorAbs ./ max(abs(T_analitica_nodos),eps);
    errorRelInf = norm(Tfinal - T_analitica_nodos,inf) / ...
                  norm(T_analitica_nodos,inf);

    disp('   x          T_num       T_exacta       error_abs       error_rel');
    disp([xnode,Tfinal,T_analitica_nodos,errorAbs,errorRelNodal]);
    fprintf('Error relativo global (norma infinito) = %e\n',errorRelInf);

    figure;
    plot(xnode,Tfinal,'o-',xnode,T_analitica_nodos,'--');
    legend('FEM','Solucion exacta');
    grid on;
endif

% 11) Temperaturas en puntos pedidos por el enunciado.
% Para puntos que no coinciden con un nodo, interp1 realiza interpolacion
% lineal dentro de la malla FEM.
% puntosPedidos = [0.25, 0.5, 0.75];
% temperaturasPedidos = interp1(xnode,Tfinal,puntosPedidos);
% disp([puntosPedidos(:), temperaturasPedidos(:)]);

% 12) Flujo de calor en un punto interno.
% Para FEM lineal, la temperatura es lineal dentro de cada elemento y su
% gradiente es constante. Fourier: q(x) = -k*dT/dx.
xFlujo = 0.5;                        % CAMBIAR: punto interno pedido
if xFlujo <= L1 || xFlujo >= L2
    error('xFlujo debe estar estrictamente dentro del dominio.');
endif
elementoFlujo = find(xnode <= xFlujo,1,'last');
if elementoFlujo == length(xnode)
    elementoFlujo = elementoFlujo - 1;
endif
gradienteFlujo = (Tfinal(elementoFlujo+1) - Tfinal(elementoFlujo)) / ...
                 (xnode(elementoFlujo+1) - xnode(elementoFlujo));
flujo = -model.k*gradienteFlujo;
fprintf('Flujo en x = %g: q = %g\n',xFlujo,flujo);
