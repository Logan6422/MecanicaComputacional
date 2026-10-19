clear;
clc;

% Problema 1: T_t = T_xx + 100, 0 < x < 1.
% T(0,t) = 0, T_x(1,t) = -T(1,t).

addpath('../1D');

% 1) Estado estacionario con 4 nodos.
model.rhoCp = 0;
model.k = 1;
model.c = 0;
model.G = 100;
model.transitorio = 0;

mesh.xi = 0;
mesh.xf = 1;
mesh.h = 1/3;

% Tipo 1: Dirichlet. Tipo 3: Robin [h Tinf].
cb = [1 0 -1;
			3 1  0];

T4 = difFinitas(mesh, model, cb);
x4 = mesh.xi:mesh.h:mesh.xf;
Tex = @(x) 25*x.*(3 - 2*x);
errorRel4 = abs((T4 - Tex(x4)') ./ Tex(x4)');
errorRel4(1) = 0;

fprintf('Punto 1: estacionario con 4 nodos\n');
disp([x4(:), T4, Tex(x4(:)), errorRel4]);

% 2) Transitorio implícito con 41 nodos y lambda = k*dt/h^2 = 40.
model.rhoCp = 1;
model.transitorio = 1;
model.tita = 1;              % Euler implícito
mesh.h = 1/40;
dt = 40*mesh.h^2/model.k;
model.dt = dt;
model.tf = dt;               % Se avanza un único paso por llamada.

x = mesh.xi:mesh.h:mesh.xf;
model.T0 = (10*x.*(1-x).^2)';

tol = 1e-4;
maxPasos = 100000;
temperatura = model.T0;
tiempo = 0;
primerTiempo24 = NaN;
historial = temperatura;

for paso = 1:maxPasos
	model.T0 = temperatura;
	temperaturaNueva = difFinitas(mesh, model, cb);
	tiempo = tiempo + dt;

	if isnan(primerTiempo24) && temperaturaNueva(21) > 24
		primerTiempo24 = tiempo;
	endif

	historial(:, paso + 1) = temperaturaNueva;
	if norm(temperaturaNueva - temperatura, inf) < tol
		break;
	endif
	temperatura = temperaturaNueva;
endfor

if paso == maxPasos
	error('No se alcanzo la tolerancia antes de maxPasos.');
endif

errorRel41 = abs((temperaturaNueva - Tex(x)') ./ Tex(x)');
errorRel41(1) = 0;

fprintf('\nPunto 2: lambda = %.1f, dt = %.5f, pasos = %d\n', ...
				model.k*dt/mesh.h^2, dt, paso);
fprintf('Primer instante con T(1/2,t) > 24: t = %.5f\n', primerTiempo24);
fprintf('Tiempo de convergencia: t = %.5f\n', tiempo);

% 3) Valores finales solicitados.
posiciones = [0; 0.25; 0.5; 0.75; 1];
indices = round(posiciones/mesh.h) + 1;
fprintf('\nPunto 3: valores finales\n');
disp([posiciones, temperaturaNueva(indices)]);

% 4) Error relativo respecto de Tex en todos los nodos pedidos.
fprintf('\nPunto 4: error relativo en los nodos pedidos\n');
disp([posiciones, errorRel41(indices)]);

figure;
plot(x, temperaturaNueva, 'o-', x, Tex(x), '--');
grid on;
xlabel('x');
ylabel('T');
legend('Diferencias finitas', 'Solucion exacta');
title('Problema 1: estado estacionario');
