function [K,F] = fem2d_heat_robin(K,F,ROB,xnode)
% Descripción: módulo para calcular y ensamblar las contribuciones de nodos
% pertenecientes a fronteras de tipo Robin.

% Entrada:
% * K: matriz del sistema (difusión + reacción)
% * F: vector de flujo térmico.
% * ROB: matriz con la información sobre la frontera de tipo Robin.
%   - Columnas 1-2: dos nodos contiguos formando un lado de un elemento.
%   - Columna 3: valor de coeficiente de calor (h).
%   - Columna 4: valor de temperatura de referencia (phi_inf).
% * xnode: matriz de nodos con pares (x,y) representando las coordenadas de cada nodo de la malla.

% Salida:
% * K: matriz del sistema (difusión + reacción) con modificaciones luego de
%   aplicar la condición de borde.
% * F: vector de flujo térmico con modificaciones luego de aplicar la 
%   condición de borde.
% ----------------------------------------------------------------------
	for i = 1 : size(ROB,1)
		side = ROB(i,1:2);
		h = ROB(i,3);
		referenceTemperature = ROB(i,4);
		edgeLength = norm(xnode(side(1),:) - xnode(side(2),:),2);

		K(side,side) = K(side,side) + (h*edgeLength/6)*[2 1; 1 2];
		F(side) = F(side) + (h*referenceTemperature*edgeLength/2)*[1; 1];
	end
end