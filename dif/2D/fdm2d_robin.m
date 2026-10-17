function [K,F] = fdm2d_robin(K,F,xnode,neighb,ROB)
% Descripción: módulo para calcular y ensamblar las contribuciones de nodos 
% pertenecientes a fronteras de tipo Robin.

% Entrada:
% * K: matriz del sistema (difusión + reacción)
% * F: vector de flujo térmico.
% * xnode: matriz de nodos con pares (x,y) representando las coordenadas de 
%   cada nodo de la malla.
% * neighb: matriz de vecindad.
% * ROB: matriz con la información sobre la frontera de tipo Robin.
%   - Columna 1: índice del nodo donde se aplica la condición de borde.
%   - Columna 2: valor de coeficiente de calor (h)
%   - Columna 3: valor de temperatura de referencia (phi_inf).
%   - Columna 4: dirección y sentido del flujo:
%     (1) Flujo en dirección eje-y, sentido negativo (S – South – Sur)
%     (2) Flujo en dirección eje-x, sentido positivo (E – East – Este)
%     (3) Flujo en dirección eje-y, sentido positivo (N – North – Norte)
%     (4) Flujo en dirección eje-x, sentido negativo (W – West – Oeste)

% Salida:
% * K: matriz del sistema (difusión + reacción) con modificaciones luego de 
%   aplicar la condición de borde.
% * F: vector de flujo térmico con modificaciones luego de aplicar la condición 
%   de borde.
% ----------------------------------------------------------------------
M = size(ROB,1);

for n = 1:M

    % Datos de la condición Robin
    P = ROB(n,1);
    h = ROB(n,2);
    phi_inf = ROB(n,3);
    direccion = ROB(n,4);

    % Busco el nodo interior y calculo la distancia
    if direccion == 1
        % Sur -> vecino interior hacia Norte
        vecino = neighb(P,3);
        dist = abs(xnode(P,2) - xnode(vecino,2));

    elseif direccion == 2
        % Este -> vecino interior hacia Oeste
        vecino = neighb(P,4);
        dist = abs(xnode(P,1) - xnode(vecino,1));

    elseif direccion == 3
        % Norte -> vecino interior hacia Sur
        vecino = neighb(P,1);
        dist = abs(xnode(P,2) - xnode(vecino,2));

    elseif direccion == 4
        % Oeste -> vecino interior hacia Este
        vecino = neighb(P,2);
        dist = abs(xnode(P,1) - xnode(vecino,1));
    end

    % Aporte de Robin a la diagonal de K
    K(P,P) = K(P,P) + 2*h/dist;

    % Aporte de Robin al vector F
    F(P) = F(P) + 2*h*phi_inf/dist;

end

end