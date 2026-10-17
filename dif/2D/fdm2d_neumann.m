function [F] = fdm2d_neumann(F,xnode,neighb,NEU)
% Descripción: módulo para calcular y ensamblar las contribuciones de nodos 
% pertenecientes a fronteras de tipo Neumann.

% Entrada:
% * F: vector de flujo térmico.
% * xnode: matriz de nodos con pares (x,y) representando las coordenadas de cada 
%   nodo de la malla.
% * neighb: matriz de vecindad.
% * NEU: matriz con la información sobre la frontera de tipo Neumann. 
%   - Columna 1: índice del nodo donde se aplica la condición de borde.
%   - Columna 2: valor de flujo térmico (q) asociado al lado del elemento.
%   - Columna 3: dirección y sentido del flujo:
%     (1) Flujo en dirección eje-y, sentido negativo (S – South - Sur)
%     (2) Flujo en dirección eje-x, sentido positivo (E – East - Este)
%     (3) Flujo en dirección eje-y, sentido positivo (N – North - Norte)
%     (4) Flujo en dirección eje-x, sentido negativo (W – West – Oeste)

% Salida:
% * F: vector de flujo térmico con modificaciones luego de aplicar la condición de borde.
% ----------------------------------------------------------------------
for i = 1:size(NEU,1)%Recorro nodos Neumann
    nodo = NEU(i,1);        
    q = NEU(i,2); 
    direccion = NEU(i,3); %Direccion neumann

    if direccion == 1            %S
        vecino = neighb(nodo,3);   %Vecino norte
        h = abs(xnode(nodo,2) - xnode(vecino,2));  %Calculo de distancia

        F(nodo) = F(nodo) - 2*q/h;

    elseif direccion == 2          %E
        vecino = neighb(nodo,4);   %Vecino oeste
        h = abs(xnode(nodo,1) - xnode(vecino,1));

        F(nodo) = F(nodo) - 2*q/h;

    elseif direccion == 3          %N
        vecino = neighb(nodo,1);   %Vecino sur
        h = abs(xnode(nodo,2) - xnode(vecino,2));

        F(nodo) = F(nodo) - 2*q/h;

    elseif direccion == 4          %O
        vecino = neighb(nodo,2);   %Vecino este
        h = abs(xnode(nodo,1) - xnode(vecino,1));

        F(nodo) = F(nodo) - 2*q/h;
    end
end