function [K,F] = fdm2d_dirichlet(K,F,DIR)
% Descripción: módulo para calcular y ensamblar las contribuciones de nodos 
% pertenecientes a fronteras de tipo Dirichlet.

% Entrada:
% * K: matriz del sistema (difusión + reacción)
% * F: vector de flujo térmico.
% * DIR: matriz con la información sobre la frontera de tipo Dirchlet.
%   - Columna 1: número de nodo.
%   - Columna 2: valor en ese nodo (escalar)

% Salida:
% * K: matriz del sistema (difusión + reacción) luego de realizar las simplificaciones
%   que surgen de aplicar la condición de borde Dirichlet.
% * F: vector de flujo térmico luego de realizar las simplificaciones que surgen de
%   aplicar la condición de borde Dirichlet.
% ----------------------------------------------------------------------
    
    for j = 1:size(DIR,1) %recorro los nodos dirichlet
        nodo = DIR(j,1); %obtengo el nodo
        valor = DIR(j,2); %obtengo el valor en el nodo
     
        K(nodo,:) = 0; %cero en todos los nodos  de la fila nodo
        K(nodo,nodo) = 1; %uno en el nodo correspondiente
    
        F(nodo) = valor; %valor en el termino independiente
    
    end
end

