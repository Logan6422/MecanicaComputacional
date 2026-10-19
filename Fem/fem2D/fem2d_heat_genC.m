function [localC] = fem2d_heat_genC(nodes)
% Descripción: módulo para calcular y evaluar de forma numérica la matriz de
% masa C. Se utilizan funciones de forma en coordenadas naturales y se 
% resuelve la integral de forma numérica utilizando cuadratura de Gauss.

% Entrada:
% * nodes: nodos (x,y) del elemento. Los elementos admisibles son de 3 o 4 nodos.

% Salida:
% * localC: matriz de masa para elemento (local).
% ----------------------------------------------------------------------

    if size(nodes,1) == 3
        area = polyarea(nodes(:,1),nodes(:,2));
        gaussPoints = [1/6 1/6; 2/3 1/6; 1/6 2/3];
        localC = zeros(3,3);

        for i = 1 : size(gaussPoints,1)
            xi = gaussPoints(i,1);
            eta = gaussPoints(i,2);
            shape = [1-xi-eta; xi; eta];
            localC = localC + (area/3)*(shape*shape');
        end
    else
        gaussCoordinate = 1/sqrt(3);
        gaussPoints = [-gaussCoordinate -gaussCoordinate;
                        gaussCoordinate -gaussCoordinate;
                        gaussCoordinate gaussCoordinate;
                       -gaussCoordinate gaussCoordinate];
        localC = zeros(4,4);

        for i = 1 : size(gaussPoints,1)
            xi = gaussPoints(i,1);
            eta = gaussPoints(i,2);
            shape = 0.25*[(1-xi)*(1-eta);
                          (1+xi)*(1-eta);
                          (1+xi)*(1+eta);
                          (1-xi)*(1+eta)];
            dShapeXi = 0.25*[-(1-eta); 1-eta; 1+eta; -(1+eta)];
            dShapeEta = 0.25*[-(1-xi); -(1+xi); 1+xi; 1-xi];
            jacobian = [dShapeXi'; dShapeEta']*nodes;
            localC = localC + (shape*shape')*abs(det(jacobian));
        end
    end

end