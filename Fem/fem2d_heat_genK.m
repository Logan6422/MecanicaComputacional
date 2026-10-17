function [localK] = fem2d_heat_genK(nodes,kx,ky)
% Descripción: módulo para calcular y evaluar de forma numérica la matriz de
% difusión K. Se utilizan funciones de forma en coordenadas naturales y se 
% resuelve la integral de forma numérica utilizando cuadratura de Gauss.

% Entrada:
% * nodes: nodos (x,y) del elemento. Los elementos admisibles son de 3 o 4 nodos.
% * kx: conductividad térmica orientada en eje-x.
% * ky: conductividad térmica orientada en eje-y.

% Salida:
% * localK: matriz de difusión del elemento (local).
% ----------------------------------------------------------------------

    conductivity = [kx 0; 0 ky];

    if size(nodes,1) == 3
        x1 = nodes(1,1);
        y1 = nodes(1,2);
        x2 = nodes(2,1);
        y2 = nodes(2,2);
        x3 = nodes(3,1);
        y3 = nodes(3,2);
        area = polyarea(nodes(:,1),nodes(:,2));

        b = [y2-y3; y3-y1; y1-y2];
        c = [x3-x2; x1-x3; x2-x1];
        localK = (kx*(b*b') + ky*(c*c'))/(4*area);
    else
        gaussCoordinate = 1/sqrt(3);
        gaussPoints = [-gaussCoordinate -gaussCoordinate;
                        gaussCoordinate -gaussCoordinate;
                        gaussCoordinate gaussCoordinate;
                       -gaussCoordinate gaussCoordinate];
        localK = zeros(4,4);

        for i = 1 : size(gaussPoints,1)
            xi = gaussPoints(i,1);
            eta = gaussPoints(i,2);
            dShapeXi = 0.25*[-(1-eta); 1-eta; 1+eta; -(1+eta)];
            dShapeEta = 0.25*[-(1-xi); -(1+xi); 1+xi; 1-xi];
            jacobian = [dShapeXi'; dShapeEta']*nodes;
            dShape = jacobian\[dShapeXi'; dShapeEta'];
            localK = localK + abs(det(jacobian))*(dShape'*conductivity*dShape);
        end
    end

end