function [Qn] = fem2d_heat_flux(xnode,icone,model,PHI)
% Descripción: módulo calcular el flujo de calor en todo el dominio. Se aplica
% la Ley de Fourier y se evalúa como fluye el calor en todos los puntos (nodos)
% del dominio.

% Entrada:
% * K: matriz del sistema (difusión + reacción)
% * F: vector de flujo térmico.
% * PHI: vector solución. Cada elemento del vector representa un valor escalar
%   asociado a cada nodo de la malla, y su posición dentro del vector depende 
%   de cómo se especificó cada nodo en xnode.

% Salida:
% * Q: vector de flujo de calor.
% ----------------------------------------------------------------------

    Q = zeros(model.nnodes,3);
    conductivity = [model.kx 0; 0 model.ky];

    for e = 1 : model.nelem
        if icone(e,4) == -1
            elementNodes = icone(e,1:3);
        else
            elementNodes = icone(e,:);
        end

        nodes = xnode(elementNodes,:);
        temperature = PHI(elementNodes);

        if size(nodes,1) == 3
            jacobian = [nodes(2,1)-nodes(1,1) nodes(2,2)-nodes(1,2);
                        nodes(3,1)-nodes(1,1) nodes(3,2)-nodes(1,2)];
            dShapeNatural = [-1 1 0; -1 0 1];
            gradient = jacobian\dShapeNatural;
            elementFlux = (-conductivity*gradient*temperature)';

            for i = 1 : 3
                node = elementNodes(i);
                Q(node,1:2) = Q(node,1:2) + elementFlux;
                Q(node,3) = Q(node,3) + 1;
            end
        else
            naturalPoints = [-1 -1; 1 -1; 1 1; -1 1];

            for i = 1 : 4
                xi = naturalPoints(i,1);
                eta = naturalPoints(i,2);
                dShapeXi = 0.25*[-(1-eta); 1-eta; 1+eta; -(1+eta)];
                dShapeEta = 0.25*[-(1-xi); -(1+xi); 1+xi; 1-xi];
                jacobian = [dShapeXi'; dShapeEta']*nodes;
                gradient = jacobian\[dShapeXi'; dShapeEta'];
                shape = 0.25*[(1-xi)*(1-eta);
                              (1+xi)*(1-eta);
                              (1+xi)*(1+eta);
                              (1-xi)*(1+eta)];
                elementFlux = (-conductivity*gradient*temperature)';
                node = elementNodes(i);
                Q(node,1:2) = Q(node,1:2) + elementFlux;
                Q(node,3) = Q(node,3) + 1;
            end
        end
    end

    Qn = zeros(model.nnodes,2);
    usedNodes = Q(:,3) > 0;
    Qn(usedNodes,:) = Q(usedNodes,1:2)./Q(usedNodes,3);
end

