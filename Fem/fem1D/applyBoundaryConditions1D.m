function [K,F,dirichletNodes,dirichletValues] = ...
        applyBoundaryConditions1D(K,F,cb)

    dirichletNodes = [];
    dirichletValues = [];

    for side = 1:2
        if side == 1
            node = 1;
        else
            node = length(F);
        end
        boundaryType = cb(side,1);
        boundaryValue = cb(side,2);
        externalTemperature = cb(side,3);

        if boundaryType == 1
            F = F - K(:,node)*boundaryValue;

            K(node,:) = 0;
            K(:,node) = 0;
            K(node,node) = 1;
            F(node) = boundaryValue;

            dirichletNodes = [dirichletNodes; node];
            dirichletValues = [dirichletValues; boundaryValue];

        elseif boundaryType == 2
            F(node) = F(node) - boundaryValue;

        elseif boundaryType == 3
            K(node,node) = K(node,node) + boundaryValue;
            F(node) = F(node) + boundaryValue*externalTemperature;

        else
            error('Tipo de condicion de borde invalido.');
        end
    end
end