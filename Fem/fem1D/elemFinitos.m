function [T] = elemFinitos(xnode, model, cb, et)
% Coordina las etapas del calculo FEM 1D.

    xnode = xnode(:);
    numberOfNodes = length(xnode);
    validateInputs(xnode,model,cb,numberOfNodes,et);
    [K,C,F] = assembleSystem1D(xnode,model);
    [K,F,dirichletNodes,dirichletValues] = ...
        applyBoundaryConditions1D(K,F,cb);

    if et ~= 0 && et ~= 1
        T = solveSteady1D(K,F);
        return;
    end

    if et == 1
        T = solveImplicit1D(C,K,F,model,dirichletNodes,dirichletValues);
    else
        T = solveExplicit1D(C,K,F,model,dirichletNodes,dirichletValues);
    end
end

function validateInputs(xnode,model,cb,numberOfNodes,et)
    if numberOfNodes < 2
        error('xnode debe contener al menos dos nodos.');
    end
    if size(cb,1) ~= 2 || size(cb,2) ~= 3
        error('cb debe ser una matriz de 2x3.');
    end
    if any(diff(xnode) <= 0)
        error('Los nodos de xnode deben estar ordenados y ser distintos.');
    end
    requiredFields = {'k','c','rho','cp','G'};
    for i = 1 : length(requiredFields)
        if ~isfield(model,requiredFields{i})
            error('Falta el campo model.%s.',requiredFields{i});
        end
    end
    if et == 0 || et == 1
        temporalFields = {'dt','maxit','tol','T0'};
        for i = 1 : length(temporalFields)
            if ~isfield(model,temporalFields{i})
                error('Para un esquema temporal se requiere model.%s.',...
                    temporalFields{i});
            end
        end
        if length(model.T0) ~= numberOfNodes
            error('model.T0 debe tener un valor por nodo.');
        end
    end
end
