function [K,C,F] = assembleSystem1D(xnode,model)
% Ensambla las matrices y el vector globales del problema 1D.

    if isfield(model,'quadraturePoints')
        numberOfPoints = model.quadraturePoints;
    else
        numberOfPoints = 2;
    end

    numberOfNodes = length(xnode);
    K = zeros(numberOfNodes);
    C = zeros(numberOfNodes);
    F = zeros(numberOfNodes,1);

    for element = 1 : numberOfNodes-1
        elementNodes = [element element+1];
        elementLength = xnode(element+1) - xnode(element);

        localK = model.k/elementLength*[1 -1; -1 1];
        localC = elementLength/6*[2 1; 1 2];
        localF = sourceVector1D(xnode(element),xnode(element+1),...
                    model.G,numberOfPoints);

        K(elementNodes,elementNodes) = ...
            K(elementNodes,elementNodes) + localK;
        C(elementNodes,elementNodes) = ...
            C(elementNodes,elementNodes) + localC;
        F(elementNodes) = F(elementNodes) + localF;
    end

    K = K + model.c*C;
end
