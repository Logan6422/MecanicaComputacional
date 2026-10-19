function [localF] = sourceVector1D(xLeft,xRight,G,numberOfPoints)
% Calcula el vector local de fuente para G constante o G(x) variable.
%
% Fe = integral(N' * G(x) dx)
% numberOfPoints indica la cantidad de puntos de Gauss-Legendre.

    if nargin < 4
        numberOfPoints = 2;
    end
    [naturalPoints,weights] = gaussLegendreRule(numberOfPoints);

    elementLength = xRight - xLeft;
    localF = zeros(2,1);

    for i = 1 : length(naturalPoints)
        naturalCoordinate = naturalPoints(i);
        shape = [(1-naturalCoordinate)/2;
                 (1+naturalCoordinate)/2];
        x = ((1-naturalCoordinate)*xLeft + ...
             (1+naturalCoordinate)*xRight)/2;

        if isa(G,'function_handle')
            sourceValue = G(x);
        else
            sourceValue = G;
        end

        if ~isscalar(sourceValue)
            error('G debe devolver un valor escalar para cada coordenada x.');
        end

        localF = localF + shape*sourceValue*weights(i)*(elementLength/2);
    end
end

function [naturalPoints,weights] = gaussLegendreRule(numberOfPoints)
    switch numberOfPoints
        case 1
            naturalPoints = 0;
            weights = 2;
        case 2
            naturalPoints = [-1/sqrt(3), 1/sqrt(3)];
            weights = [1, 1];
        case 3
            naturalPoints = [-sqrt(3/5), 0, sqrt(3/5)];
            weights = [5/9, 8/9, 5/9];
        case 4
            naturalPoints = [-sqrt((3+2*sqrt(6/5))/7), ...
                             -sqrt((3-2*sqrt(6/5))/7), ...
                              sqrt((3-2*sqrt(6/5))/7), ...
                              sqrt((3+2*sqrt(6/5))/7)];
            weights = [(18-sqrt(30))/36, (18+sqrt(30))/36, ...
                       (18+sqrt(30))/36, (18-sqrt(30))/36];
        otherwise
            error('La cuadratura de Gauss admite entre 1 y 4 puntos.');
    end
end
