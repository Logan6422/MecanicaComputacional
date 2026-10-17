function [localF] = sourceVector1D(xLeft,xRight,G)
% Calcula el vector local de fuente para G constante o G(x) variable.
%
% Fe = integral(N' * G(x) dx)
% Se usa cuadratura de Gauss de dos puntos.

    elementLength = xRight - xLeft;
    gaussPoint = 1/sqrt(3);
    naturalPoints = [-gaussPoint, gaussPoint];
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

        localF = localF + shape*sourceValue*(elementLength/2);
    end
end
