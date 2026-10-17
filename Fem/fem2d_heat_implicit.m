function [PHI_vec,Q_vec] = fem2d_heat_implicit(K,C,F,xnode,icone,model,dt)
% Descripción: módulo para resolver el sistema lineal de ecuaciones utilizando
% esquema temporal implícito.

% Entrada:
% * K: matriz del sistema (difusión + reacción).
% * C: matriz de masa del sistema (sin escalar por la constante de reacción).
% * F: vector de flujo térmico.
% * xnode: matriz de nodos con pares (x,y) representando las coordenadas de 
%   cada nodo de la malla.
% * icone: matriz de conectividad. Indica los 3 ó 4 nodos que integran el 
%   elemento, recorridos en cualquier orden pero en sentido antihorario. 
%   En caso de elementos triangulares, la cuarta columna siempre es -1.
% * model: struct con todos los datos del modelo (constantes, esquema numérico, etc.)
% * dt: paso temporal arbitrario para método implícito.

% Salida:
% * PHI: matriz solución. Cada elemento del vector representa un valor 
%   escalar asociado a cada nodo de la malla, y su posición dentro del vector 
%   depende de cómo se especificó cada nodo en xnode. Cada columna representa 
%   una iteración del esquema temporal (en total nit columnas).
% * Q: matriz de flujo de calor. Para cada nodo se halla un vector bidimensional
%   de flujo de calor, representado por un par (Qx,Qy). Cada par de columnas 
%   representa una iteración del esquema temporal (en total 2×nit columnas).
% ----------------------------------------------------------------------

    auxA = dt/(model.rho*model.cp);
    systemMatrix = C + auxA*K;
    PHI = model.PHI_n;
    PHI_n = model.PHI_n;
    PHI_vec = PHI;
    Q_vec = zeros(model.nnodes,2);

    for n = 1 : model.maxit
        PHI = systemMatrix\(C*PHI_n + auxA*F);
        denominator = norm(PHI,2);

        if denominator == 0
            err = norm(PHI-PHI_n,2);
        else
            err = norm(PHI-PHI_n,2)/denominator;
        end

        PHI_n = PHI;
        PHI_vec = [PHI_vec PHI];
        Q = fem2d_heat_flux(xnode,icone,model,PHI);
        Q_vec = [Q_vec Q];

        if err < model.tol
            disp('Metodo terminado por tolerancia de error.');
            return;
        end
    end
end
