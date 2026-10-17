function [PHI, Q] = fdm2d_implicit(K,F,xnode,neighb,model,dt)
% Descripción: módulo para resolver el sistema lineal de ecuaciones utilizando esquema
% temporal implícito.

% Entrada:
% * K: matriz del sistema (difusión + reacción)
% * F: vector de flujo térmico.
% * xnode: matriz de nodos con pares (x,y) representando las coordenadas de cada nodo 
% de la malla.
% * neighb: matriz de vecindad.
% * model: struct con todos los datos del modelo (constantes, esquema numérico, etc.)
% * dt: paso temporal crítico para método explícito.

% Salida:
% * PHI: matriz solución. Cada elemento del vector representa un valor escalar 
%   asociado a cada nodo de la malla, y su posición dentro del vector depende de
%   cómo se especificó cada nodo en xnode. Cada columna representa una iteración
%   del esquema temporal (en total nit columnas).
% * Q: matriz de flujo de calor. Para cada nodo se halla un vector bidimensional
%   de flujo de calor, representado por un par (Qx,Qy). Cada par de columnas 
%   representa una iteración del esquema temporal (en total 2×nit columnas).
% ----------------------------------------------------------------------
    % Matriz identidad
    I = eye(model.nnodes,model.nnodes);

    % Factor temporal
    A = model.rho*model.cp/dt;

    % Matriz del sistema implícito
    M = A*I + K;

    % Condición inicial
    PHI_n = model.PHI_n;

    % Guardamos la solución inicial
    PHI_vec = PHI_n;

    % Flujo correspondiente a la condición inicial
    Q = fdm2d_flux(PHI_n,neighb,xnode,model.k);
    Q_vec = Q;

    for n = 1:model.maxit

        % Resolvemos el sistema implícito
        PHI = M\(F + A*PHI_n);

        % Error relativo
        err = norm(PHI-PHI_n,2)/norm(PHI,2);

        % Actualizamos
        PHI_n = PHI;

        % Guardamos solución
        PHI_vec = [PHI_vec PHI];

        % Calculamos flujo
        Q = fdm2d_flux(PHI,neighb,xnode,model.k);

        % Guardamos flujo
        Q_vec = [Q_vec Q];

        % Criterio de parada
        if err < model.tol
            disp('Método terminado por tolerancia de error.');
            break;
        end

    end

    if n == model.maxit && err >= model.tol
        disp('Método terminado por límite de iteraciones.');
    end

    PHI = PHI_vec;
    Q = Q_vec;

end