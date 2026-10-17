function [K,F] = fdm2d_gen_system(K,F,xnode,neighb,k,c,G)
% Descripción: módulo para ensamblar los términos difusivo, reactivo y fuente de
% todos los nodos de la malla, generando el stencil adecuado dependiendo de si es
% un nodo interior o de frontera.

% Entrada:
% * K: matriz del sistema (difusión + reacción)
% * F: vector de flujo térmico.
% * xnode: matriz de pares (x,y) representando cada nodo de la malla.
% * neighb: matriz de vecindad.
% * k: conductividad térmica del material. Es un vector que permite representar k(x,y).
% * c: constante de reacción del material. Es un vector que permite representar c(x,y).
% * G: fuente de calor. Es un vector que permite representar G(x,y).

% Salida:
% * K: matriz del sistema (difusión + reacción) con modificaciones luego del ensamble.
% * F: vector de flujo térmico con modificaciones luego del ensamble.
%----------------------------------------------------------------------
for P = 1:size(xnode,1)
        %nodos vecinos
        S = neighb(P,1);
        E = neighb(P,2);
        N = neighb(P,3);
        W = neighb(P,4);
        %distancia a cada nodo (malla no uniforme)
        ds = 0;
        de = 0;
        dn = 0;
        dw = 0;
       
        %si existe el nodo vecino (!= -1) calculo las distancias
        if(S ~= -1)
            ds = abs(xnode(S,2) - xnode(P,2));
        end
        if(E ~= -1)
            de = abs(xnode(E,1) - xnode(P,1));
        end
        if(N ~= -1)
            dn = abs(xnode(N,2) - xnode(P,2));
        end
        if(W ~= -1)
            dw = abs(xnode(W,1) - xnode(P,1));
        end

        %si no existe el nodo vecino (= -1) 
        if(W == -1)
            ax = 2/(de*de);
            bx = -2/(de*de);
            cx = 0;
        elseif(E == -1)
            ax = 0;
            bx = -2/(dw*dw);
            cx = 2/(dw*dw);
        else
            ax = 2/(de*(de+dw));
            bx = -2/(de*dw);
            cx = 2/(dw*(de+dw));
        end
        if(S == -1)
            ay = 2/(dn*dn);
            by = -2/(dn*dn);
            cy = 0;
        elseif(N == -1)
            ay = 0;
            by = -2/(ds*ds);
            cy = 2/(ds*ds);
        else
            ay = 2/(dn*(dn+ds));
            by = -2/(dn*ds);
            cy = 2/(ds*(dn+ds));
        end
        if(E ~= -1)
            K(P,E) = -ax;
        end

        if(W ~= -1)
            K(P,W) = -cx;
        end
        if(N ~= -1)
            K(P,N) = -ay;
        end
        if(S ~= -1)
            K(P,S) = -cy;
        end

        %agrego los coef a K y agrego la fuente
        K(P,P) = -bx - by + c(P)/k(P);
        F(P) = G(P)/k(P);
end
