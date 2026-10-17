function [components] = mostrar_componentes(xnode,model,cb)
% Muestra y devuelve los componentes principales del FEM 1D.
%
% Uso:
%   components = mostrar_componentes(xnode,model,cb)
%
% La estructura devuelta permite consultar, por ejemplo:
%   components.Kdifusion
%   components.Cmasa
%   components.F
%   components.KconBorde
%   components.FconBorde

    xnode = xnode(:);
    numberOfNodes = length(xnode);

    if numberOfNodes < 2 || any(diff(xnode) <= 0)
        error('xnode debe tener al menos dos nodos ordenados.');
    end
    if size(cb,1) ~= 2 || size(cb,2) ~= 3
        error('cb debe ser una matriz de 2x3.');
    end

    [Ktotal,Cmasa,F] = assembleSystem1D(xnode,model);
    Kdifusion = Ktotal - model.c*Cmasa;

    [KconBorde,FconBorde,dirichletNodes,dirichletValues] = ...
        applyBoundaryConditions1D(Ktotal,F,cb);

    elementLength = xnode(2) - xnode(1);
    localK = model.k/elementLength*[1 -1; -1 1];
    localC = elementLength/6*[2 1; 1 2];
    localF = sourceVector1D(xnode(1),xnode(2),model.G);

    components = struct();
    components.xnode = xnode;
    components.cb = cb;
    components.Kdifusion = Kdifusion;
    components.Cmasa = Cmasa;
    components.Kreaccion = model.c*Cmasa;
    components.Ktotal = Ktotal;
    components.F = F;
    components.KconBorde = KconBorde;
    components.FconBorde = FconBorde;
    components.localK = localK;
    components.localC = localC;
    components.localF = localF;
    components.dirichletNodes = dirichletNodes;
    components.dirichletValues = dirichletValues;

    fprintf('\n========== COMPONENTES FEM 1D ==========\n');
    fprintf('Nodos: %d | Elementos: %d\n',numberOfNodes,numberOfNodes-1);
    fprintf('xnode =\n');
    disp(xnode');
    fprintf('model =\n');
    disp(model);
    fprintf('cb =\n');
    disp(cb);
    fprintf('Matriz local de difusion Ke = (k/he)*[1 -1; -1 1]\n');
    disp(localK);
    fprintf('Matriz local de masa Ce = (he/6)*[2 1; 1 2]\n');
    disp(localC);
    fprintf('Vector local de fuente Fe = integral(N''*G(x)*dx)\n');
    disp(localF);
    fprintf('Matriz global de difusion Kdifusion = ensamblaje de Ke\n');
    disp(Kdifusion);
    fprintf('Matriz global de masa Cmasa = ensamblaje de Ce\n');
    disp(Cmasa);
    fprintf('Matriz global de reaccion Kreaccion = c*Cmasa\n');
    disp(components.Kreaccion);
    fprintf('Matriz global antes de bordes Ktotal = Kdifusion + c*Cmasa\n');
    disp(Ktotal);
    fprintf('Vector global de fuente F = ensamblaje de Fe\n');
    disp(F);
    fprintf('Matriz despues de bordes KconBorde = aplicar condiciones a Ktotal\n');
    disp(KconBorde);
    fprintf('Vector despues de bordes FconBorde = aplicar condiciones a F\n');
    disp(FconBorde);
    fprintf('==========================================\n\n');
end
