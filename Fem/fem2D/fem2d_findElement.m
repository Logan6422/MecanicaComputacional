function element = fem2d_findElement(xnode,icone,xp,yp)
    element = 0;
    tolerance = 1e-12;

    for e = 1:size(icone,1)
        nodes = xnode(icone(e,:),:);
        if xp >= min(nodes(:,1))-tolerance && ...
           xp <= max(nodes(:,1))+tolerance && ...
           yp >= min(nodes(:,2))-tolerance && ...
           yp <= max(nodes(:,2))+tolerance
            element = e;
            return;
        endif
    endfor

    error('El punto (%g,%g) no pertenece a la malla.',xp,yp);
endfunction