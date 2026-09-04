# Si G es constante:
# model.G = 75;

# Si G depende de x:
# model.G = @(x) 10*x + 5;

# Datos del problema de calor
# model.k = 2;
# model.c = 2;
# model.rho = 1;
# model.cp = 1;
# model.rhoCp = model.rho * model.cp;
# model.G = 75;

# Datos de la subdivision
# mesh.xi = 0;
# mesh.xf = 10;
# mesh.h = 0.5;

# cb_Ej = [1 10 -1]


function T = difFinitas(mesh, model, cb, et)
  xnode = mesh.xi:mesh.h:mesh.xf; #subdivision 1D

  N = length(xnode);
  A = zeros(N,N);
  b = zeros(N,1);

  if model.transitorio == 0
    for i = 2:N-1
      #nodo interior
      A(i,i-1) = -1;
      A(i,i) = 2 + (mesh.h^2*model.c/model.k);
      A(i,i+1) = -1;

      # G puede ser constante o funcion de x
      if isa(model.G, "function_handle")
        G = model.G(xnode(i));
      else
        G = model.G;
      endif

      b(i) = mesh.h^2 * G / model.k;
    endfor

    #CONDICIONES DE BORDE
    #lado izquierdo
    switch cb(1,1)
      case 1
        #Dirichlet
        A(1,1) = 1;
        b(1) = cb(1,2);

      case 2
        #Neumann
        q = cb(1,2);

        A(1,1) = 2 + (mesh.h^2 * model.c / model.k);
        A(1,2) = -2;

        if isa(model.G, "function_handle")
          G = model.G(xnode(1));
        else
          G = model.G;
        endif

        b(1) = (2*mesh.h * q / model.k) ...
            + (mesh.h^2 * G / model.k);

      case 3
        #Robin
        hc = cb(1,2);
        Tinf = cb(1,3);

        A(1,1) = 2 + (mesh.h^2 * model.c / model.k) ...
              + (2 * mesh.h * hc / model.k);

        A(1,2) = -2;

        if isa(model.G, "function_handle")
          G = model.G(xnode(1));
        else
          G = model.G;
        endif

        b(1) = (mesh.h^2 * G / model.k) ...
              + (2*mesh.h * hc * Tinf / model.k);
    endswitch


    #lado derecho
    switch cb(2,1)
      case 1
        #Dirichlet
        A(N,N) = 1;
        b(N) = cb(2,2);

      case 2
        #Neumann
        q = cb(2,2);

        A(N,N-1) = -2;
        A(N,N) = 2 + (mesh.h^2 * model.c / model.k);

        if isa(model.G, "function_handle")
          G = model.G(xnode(N));
        else
          G = model.G;
        endif

        b(N) = -(2*mesh.h*q/model.k) ...
             + (mesh.h^2*G/model.k);

      case 3
        #Robin
        hc = cb(2,2);
        Tinf = cb(2,3);

        A(N,N-1) = -2;

        A(N,N) = 2 + (mesh.h^2 * model.c / model.k) ...
                 + (2*mesh.h*hc/model.k);

        if isa(model.G, "function_handle")
          G = model.G(xnode(N));
        else
          G = model.G;
        endif

        b(N) = (mesh.h^2 * G / model.k) ...
             + (2*mesh.h*hc*Tinf/model.k);
    endswitch

    T = A\b;

  else
    # PARTE TRANSITORIA
    C = model.rhoCp;
    dt = model.dt;
    theta = model.tita;
    T = model.T0;

    # Cantidad de pasos temporales
    cantPasos = round(model.tf / dt);

    for paso = 1:cantPasos
      # r representa la difusion durante un paso temporal
      r = model.k * dt / (C * mesh.h^2);

      # s representa el termino c durante un paso temporal
      s = model.c * dt / C;

      # Reiniciamos la matriz y el vector en cada paso
      A = zeros(N,N);
      b = zeros(N,1);


      # NODOS INTERIORES
      # theta = 0   -> Euler explicito
      # theta = 1   -> Euler hacia atras
      # theta = 1/2 -> Crank-Nicolson
      for i = 2:N-1
        A(i,i-1) = -theta * r;
        A(i,i) = 1 + theta * (2*r + s);
        A(i,i+1) = -theta * r;

        b(i) = (1-theta) * r * T(i-1) ...
             + (1 - (1-theta)*(2*r+s)) * T(i) ...
             + (1-theta) * r * T(i+1);

        # G puede ser constante o funcion de x
        if isa(model.G, "function_handle")
          G = model.G(xnode(i));
        else
          G = model.G;
        endif

        b(i) = b(i) + dt * G / C;
      endfor


      # CONDICIONES DE BORDE
      # LADO IZQUIERDO
      switch cb(1,1)
        case 1
          #Dirichlet
          A(1,1) = 1;
          b(1) = cb(1,2);
        case 2
          #Neumann
          q = cb(1,2);
          A(1,1) = 1 + theta * (2*r + s);
          A(1,2) = -2 * theta * r;

          b(1) = (1 - theta) * 2*r * T(2) ...
               + (1 - (1-theta)*(2*r+s)) * T(1);
          if isa(model.G, "function_handle")
            G = model.G(xnode(1));
          else
            G = model.G;
          endif

          b(1) = b(1) ...
               + dt * G / C ...
               + 2 * dt * q / (C * mesh.h);

        case 3
          #Robin
          hc = cb(1,2);
          Tinf = cb(1,3);

          A(1,1) = 1 + theta * (2*r + s) ...
                  + theta * (2*dt*hc/(C*mesh.h));

          A(1,2) = -2 * theta * r;

          b(1) = (1-theta) * 2*r * T(2) ...
               + (1 - (1-theta)*(2*r+s)) * T(1);

          if isa(model.G, "function_handle")
            G = model.G(xnode(1));
          else
            G = model.G;
          endif

          b(1) = b(1) ...
               + dt * G / C ...
               + (1-theta) * 2*dt*hc/(C*mesh.h) * Tinf;
      endswitch


      # LADO DERECHO
      switch cb(2,1)

        case 1
          #Dirichlet
          A(N,N) = 1;
          b(N) = cb(2,2);

        case 2
          #Neumann
          q = cb(2,2);

          A(N,N-1) = -2 * theta * r;
          A(N,N) = 1 + theta * (2*r + s);

          b(N) = (1-theta) * 2*r * T(N-1) ...
               + (1 - (1-theta)*(2*r+s)) * T(N);

          if isa(model.G, "function_handle")
            G = model.G(xnode(N));
          else
            G = model.G;
          endif

          b(N) = b(N) ...
               + dt * G / C ...
               - 2 * dt * q / (C * mesh.h);

        case 3
          #Robin
          hc = cb(2,2);
          Tinf = cb(2,3);

          A(N,N-1) = -2 * theta * r;

          A(N,N) = 1 + theta * (2*r + s) ...
                  + theta * (2*dt*hc/(C*mesh.h));

          b(N) = (1-theta) * 2*r * T(N-1) ...
               + (1 - (1-theta)*(2*r+s)) * T(N);

          if isa(model.G, "function_handle")
            G = model.G(xnode(N));
          else
            G = model.G;
          endif

          b(N) = b(N) ...
               + dt * G / C ...
               + (1-theta) * 2*dt*hc/(C*mesh.h) * Tinf;
      endswitch


      # RESOLVEMOS PARA OBTENER T(n+1)
      T = A\b;

      # Forzamos nuevamente Dirichlet.
      # Esto corrige cualquier pequeno error numerico.
      if cb(1,1) == 1
        T(1) = cb(1,2);
      endif

      if cb(2,1) == 1
        T(N) = cb(2,2);
      endif

    endfor

  endif

endfunction
