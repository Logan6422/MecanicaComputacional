#Datos del problema de calor
#model.k = 2;
#model.c = 2;
#model.rho = 1;
#model.cp = 1;
#model.G = 75;

#Datos de la subdivision
#mesh.xi = 0;
#mesh.xf = 10;
#mesh.h = 0.5;

#cb_Ej = [1 10 -1]

function T = difFinitas(mesh, model, cb, et)
  xnode = mesh.xi:mesh.h:mesh.xf; #subdivision 1D

  N = length(xnode);
  A = zeros(N,N);
  b = zeros(N,1);
  for i = 2:N-1
    #nodo interior
    A(i,i-1) = -1;
    A(i,i) = 2 + (mesh.h^2*model.c/model.k);
    A(i,i+1) = -1;

    b(i) = mesh.h^2 * model.G / model.k;
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

      b(1) = (2*mesh.h * q / model.k) ...
          + (mesh.h^2 * model.G / model.k);
    case 3
      #Robin
      hc = cb(1,2);
      Tinf = cb(1,3);

      A(1,1) = 2 + (mesh.h^2 * model.c / model.k) ...
            + (2 * mesh.h * hc / model.k);

      A(1,2) = -2;

      b(1) = (mesh.h^2 * model.G / model.k) ...
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

      b(N) = -(2*mesh.h*q/model.k) ...
             + (mesh.h^2*model.G/model.k);

    case 3
      #Robin
      hc = cb(2,2);
      Tinf = cb(2,3);

      A(N,N-1) = -2;

      A(N,N) = 2 + (mesh.h^2 * model.c / model.k) ...
               + (2*mesh.h*hc/model.k);

      b(N) = (mesh.h^2*model.G/model.k) ...
             + (2*mesh.h*hc*Tinf/model.k);
  endswitch

  T = A\b;
endfunction

