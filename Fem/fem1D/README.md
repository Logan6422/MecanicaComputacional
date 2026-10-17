# FEM 1D para transferencia de calor

La funcion principal es:

```matlab
T = elemFinitos(xnode,model,cb,et);
```

El codigo esta separado en etapas:

- `elemFinitos.m`: coordina el calculo.
- `assembleSystem1D.m`: ensambla `K`, `C` y `F`.
- `applyBoundaryConditions1D.m`: aplica Dirichlet, Neumann y Robin.
- `solveSteady1D.m`: resuelve el caso estacionario.
- `solveExplicit1D.m`: avanza con Euler explicito.
- `solveImplicit1D.m`: avanza con Euler implicito.
- `graficaEvolucion.m`: muestra las columnas de la solucion temporal.

## Malla y modelo

`xnode` es un vector ordenado de coordenadas. Los elementos son los segmentos:

```text
[xnode(1), xnode(2)], [xnode(2), xnode(3)], ...
```

Campos obligatorios de `model`:

```matlab
model.k   % conductividad
model.c   % coeficiente del termino de reaccion -c*T
model.rho % densidad
model.cp  % calor especifico
model.G   % fuente constante o funcion G(x)
```

Para una fuente variable se usa una funcion anonima de Octave:

```matlab
model.G = @(x) 100*(x-3).^2;
```

El vector local de fuente se calcula como `Fe = integral(N''*G(x)*dx)` mediante
cuadratura de Gauss de dos puntos. Para una fuente constante, tambien se puede
usar simplemente `model.G = 100`.

Para los esquemas temporales tambien se necesitan:

```matlab
model.dt
model.maxit
model.tol
model.T0   % vector columna o fila con una temperatura por nodo
```

## Condiciones de borde

`cb` tiene dos filas: izquierda y derecha.

```matlab
cb = [tipoIzquierda valorIzquierda dato3Izquierda;
      tipoDerecha   valorDerecha   dato3Derecha];
```

Tipos:

- `1`: Dirichlet. La columna 2 es la temperatura prescrita y la columna 3 vale `-1`.
- `2`: Neumann. La columna 2 es el flujo impuesto y la columna 3 vale `-1`.
- `3`: Robin. La columna 2 es `h` y la columna 3 es la temperatura externa `Tinf`.

En esta implementacion el aporte Neumann es `-q` en el nodo extremo. Si la convencion de signos de la catedra define el flujo con la normal opuesta, basta cambiar el signo en `applyBoundaryConditions`.

## Esquemas

```matlab
et = 0; % Euler explicito
et = 1; % Euler implicito
et = 2; % estacionario
```

En un problema temporal, `T(:,1)` es la condicion inicial y cada columna siguiente es un paso temporal. `T(:,end)` es el estado final calculado.

Ejemplo estacionario con Dirichlet en ambos extremos:

```matlab
xnode = linspace(0,1,11);
model.k = 1;
model.c = 0;
model.rho = 1;
model.cp = 1;
model.G = 0;
cb = [1 100 -1; 1 0 -1];
T = elemFinitos(xnode,model,cb,2);
plot(xnode,T,'o-'); grid on;
```

Ejemplo temporal:

```matlab
model.dt = 0.001;
model.maxit = 1000;
model.tol = 1e-8;
model.T0 = zeros(size(xnode));
T = elemFinitos(xnode,model,cb,1);
graficaEvolucion(xnode,T,model.dt);
```
