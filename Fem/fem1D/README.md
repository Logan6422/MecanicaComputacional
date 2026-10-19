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
- `solveCrankNicolson1D.m`: avanza con Crank-Nicolson.
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

El vector local de fuente se calcula como
`Fe = integral(N' * G(x) dx)` mediante cuadratura de Gauss-Legendre. La
cantidad de puntos se elige en el modelo:

```matlab
model.quadraturePoints = 2;
```

Se admiten de 1 a 4 puntos:

| Puntos de Gauss | Grado maximo integrado exactamente | Uso habitual |
| ---: | ---: | --- |
| 1 | 1 | `G` aproximadamente constante o muy suave |
| 2 | 3 | Fuente lineal/cuadratica y opcion general economica |
| 3 | 5 | Fuente cubica o mas variable |
| 4 | 7 | Fuente polinomica de mayor grado |

En cada elemento se transforma el intervalo local a `[-1,1]`, se evalua `G`
en los puntos y se ponderan sus valores. El grado indicado corresponde a la
integranda completa `N' * G(x)`, no solamente a `G`.

## Relacion con otras cuadraturas

Estas son reglas distintas de integracion, aunque todas aproximan una
integral en un intervalo:

| Metodo | Puntos/nodos | Precision polinomica tipica |
| --- | ---: | ---: |
| Trapecio | 2 nodos extremos | grado 1 |
| Simpson 1/3 | 3 nodos: extremo, medio, extremo | grado 3 |
| Gauss-Legendre | 1 a 4 puntos interiores | grado `2*n-1` |

Por ejemplo, `model.quadraturePoints = 2` significa **Gauss de 2 puntos**,
no trapecio. Gauss usa puntos interiores y pesos optimizados; trapecio y
Simpson usan nodos que incluyen los extremos. En esta implementacion se usa
Gauss-Legendre porque se integra naturalmente con las funciones de forma FEM.

Para una fuente constante, tambien se puede usar simplemente `model.G = 100`.

La cuadratura se ejecuta durante el ensamblaje de `F`, antes de aplicar las
condiciones de borde. Por eso se usa igual en los problemas estacionarios
(`et = 2`) y transitorios (`et = 0`, `et = 1` y `et = 3`): los metodos
temporales reutilizan el mismo vector `F` en cada paso. Si `G` depende del
tiempo, esta implementacion no lo actualiza y habria que ensamblar `F` dentro
del bucle temporal.

Para los esquemas temporales tambien se necesitan:

```matlab
model.dt
model.maxit
model.tol
model.T0   % vector columna o fila con una temperatura por nodo
```

Si el enunciado proporciona un parametro de malla

$$
\lambda=\frac{k\,\Delta t}{\rho c_p h^2},
$$

se puede calcular el paso como

```matlab
h = (xnode(end) - xnode(1))/(length(xnode) - 1);
model.dt = lambda*model.rho*model.cp*h^2/model.k;
```

El template deja esta opcion disponible con `usarLambda = true`. Para FEM
explicito, la estabilidad depende de los autovalores de la matriz
`C\K`, por lo que esta formula sirve como parametro orientativo y debe
usarse con un valor conservador de `lambda`. La tolerancia se aplica al
cambio relativo entre dos soluciones consecutivas; `T(:,1)` es la condicion
inicial y `T(:,end)` es el ultimo estado calculado.

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

Para Robin se usa la condicion convectiva

$$
-k\frac{dT}{dn}=h(T-T_{inf}),
$$

por lo que la columna 2 contiene el coeficiente convectivo `h` completo,
no el cociente `h/k`.

En esta implementacion el aporte Neumann es `-q` en el nodo extremo. Si la convencion de signos de la catedra define el flujo con la normal opuesta, basta cambiar el signo en `applyBoundaryConditions`.

## Esquemas

```matlab
et = 0; % Euler explicito
et = 1; % Euler implicito
et = 2; % estacionario
et = 3; % Crank-Nicolson
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
