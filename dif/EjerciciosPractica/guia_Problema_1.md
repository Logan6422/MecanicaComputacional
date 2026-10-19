# Guia de estudio - Problema 1

## 1. Ecuacion general del calor

La ecuacion unidimensional del calor, incluyendo difusion, reaccion y fuente, puede escribirse como

$$
\rho c_p \frac{\partial T}{\partial t}
= k \frac{\partial^2 T}{\partial x^2} - cT + G(x,t),
\qquad 0 < x < L,
$$

donde:

- $T(x,t)$ es la temperatura.
- $k$ es la conductividad termica.
- $\rho$ es la densidad.
- $c_p$ es el calor especifico.
- $c$ es el coeficiente del termino de reaccion.
- $G$ es una fuente de calor.

Si dividimos por $\rho c_p$:

$$
\frac{\partial T}{\partial t}
= \frac{k}{\rho c_p}\frac{\partial^2T}{\partial x^2}
- \frac{c}{\rho c_p}T + \frac{G}{\rho c_p}.
$$

En este problema:

$$
\frac{\partial T}{\partial t}
= \frac{\partial^2T}{\partial x^2} + 100,
\qquad 0 < x < 1,\quad t > 0.
$$

Por lo tanto, se trabaja con $L=1$, $k=1$, $\rho c_p=1$, $c=0$ y $G=100$.

## 2. Condiciones del problema

Condicion de borde izquierda, tipo Dirichlet:

$$
T(0,t)=0.
$$

Condicion de borde derecha, tipo Robin:

$$
\frac{\partial T}{\partial x}(1,t)=-T(1,t).
$$

La ultima condicion tambien puede escribirse como

$$
\frac{\partial T}{\partial x}(1,t) + T(1,t)=0.
$$

Condicion inicial:

$$
T(x,0)=10x(1-x)^2.
$$

---

# Punto 1 - Estado estacionario con 4 nodos

## 3. Ecuacion estacionaria

En estado estacionario la temperatura no depende del tiempo, por lo que

$$
\frac{\partial T}{\partial t}=0.
$$

La ecuacion queda

$$
\frac{d^2T}{dx^2}+100=0,
$$
o equivalentemente,
$$
\frac{d^2T}{dx^2}=-100.
$$

Para una malla uniforme de 4 nodos:

$$
N=4,\qquad h=\frac{L}{N-1}=\frac{1}{3}.
$$

Los nodos son

$$
\begin{array}{c|c|c|c|c}
\text{nodo} & 0 & 1 & 2 & 3\\ \hline
x_i & 0 & \frac{1}{3} & \frac{2}{3} & 1\\
T_i & T_0 & T_1 & T_2 & T_3
\end{array}
$$

## 4. Nodos interiores: diferencias finitas

Para un nodo interior se usa la diferencia central:

$$
\frac{d^2T}{dx^2}(x_i)
\approx \frac{T_{i-1}-2T_i+T_{i+1}}{h^2}.
$$

Como $T''=-100$:

$$
\frac{T_{i-1}-2T_i+T_{i+1}}{h^2}=-100.
$$

Reordenando:

$$
-T_{i-1}+2T_i-T_{i+1}=100h^2.
$$

Dado que $h=1/3$:

$$
100h^2=\frac{100}{9}.
$$

Para el nodo $1$:

$$
-T_0+2T_1-T_2=\frac{100}{9}.
$$

Para el nodo $2$:

$$
-T_1+2T_2-T_3=\frac{100}{9}.
$$

## 5. Condiciones de borde discretizadas

### Borde izquierdo

La condicion es directamente

$$
T_0=0.
$$

### Borde derecho Robin

La derivada se aproxima mediante una diferencia central usando un nodo fantasma $T_4$:

$$
\frac{T_4-T_2}{2h}=-T_3.
$$

De aqui:

$$
T_4=T_2-2hT_3.
$$

Aplicamos la ecuacion estacionaria en el nodo $3$:

$$
\frac{T_2-2T_3+T_4}{h^2}=-100.
$$

Sustituyendo $T_4$:

$$
\frac{2T_2-2(1+h)T_3}{h^2}=-100.
$$

Reordenando:

$$
-2T_2+2(1+h)T_3=100h^2.
$$

Como $h=1/3$:

$$
-2T_2+\frac{8}{3}T_3=\frac{100}{9}.
$$

## 6. Sistema lineal completo

Juntando las cuatro ecuaciones:

$$
\begin{bmatrix}
1 & 0 & 0 & 0\\
-1 & 2 & -1 & 0\\
0 & -1 & 2 & -1\\
0 & 0 & -2 & \frac{8}{3}
\end{bmatrix}
\begin{bmatrix}
T_0\\T_1\\T_2\\T_3
\end{bmatrix}
=
\begin{bmatrix}
0\\\frac{100}{9}\\\frac{100}{9}\\\frac{100}{9}
\end{bmatrix}.
$$

En forma abreviada:

$$
A\mathbf{T}=\mathbf{b}.
$$

La solucion numerica se obtiene mediante

```matlab
T = A\\b;
```

La implementacion equivalente para este punto es:

```matlab
x = [0, 1/3, 2/3, 1];
h = 1/3;

A = [ 1,  0,  0,       0;
     -1,  2, -1,       0;
      0, -1,  2,      -1;
      0,  0, -2,  2+2*h ];

b = [0; 100*h^2; 100*h^2; 100*h^2];
T = A\\b;
```

## 7. Solucion exacta y error

Integrando dos veces la ecuacion $T''=-100$ y aplicando las condiciones de borde:

$$
T_{ex}(x)=25x(3-2x).
$$

Para comparar cada nodo se calcula el error relativo:

$$
\varepsilon_i
=\left|\frac{T_i-T_{ex}(x_i)}{T_{ex}(x_i)}\right|.
$$

En el nodo $x=0$, la solucion exacta vale cero, por lo que el error relativo no esta definido. En ese nodo se informa directamente el error absoluto o se fija el error relativo como cero si ambas soluciones valen cero.

---

# Punto 2 - Problema transitorio con 41 nodos

## 8. Malla y parametro temporal

Para 41 nodos:

$$
N=41,\qquad h=\frac{1}{40}.
$$

El parametro pedido es

$$
\lambda=\frac{k\Delta t}{\rho c_p h^2}=40.
$$

Como $k=1$ y $\rho c_p=1$:

$$
\Delta t=40h^2=40\left(\frac{1}{40}\right)^2=\frac{1}{40}=0.025.
$$

## 9. Esquema implícito

Para un nodo interior, el esquema de Euler implicito es

$$
\frac{T_i^{n+1}-T_i^n}{\Delta t}
=\frac{T_{i-1}^{n+1}-2T_i^{n+1}+T_{i+1}^{n+1}}{h^2}+100.
$$

Multiplicando por $\Delta t$ y usando $\lambda=\Delta t/h^2$:

$$
-\lambda T_{i-1}^{n+1}
+(1+2\lambda)T_i^{n+1}
-\lambda T_{i+1}^{n+1}
=T_i^n+100\Delta t.
$$

Con $\lambda=40$:

$$
-40T_{i-1}^{n+1}+81T_i^{n+1}-40T_{i+1}^{n+1}
=T_i^n+100\Delta t.
$$

En cada paso temporal se resuelve nuevamente un sistema lineal:

$$
A_{imp}\mathbf{T}^{n+1}=\mathbf{b}(\mathbf{T}^n).
$$

La condicion inicial se evalua en los nodos:

$$
T_i^0=10x_i(1-x_i)^2.
$$

La condicion de borde izquierda sigue siendo

$$
T_0^{n+1}=0.
$$

Para el borde derecho se usa un nodo fantasma en cada paso. La condicion
Robin en el instante $n+1$ es

$$
\frac{T_N^{n+1}-T_{N-2}^{n+1}}{2h}=-T_{N-1}^{n+1},
$$

por lo tanto

$$
T_N^{n+1}=T_{N-2}^{n+1}-2hT_{N-1}^{n+1}.
$$

Reemplazando este valor en la ecuacion del calor aplicada al ultimo nodo:

$$
\frac{T_{N-1}^{n+1}-T_{N-1}^{n}}{\Delta t}
=\frac{2T_{N-2}^{n+1}-2(1+h)T_{N-1}^{n+1}}{h^2}+100.
$$

Con $\lambda=\Delta t/h^2$, la fila de borde para Euler implicito queda

$$
-2\lambda T_{N-2}^{n+1}
+\left[1+2\lambda(1+h)\right]T_{N-1}^{n+1}
=T_{N-1}^{n}+100\Delta t.
$$

La fila izquierda sigue siendo $T_0^{n+1}=0$. Esta fila de borde es la que
debe completar las filas interiores de $A_{imp}$.

## 10. Esquema de Crank-Nicolson

Crank-Nicolson promedia el termino espacial entre los instantes $n$ y
$n+1$:

$$
\frac{T_i^{n+1}-T_i^n}{\Delta t}
=\frac{1}{2}\left[
\frac{T_{i-1}^{n+1}-2T_i^{n+1}+T_{i+1}^{n+1}}{h^2}
+\frac{T_{i-1}^{n}-2T_i^{n}+T_{i+1}^{n}}{h^2}
\right]+100.
$$

Para un nodo interior, la ecuacion queda

$$
-\frac{\lambda}{2}T_{i-1}^{n+1}
+(1+\lambda)T_i^{n+1}
-\frac{\lambda}{2}T_{i+1}^{n+1}
=\frac{\lambda}{2}T_{i-1}^{n}
+(1-\lambda)T_i^{n}
+\frac{\lambda}{2}T_{i+1}^{n}
+100\Delta t.
$$

Para el ultimo nodo se usa la misma eliminacion del nodo fantasma, pero en
ambos instantes. La fila Robin resulta

$$
-\lambda T_{N-2}^{n+1}
+\left[1+\lambda(1+h)\right]T_{N-1}^{n+1}
=\lambda T_{N-2}^{n}
+\left[1-\lambda(1+h)\right]T_{N-1}^{n}
+100\Delta t.
$$

Por lo tanto, en cada paso se resuelve

$$
A_{CN}\mathbf{T}^{n+1}=B_{CN}\mathbf{T}^{n}+\mathbf{f}_{CN}.
$$

La diferencia principal es que el implícito evalua toda la difusion en el
nuevo instante, mientras que Crank-Nicolson la promedia. Ambos requieren
resolver un sistema lineal en cada paso.

## 11. Criterio de estado estacionario

Se considera que se alcanzo el estado estacionario cuando el cambio entre dos pasos consecutivos es menor que la tolerancia:

$$
\left\|\mathbf{T}^{n+1}-\mathbf{T}^{n}\right\|_\infty<10^{-4}.
$$

El instante correspondiente es

$$
 t_{conv}=n\Delta t.
$$

Ademas, para encontrar el primer instante en que $T(1/2,t)>24$, se busca el primer paso que cumpla

$$
T_{21}^{n}>24,
$$

porque con 41 nodos el punto $x=0.5$ corresponde al nodo $21$ si se empieza a contar desde uno.

## 12. Valores solicitados al final

Los puntos pedidos son

$$
 x=0,\;0.25,\;0.5,\;0.75,\;1.
$$

Con $h=1/40$, corresponden a los indices

$$
1,\;11,\;21,\;31,\;41.
$$

Se informa en cada uno:

1. La temperatura numerica final.
2. La temperatura exacta $T_{ex}(x)=25x(3-2x)$.
3. El error relativo.

La estructura de calculo en Octave es:

```matlab
x = linspace(0, 1, 41);
T0 = 10*x.*(1-x).^2;

lambda = 40;
h = 1/40;
dt = lambda*h^2;

indices = [1, 11, 21, 31, 41];
posiciones = x(indices);
```

El archivo `Ejercicio_a.m` contiene la ejecucion numerica de estos pasos usando la funcion general `difFinitas`.

---

# Condicion de frontera Neumann en el problema estacionario

## 13. Que significa una condicion Neumann

Una condicion Neumann prescribe la derivada normal o el flujo en el borde.
En un extremo derecho podemos escribir

$$
\frac{dT}{dx}(L)=g.
$$

Si se prescribe el flujo fisico $q$ y la ley de Fourier es

$$
q=-k\frac{dT}{dx},
$$

entonces

$$
g=-\frac{q}{k}.
$$

Es importante fijar esta convencion antes de armar el sistema, porque el
signo depende de si se da la derivada $T_x$ o el flujo $q$.

## 14. Neumann con nodo fantasma en el extremo derecho

Para el problema estacionario general

$$
T''+G=0,
$$

la derivada central en el extremo derecho $x=L$ es

$$
\frac{T_N-T_{N-2}}{2h}=g,
$$

donde $T_N$ es el nodo fantasma. Se despeja:

$$
T_N=T_{N-2}+2hg.
$$

Luego se aplica la ecuacion diferencial en el ultimo nodo real:

$$
\frac{T_{N-2}-2T_{N-1}+T_N}{h^2}=-G.
$$

Sustituyendo el fantasma:

$$
\frac{2T_{N-2}-2T_{N-1}+2hg}{h^2}=-G.
$$

La fila que se incorpora al sistema es

$$
-2T_{N-2}+2T_{N-1}=Gh^2+2hg.
$$

Para flujo prescrito $q=-kT_x$, la misma fila se escribe como

$$
-2T_{N-2}+2T_{N-1}=Gh^2-\frac{2hq}{k}.
$$

Por ejemplo, con $G=100$, $k=1$ y flujo nulo ($q=0$, por lo que $g=0$):

$$
-2T_{N-2}+2T_{N-1}=100h^2.
$$

Esta es la fila que reemplaza a la fila Robin del sistema del punto 1.

## 15. Neumann en el extremo izquierdo

Si la condicion esta en $x=0$ y se prescribe $T_x(0)=g_0$, usando un fantasma
$T_{-1}$:

$$
\frac{T_1-T_{-1}}{2h}=g_0,
\qquad T_{-1}=T_1-2hg_0.
$$

Al aplicar la ecuacion diferencial en $x=0$ se obtiene

$$
2T_0-2T_1=Gh^2-2hg_0.
$$

La expresion cambia de posicion, pero el procedimiento es el mismo:
escribir la ecuacion de flujo, despejar el nodo fantasma y reemplazarlo en
la ecuacion diferencial del nodo de borde.

## 16. Robin como caso relacionado

Una condicion Robin combina derivada y temperatura:

$$
T_x(L)+\beta T(L)=r.
$$

Con diferencia central:

$$
\frac{T_N-T_{N-2}}{2h}+\beta T_{N-1}=r,
$$

de donde

$$
T_N=T_{N-2}+2h\left(r-\beta T_{N-1}\right).
$$

Para el problema del examen, $\beta=1$ y $r=0$, entonces

$$
T_N=T_{N-2}-2hT_{N-1},
$$

que es exactamente la expresion utilizada en el punto 1.

## 17. Como identificar $h$ y $T_{inf}$

La condicion Robin de conveccion suele escribirse como

$$
-k\frac{dT}{dx}=h\left(T-T_{inf}\right).
$$

Despejando la derivada:

$$
\frac{dT}{dx}+\frac{h}{k}T=\frac{h}{k}T_{inf}.
$$

En el enunciado se tiene

$$
\frac{dT}{dx}+T=0.
$$

Comparando ambas expresiones:

$$
\frac{h}{k}=1,
\qquad
\frac{h}{k}T_{inf}=0.
$$

Como en este problema $k=1$:

$$
h=1,
\qquad
T_{inf}=0.
$$

Por eso, usando la convencion de la funcion `difFinitas`, la condicion
Robin se carga como

```matlab
cb = [1, 0, -1;    % Dirichlet izquierdo: T = 0
      3, 1,  0];   % Robin derecho: h = 1, Tinf = 0
```

En general, si $k$ no vale $1$, la comparacion identifica primero
$h/k=1$, por lo que se puede tomar $h=k$, mientras que $T_{inf}=0$.
La funcion necesita el coeficiente convectivo $h$ y no solamente el
coeficiente adimensional $h/k$, por eso hay que tener en cuenta el valor de
$k$ del problema.

## 18. Observacion sobre el sistema con Neumann

Si los dos extremos tienen condiciones Neumann, el problema estacionario
puede no tener solucion unica: una constante puede sumarse a $T$ sin cambiar
$T'$. Ademas, debe cumplirse una condicion de compatibilidad entre la fuente
y los flujos. Si un extremo es Dirichlet, como en el examen, la temperatura
queda fijada y el sistema si tiene solucion unica.
