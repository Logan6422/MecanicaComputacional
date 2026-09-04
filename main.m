clear all;clc;
###Main para probar el codigo de diferencias finitas practica 1 mecCo
##
###Ej 1a)
##model.rhoCp = 0;
##model.k = 2;
##model.c = 0;
##model.G = 100;
##
##mesh.xi = 0;
##mesh.xf = 1;
##mesh.h = 0.1;
##
##cb = [1 10 -1;
##      1 50 -1];
###T = difFinitas(mesh, model, cb, et)
##T = difFinitas(mesh, model, cb, 0)

###Ej 1b)
##model.rhoCp = 0;
##model.k = 1;
##model.c = 1;
##model.G = 0;
##
##mesh.xi = 0;
##mesh.xf = 2;
##mesh.h = 0.1;
##
##cb = [1 100 -1;
##      2  0 -1];
###T = difFinitas(mesh, model, cb, et)
##T = difFinitas(mesh, model, cb, 0)

###Ej 1c)
##model.rhoCp = 0;
##model.k = 1;
##model.c = 0;
##model.G = @(x)(x - 3)^2;
##model.transitorio = 0;
##
##mesh.xi = 1;
##mesh.xf = 5;
##mesh.h = 0.1;
##
##cb = [2 2 -1;
##      1  0 -1];
##T = difFinitas(mesh, model, cb)


#Ej 1d)
##model.rhoCp = 0;
##model.k = 1;
##model.c = 1;
##model.G = 50;
##
##mesh.xi = 0;
##mesh.xf = 1;
##mesh.h = 0.1;
##
##cb = [1 10 -1;
##      3  0.2 50];
###T = difFinitas(mesh, model, cb, et)
##T = difFinitas(mesh, model, cb, 0)

###Ej 1e)
##model.rhoCp = 1;
##model.k = 2;
##model.c = 0;
##model.G = @(x) x.^3;
##model.transitorio = 1;
##model.tita = 0;
##model.dt = 0.1;
##
##mesh.xi = 5;
##mesh.xf = 10;
##mesh.h = 0.1;
##
###temperaturas iniciales
##temp0 = 0;
##xnode = mesh.xi:mesh.h:mesh.xf;
##model.T0 = temp0 * ones(length(xnode),1);
###//////////////////////
##model.tf = 10;
##
##
##cb = [3 2 100;
##      1 50 -1];
##T = difFinitas(mesh, model, cb)


###Ej 1f)
##model.rhoCp = 2;
##model.k = 2;
##model.c = 2;
##model.G = 75;
##model.transitorio = 1;
##model.tita = 0;
##model.dt = 0.1;
##
##mesh.xi = 0;
##mesh.xf = 1;
##mesh.h = 0.1;
##
###temperaturas iniciales
##temp0 = 0;
##xnode = mesh.xi:mesh.h:mesh.xf;
##model.T0 = temp0 * ones(length(xnode),1);
###//////////////////////
##model.tf = 10;
##
##
##cb = [1 0 -1;
##      3 2 10];
##T = difFinitas(mesh, model, cb)


#Ej 1g)
model.rhoCp = 1;
model.k = 2;
model.c = -2;
model.G = 0;
model.transitorio = 1;
model.tita = 0;
model.dt = 0.1;

mesh.xi = 0;
mesh.xf = 1;
mesh.h = 0.1;

#temperaturas iniciales
temp0 = 0;
xnode = mesh.xi:mesh.h:mesh.xf;
model.T0 = temp0 * ones(length(xnode),1);
#//////////////////////
model.tf = 10;


cb = [1 50 -1;
      2 50 -1];
T = difFinitas(mesh, model, cb)

