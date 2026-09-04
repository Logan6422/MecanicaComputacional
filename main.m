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


#Ej 1d)
model.rhoCp = 0;
model.k = 1;
model.c = 1;
model.G = 50;

mesh.xi = 0;
mesh.xf = 1;
mesh.h = 0.1;

cb = [1 10 -1;
      3  0.2 50];
#T = difFinitas(mesh, model, cb, et)
T = difFinitas(mesh, model, cb, 0)
