clear all;

load VNet7000.txt

Nx = 600;
Ny = 320;

x = 1:Nx;
y = 1:Ny;

[X,Y] = meshgrid(x,y);

contourf(X,Y,VNet7000, 10);