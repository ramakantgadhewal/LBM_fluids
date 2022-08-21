clear all;
load Ux20000.txt
Nx = 40;
Ny = 40;

x = 1:Nx;
y = 1:Ny;

[X,Y] = meshgrid(x,y);

figure("Name", "mainC")
contourf(X,Y,Ux20000);

%ux-main0.txt