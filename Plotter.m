clear all;
load Ux940.txt
Nx = 70;
Ny = 40;

x = 1:Nx;
y = 1:Ny;

[X,Y] = meshgrid(x,y);

figure("Name", "mainC")
contourf(X,Y,Ux940, 16);

%ux-main0.txt