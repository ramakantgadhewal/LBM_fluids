clear all;

%{
fileinfo  = dir('*.txt');
fnames = {fileinfo.name};
numFiles = numel(fnames);
%}

load Cd0000.txt

figure()

NSTEPS = 1:1:length(Cd0000);

plot(NSTEPS, Cd0000); 

files = dir('*.txt');
fnames = {files.name};

Nx = 200;
Ny = 120;

x = 1:Nx;
y = 1:Ny;

[X,Y] = meshgrid(x,y);

figure("Name", "mainC")

for i = 1:1:length(files)
    
    V = load(files(i).name, '-ascii');
    contourf(X, Y, V, 12);
    pause(0.01);

    
    % Clears frame unless its the last iteration
    if i~=length(files)
        clf
    end 
    
     
end



