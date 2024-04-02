% debug loglikelihood
clear all
close all
%% parameters
L = 200;
H = 200;
N = 1000;
ux = 10;
nx = 30;
uy = 0;
ny = 5;
n = 100;
teta = 20;

%% generate points with mean=ux,uy and n = nx,ny
[x,y] = GenerateRectangle_points (ux,nx,uy,ny,n);
[Xr, Yr] = rotate_array_vector(teta,x,y,ux,uy);

figure
hold on
scatter(x,y,'filled','o','b')
scatter(Xr,Yr,'filled','o','r')
legend('original',sprintf('rotated teta = %.1f',teta))

