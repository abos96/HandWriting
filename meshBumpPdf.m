function [Xrn,Yrn] =  meshBumpPdf(params,H,L,n,cl)

ux = params(1);
nx = params(2);
uy = params(3);
ny = params(4);
teta = params(5);

x2plot = linspace(- L/2, + L/2, n);
y2plot = linspace(- H/2, + H/2, n);

[X, Y] = meshgrid(x2plot, y2plot);

[Xr,Yr] = rotate_array_vector(teta,X,Y,ux,uy,0);

Xrn = reshape(Xr,n,n);
Yrn = reshape(Yr,n,n);


Bump1_handlex = @(xx) Bump1(xx, ux, nx);
Bump1_handley = @(yy) Bump1(yy, uy, ny);
Bnx = integral(Bump1_handlex, ux-nx, ux+nx);
Bny = integral(Bump1_handley, uy-ny, uy+ny);

Z = Bump1(Xrn,ux,nx).*Bump1(Yrn,uy,ny)./(Bnx.*Bny);

% Create a contour plot
contour(X, Y, Z, 'LineWidth', 2, 'LineColor', cl);

