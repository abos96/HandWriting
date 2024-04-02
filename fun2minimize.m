function loglikelihood = fun2minimize (params,x,y)

debug = 0;

% extract parameters
ux = params(1);
nx = params(2);
uy = params(3);
ny = params(4);
teta = params(5);

% rotate and traslation
[Xr, Yr] = rotate_array_vector(teta,x,y,ux,uy,1);

if debug
    figure(5)
    hold on
    scatter(Xr, Yr, 'b', 'filled');  % 'b' specifies blue color, 'filled' fills the markers
    [a] = meshBumpPdf([0 nx 0 ny 0],200,200,500,'black');
    title('Optimized initial points')
    legend('Initial Points',sprintf('Optimal Bump Function teta=%.1f',teta))
    pause(0.2)
end

% compute bn
Bump1_handlex = @(xx) Bump1(xx, 0, nx);
Bump1_handley = @(yy) Bump1(yy, 0, ny);
% figure(6)
% plot(-nx:0.1:nx,Bump1(-nx:0.1:nx,0,nx))
Bnx = integral(Bump1_handlex, -nx, +nx);
Bny = integral(Bump1_handley, -ny, +ny);

N = length(y);

% compute loglikelihood
loglikelihood = N.*log(Bny) + sum(1./(nx.^2 - (Xr).^2)) - sum(log(double((abs(Xr) <= nx)))) +...
                N.*log(Bnx) + sum(1./(ny.^2 - (Yr).^2)) - sum(log(double((abs(Yr) <= ny))));
% pause(1)
% close all