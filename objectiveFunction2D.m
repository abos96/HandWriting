function loglikelihood = objectiveFunction2D (params,x,y)
ux = params(1);
nx = params(2);
uy = params(3);
ny = params(4);

Bump1_handlex = @(x) Bump1(x, ux, nx);
Bump1_handley = @(y) Bump1(y, uy, ny);
Bnx = integral(Bump1_handlex, ux-nx, ux+nx);
Bny = integral(Bump1_handley, uy-ny, uy+ny);

N = length(y);

loglikelihood = N.*log(Bny) + sum(1./(ny.^2 - (y-uy).^2)) + ...
                N.*log(Bnx) + sum(1./(nx.^2 - (x-ux).^2)) ;





