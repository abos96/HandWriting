function loglikelihood = objectiveFunction1D (params,y)
uy = params(1);
ny = params(2);


Bump1_handley = @(y) Bump1(y, uy, ny);
Bny = integral(Bump1_handley, uy-ny, uy+ny);

% define loglikelihood 
Pyi = Bump1(y,uy,ny)./Bny;
N = length(y);

%Likelihood = -sum(log((Pyi)));
loglikelihood = N.*log(Bny) + sum(1./(ny.^2 - (y-uy).^2)) ;




