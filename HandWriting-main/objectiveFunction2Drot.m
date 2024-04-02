function loglikelihood = objectiveFunction2Drot (params,x,y)
debug = 0;
% extract parameters
ux = params(1);
nx = params(2);
uy = params(3);
ny = params(4);
teta = params(5);


[Bnx, Bny, Xr, Yr, N] = rotateBn(teta,x,y,ux,uy,nx,ny);


loglikelihood = N.*log(Bny) + sum(1./(nx.^2 - (Xr).^2)) + ...
                N.*log(Bnx) + sum(1./(ny.^2 - (Yr).^2)) ;

if debug

% compute initial pdf
[Xi,Yi,Zi]= NormBump2rot ([ux nx uy ny teta],200,200,1000);
% Plot initial bump functio
hold on
contour(Xi,Yi,Zi,'LineWidth', 1 , 'LineColor', 'black')
pause(0.5)

end