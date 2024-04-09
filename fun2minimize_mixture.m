function loglikelihood = fun2minimize_mixture (params,x,y)

debug = 0;
number_of_bump = size(params,2);

for i = 1 : number_of_bump
    % extract parameters for a single bump
    ux = params(1,i);
    nx = params(2,i);
    uy = params(3,i);
    ny = params(4,i);
    teta = params(5,i);
    W_k = params(6,i); 
    
    % rotate and traslation of dataset
    [Xr, Yr] = rotate_array_vector(teta,x,y,ux,uy,1);

    if debug
    figure(5)
    hold on
    scatter(Xr, Yr, 'b', 'filled');  % 'b' specifies blue color, 'filled' fills the markers
    [a] = meshBumpPdf([0 nx 0 ny 0]',200,200,500,'black');
    title('Optimized initial points')
    legend('Initial Points',sprintf('Optimal Bump Function teta=%.1f',teta))
    pause(0.2)
    end

    % compute Bn
    Bump1_handlex = @(xx) Bump1(xx, 0, nx);
    Bump1_handley = @(yy) Bump1(yy, 0, ny);

    Bnx = integral(Bump1_handlex, -nx, +nx);
    Bny = integral(Bump1_handley, -ny, +ny);

    % compute pdf of i-th bump
    P_i(i,:) = W_k.*Bump1(Xr,0,nx).*Bump1(Yr,0,ny)./(Bnx.*Bny);
   
end

% compute sum of the bumps
mix_pdf = sum(P_i,1);

% compute log of the sum
Log_sum = -log(mix_pdf);

% sum over each single point
loglikelihood = sum(Log_sum);

