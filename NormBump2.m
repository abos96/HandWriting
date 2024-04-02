function [z, norm]= NormBump2 (params,x,y)

    % Extract parameters from params
    ux = params(1);
    nx = params(2);
    uy = params(3);
    ny = params(4);
    
     % define bump 
    z = Bump1(x,ux,nx).*Bump1(y,uy,ny);
 
    % Define the limits of integration
    x_min = (ux-nx);
    x_max = (ux+nx);
    y_min = (uy+ny);
    y_max = (uy-ny);
    
    % Compute the integral
    % Define the function handle for Bump1
    Bump1_handlex = @(x) Bump1(x, ux, nx);
    Bump1_handley = @(y) Bump1(y, uy, ny);
    normX = integral(Bump1_handlex, x_min, x_max);
    normY = integral(Bump1_handley, y_min, y_max);

    %compute normalized bump function
    z = z ./ (normX.* normY);
    norm = (normX.* normY);


