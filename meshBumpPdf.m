function [Xrn,Yrn] =  meshBumpPdf(params,H,L,n,cl)

if size(params,2)>1
    
    mode = 'multi';

else

    mode = 'mono';

end

switch mode
    case 'mono'
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


    case 'multi'

        for i = 1:size(params,2)
            % extract parameters
            ux = params(1,i);
            nx = params(2,i);
            uy = params(3,i);
            ny = params(4,i);
            teta = params(5,i);
            W_k = params(6,i); 
            
            if size(L)==1

                x2plot = linspace(- L/2, + L/2, n);
                y2plot = linspace(- H/2, + H/2, n);

            elseif size(L)==[1,2] & size(H)==[1,2]

                x2plot = linspace( L(1), L(2), n);
                y2plot = linspace( H(1), H(2), n);

            end

            
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
            % Define the colormap
            colormap_name = cl; % Change this to the name of your desired colormap
            cmap = colormap(colormap_name);
            
            % Generate a random index
            num_colors = size(cmap, 1);
            random_index = randi(num_colors);
            
            % Retrieve the random color from the colormap
            random_color = cmap(random_index, :);
            contour(X, Y, Z, 'LineWidth', 2, 'LineColor', random_color);
        end


end
