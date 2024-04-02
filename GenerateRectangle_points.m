function [x,y] =  GenerateRectangle_points (ux,nx,uy,ny,n)

Bnx = integral(@(x) Bump1(x, ux,nx),ux-nx,ux+nx);
Bny = integral(@(y) Bump1(y, uy,ny),uy-ny,uy+ny);

% Define the range for generating random points
xmin = ux-nx;
xmax = ux+nx;
ymin = uy-ny;
ymax = uy+ny;

zmax = max( ...
       (Bump1(linspace(xmin, xmax, n),ux,nx)./Bnx)...
       .*...
       (Bump1(linspace(ymin, ymax, n),uy,ny)./Bny)); % Upper bound for y

% Initialize variables to store accepted points
accepted_points = zeros(2, 0)';

% Generate random points and accept/reject based on the PDF
while length(accepted_points) < n
    % Generate random point coordinates
    x_rand = (xmax - xmin) * rand() + xmin;
    y_rand = (ymax - ymin) * rand() + ymin;
    
    z_rand = ymax * xmax * rand();
    % Evaluate the PDF at the random point
    pdf_value = (Bump1(x_rand,ux,nx)./Bnx).*(Bump1(y_rand,uy,ny)./Bny);
    if z_rand < pdf_value
        % Accept the point
         accepted_points(end + 1,:) = [x_rand; y_rand];
    end
   
   
end
% rotate generated points
% plot generated data
x = accepted_points(:,1);
y = accepted_points(:,2);


