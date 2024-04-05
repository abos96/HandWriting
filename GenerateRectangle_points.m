function [x,y] =  GenerateRectangle_points (ux,nx,uy,ny,n)


%check length of ux
x = [];
y = [];
number_of_chunks = length(ux); 

for i = 1 : number_of_chunks


    Bnx = integral(@(x) Bump1(x, ux(i),nx(i)),ux(i) - nx(i) , ux(i) + nx(i));
    Bny = integral(@(y) Bump1(y, uy(i),ny(i)),uy(i) - ny(i) , uy(i) + ny(i));
    
    % Define the range for generating random points
    xmin = ux(i)-nx(i);
    xmax = ux(i)+nx(i);
    ymin = uy(i)-ny(i);
    ymax = uy(i)+ny(i);
    
    zmax = max( ...
           (Bump1(linspace(xmin, xmax, n(i)),ux(i),nx(i))./Bnx)...
           .*...
           (Bump1(linspace(ymin, ymax, n(i)),uy(i),ny(i))./Bny)); % Upper bound for y
    
    % Initialize variables to store accepted points
    accepted_points = zeros(2, 0)';
    
    % Generate random points and accept/reject based on the PDF
    while length(accepted_points) < n(i)
        % Generate random point coordinates
        x_rand = (xmax - xmin) * rand() + xmin;
        y_rand = (ymax - ymin) * rand() + ymin;
        
        z_rand = ymax * xmax * rand();
        % Evaluate the PDF at the random point
        pdf_value = (Bump1(x_rand,ux(i),nx(i))./Bnx).*(Bump1(y_rand,uy(i),ny(i))./Bny);
        if z_rand < pdf_value
            % Accept the point
             accepted_points(end + 1,:) = [x_rand; y_rand];
        end
       
       
    end
    % rotate generated points
    % plot generated data
    x = [x accepted_points(:,1)];
    y = [y accepted_points(:,2)];


end