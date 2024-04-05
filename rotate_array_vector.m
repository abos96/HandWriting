function [Xr, Yr] = rotate_array_vector(teta,x,y,ux,uy,traslate2zero)

[rows, columns] = size(x);
if rows==columns

    dimension = 'mono';

elseif rows==1 || columns == 1

    dimension = 'mono';

elseif columns>1 && columns ~= rows
    
    dimension = 'multi';

end


switch dimension
    case 'mono'
        % Compute the rotation matrix
        R = [cosd(teta), -sind(teta); sind(teta), cosd(teta)];
        
        
        if traslate2zero
        
            % Apply the rotation to each point in the grid
            coord = R*[(x(:)-ux)';(y(:)-uy)'];
            Xr = coord(1,:);
            Yr = coord(2,:);
        
        else
        
            % Apply the rotation to each point in the grid
            coord = R*[(x(:)-ux)';(y(:)-uy)'];
            Xr = coord(1,:)+ux;
            Yr = coord(2,:)+uy;
        
        end

    case 'multi'
        
        Xr = [];
        Yr = [];
        for i = 1 : columns

            % Compute the rotation matrix
            R = [cosd(teta(i)), -sind(teta(i)); sind(teta(i)), cosd(teta(i))];
            
            if traslate2zero
            
                % Apply the rotation to each point in the grid
                coord = R*[(x(:,i)-ux(i))';(y(:,i)-uy(i))'];
                Xr = [Xr coord(1,:)];
                Yr = [Yr coord(2,:)];
            
            else
            
                % Apply the rotation to each point in the grid
                coord = R*[(x(:,i)-ux(i))';(y(:,i)-uy(i))'];
                Xr = [Xr coord(1,:)+ux(i)];
                Yr = [Yr coord(2,:)+uy(i)];
            
            end

        end

end