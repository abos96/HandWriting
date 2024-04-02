function [Xr, Yr] = rotate_array_vector(teta,x,y,ux,uy,traslate2zero)

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