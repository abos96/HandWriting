function [c, ceq] = nonlinear_constraints(params,x,y)
  
    teta = params(5);
    % Compute other nonlinear constraints (if any)
    % inequality const
    A = [-1 -1  0  0 0;... %  ux+nx >= max(x)
          1 -1  0  0 0;... %  ux-nx <= min(x)
          0 -1  0  0 0;... %  nx >= 1
          0  0 -1 -1 0;... %  uy-ny >= max(y)
          0  0  1 -1 0;... %  ux-nx <= min(y)
          0  0  0 -1 0];   %  ny >= 1
    
    b = [-max(x)+0.001 min(x)-0.001  -1,...
         -max(y)+0.001 min(y)-0.001  -1];

    const = A*params - b' ;
%     % Use PCA to compute the principal components
%     [coeff, score, latent] = pca([x; y]');
%     
%     % The first principal component (coeff(:,1)) represents the direction of maximum variance
%     tilt_vector = coeff(:,1);  % Tilt vector
%     
%     % Compute the tilt angle from the tilt vector
%     tilt_angle = atan2d(tilt_vector(2), tilt_vector(1));  % Angle in degrees
%     
%     % Return the modified parameters and other constraints
%     c = abs(-teta-tilt_angle)-1;
    c = const;
    ceq = [];
end
