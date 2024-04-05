function plot_gaussian_2d(mu, Sigma, varargin)
    % Plot 2D Gaussian distribution
    
    % Define the number of points to generate
    n = 100;
    
    % Generate points on a unit circle
    t = linspace(0, 2*pi, n);
    xy = [cos(t); sin(t)];
    
    % Scale and rotate the points based on the covariance matrix
    d = sqrtm(Sigma);
    xy = d * xy;
    
    % Translate the points based on the mean
    xy = xy + mu';
    
    % Plot the points with specified line width
    plot(xy(1, :), xy(2, :), 'LineWidth', 3, varargin{:});
end