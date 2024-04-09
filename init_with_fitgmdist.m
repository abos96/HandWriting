function iparams = init_with_fitgmdist(x_data,y_data,chunks_num,n_fit)

close all

%fit the data with GMM n_fit time
gmm = {};

for i = 1 : n_fit

    gmm{i} = fitgmdist([x_data;y_data]',chunks_num);
    loglike(i) = gmm{i}.NegativeLogLikelihood;
    
end

% find the best gmm fitting
[val,idx] = min(loglike);

% Plot GMM components
% Plot data points
figure;
scatter(x_data, y_data, 'filled' ,'b','DisplayName', 'Data Points');
hold on
for ii = 1:chunks_num
    mu = gmm{idx}.mu(ii, :);  % Mean of the i-th component
    Sigma = gmm{idx}.Sigma(:, :, ii);  % Covariance matrix of the i-th component
    plot_gaussian_2d(mu, Sigma, 'DisplayName', sprintf('Component %d', i));
end
title('GMM initialization')
% compute the mean ux and ux of the bumps
ux = gmm{idx}.mu(:, 1);
uy = gmm{idx}.mu(:, 2);

% compute covariance matrix
min_sigma = gmm{idx}.Sigma;

% Compute the standard deviation and tilt of each gaussian
for i = 1 : chunks_num

covariance_matrix = min_sigma(:,:,i);

nx(i) = sqrt(diag(covariance_matrix(1,1)));
ny(i) = sqrt(diag(covariance_matrix(2,2)));

% Compute the eigenvectors and eigenvalues of the covariance matrix
[eigenvectors, eigenvalues] = eig(covariance_matrix);

% Find the index of the eigenvalue corresponding to the largest magnitude
[~, max_index] = max(abs(diag(eigenvalues)));

% Get the corresponding eigenvector (principal axis)
principal_axis = eigenvectors(:, max_index);

% Compute the tilt angle (angle between principal axis and x-axis)
tilt(i) = atan2d(-principal_axis(2), -principal_axis(1));


end

% build iparams 
for i = 1 : chunks_num

iparams(1,i) = ux(i);
iparams(2,i) = nx(i);
iparams(3,i) = uy(i);
iparams(4,i) = ny(i);
iparams(5,i) = tilt(i);

end