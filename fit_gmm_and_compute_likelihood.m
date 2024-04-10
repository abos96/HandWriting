function [min_aic,min_bic] = fit_gmm_and_compute_likelihood(x_data, y_data, chunks_num,n_fit,reg)
    %fit the data with GMM n_fit time
gmm = {};

for i = 1 : n_fit
    if reg
        gmm{i} = fitgmdist([x_data;y_data]',round(chunks_num),'RegularizationValue',0.1);
    else
        gmm{i} = fitgmdist([x_data;y_data]',round(chunks_num));
    end

    aic(i) = gmm{i}.AIC;
    bic(i) = gmm{i}.BIC;
    
end

% find the best gmm fitting
min_aic = min(aic);
min_bic = min(bic);


end