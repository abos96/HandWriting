%% This example fit a Bump Function over 2d x,y coordinates handwriting dataset
% The Bump-Function is defined in 'Bump1.m' as y = exp(-1./(nx.^2-(x-ux).^2)).*(abs(x-ux) <= nx);

% Input of Bump: 
% x: input data with bump pdf
% nx: half length of the bump
% ux: center of the bump
%--------------------------------------------------------------------------
% 2 do:
% 



clear all
close all
clc

case_variable = 'real'; % 'real' or 'generated'

switch case_variable
    case 'real'
        %% plot data
        % Open a dialog box for selecting a file
        [file,pathdir] = uigetfile('*.*', 'Select a file');
        
        % Check if the user canceled the selection
        if isequal(file, 0)
            disp('User canceled the file selection.');
        else
            % Display the selected file name and path
            disp(['Selected file: ', file]);
            disp(['File path: ', pathdir]);
        end

        orient = 'l';  % Variable defining tablet orientation

        % Read Data
        wd = read(wacomdata,fullfile(pathdir,file),orient); % Call the wacomdata class constructor
        
        %define inair and incontact
        incontact = wd.pressure>0;
        x = wd.position(incontact,1);
        y = wd.position(incontact,2);

        %plot trajectory in contact
        marker_size = 2;
        % Plot x-y data with filled blue markers
        plot(x, y, 'bo', 'MarkerFaceColor', 'blue', 'MarkerSize', marker_size);
        xlabel('X');
        ylabel('Y');
        title('Plot of X-Y Data with Blue Markers');

        % paper parameters
        L = [0 250]; 
        H = [-120 30];
        n = 500; % mesh definition for "contour"
        n_fit = 50;
        max_chunks = 20;
        %% Init GMM
        %Initialization of GMM defining best chunk_num

        xplot = [];
        min_aic = [];
        min_bic = [];
        figure
        parfor i = 2:max_chunks
        
            try
                [min_aic(i-1),min_bic(i-1)] = fit_gmm_and_compute_likelihood(x', y', i, n_fit,0);
            catch exception
                disp('There was an error fitting the Gaussian mixture model - try regularization...')
                error = exception.message
                [min_aic(i-1),min_bic(i-1)] = fit_gmm_and_compute_likelihood(x', y', i, n_fit,1);
            end
            xplot = [xplot i];
        end
        
        
        plot(xplot,min_aic,'bo')
        hold on
        plot(xplot,min_aic,'b')
        
        plot(xplot,min_bic,'ro')
        plot(xplot,min_bic,'r')
        hold off
        legend('aic','','bic','')
        title('Optimization GMM')
        xlabel('Number of Gaussian')
        ylabel('Optimiztion Value BIC-AIC')
        
        [val,chunks_num] = min(min_bic);
        best_chunk = chunks_num + 1;
        hold on
        % Plot the big green circle around the best chunks num
        xl = xline(best_chunk,'-.','Best Cluster Number','DisplayName','Best Cluster');

        % Create a dialog box to input the number of clusters
        prompt = {'Enter the desired number of clusters:'};
        dlg_title = 'Number of Clusters';
        num_lines = 1;
        default_answer = {'5'};  % Default value
        answer = inputdlg(prompt, dlg_title, num_lines, default_answer);
        
        % Convert the user input to a numeric value
        best_chunk = str2double(answer{1});
        
        % Display the entered number of clusters
        fprintf('Number of clusters selected: %d\n', best_chunk);

        %% RUN GMM to initialize the bumps
        % add vargin gain to std
        iparams = init_with_fitgmdist(x',y',best_chunk,n_fit);

        % refine step - delete gm with tilt>45
        ToDelete = find(abs(iparams(5,:))>45);
        best_chunk = best_chunk - length(ToDelete);
        iparams(:,ToDelete) = [];
        
        
        % boost std of iparams
        gain_x = 4;
        gain_y = 4;
        iparams(2,:) = gain_x.*iparams(2,:) ;
        iparams(4,:) = gain_y.*iparams(4,:) ;
        
        % define Wk initial and bounds
        W_k_n = (ones(best_chunk,1)./best_chunk)';
        x0_n = [iparams;W_k_n];

        % Plot initial points
        figure(4)
        hold on
        %initial point 2 fit
        scatter(x, y, 'b', 'filled');  % 'b' specifies blue color, 'filled' fills the markers
        % compute and plot mesh with initial parameters
        x0_n(5,:) = -x0_n(5,:);
        [Xi,Yi] = meshBumpPdf(x0_n,H,L,n,'parula');
        title('Initial bumps')
        %legend('Initial Points','Optimal Bump Function','Initial Bump Function')

        %% define mixture parameters
        % replicate x0 for each bump function
        % x0_n = repmat(x0,1,chunks_num);
        % Bounds on  
          %ux   %nx   %uy    %ny  %teta
        lb_linear = [-inf  0    -inf    0    -30]; % ux: Lmin <= ux <= Lmax % nx: Lmin/2 <= nx <= Lmax/2
        ub_linear = [inf   150   inf   50   30];  % uy: Hmin <= uy <= Hmax % ny: Hmin/2 <= ny <= Hmax/2
        lb_n = repmat(lb_linear',1,best_chunk); 
        ub_n = repmat(ub_linear',1,best_chunk); 
        
        % define Wk initial and bounds
        lb_n = [lb_n; zeros(best_chunk,1)'];
        ub_n = [ub_n; ones(best_chunk,1)'];

        %% solving optimization problem
        % Objective function: 
        obj = @(params) fun2minimize_mixture (params,x,y);
        
        % Define the problem structure
        fun = obj;
        
        % Define the nonlinear constraint function
        nonlinear_constraint_function = @(params) nonlinear_constraints_mixture(params,x,y);
        
        
        % Options for the solver (optional)
        options = optimoptions('fmincon', 'Display', 'iter','EnableFeasibilityMode',false,'MaxFunctionEvaluations',10^5);
        
        % Solve the problem
        [optimal_params, ~, exitflag, ~] = fmincon(fun, x0_n, [], [], [], [], lb_n, ub_n, nonlinear_constraint_function, options);
        
        %results
        optimal_ux = optimal_params(1);
        optimal_nx = optimal_params(2);
        optimal_uy = optimal_params(3);
        optimal_ny = optimal_params(4);
        optimal_teta = optimal_params(5);
        
        %% Plot final PDF guessing
        %
        % Plot initial points
        figure(3)
        hold on
        %initial point 2 fit
        scatter(x, y, 'b', 'filled');  % 'b' specifies blue color, 'filled' fills the markers
        % compute and plot mesh with optimal parameters
        [Xo,Yo] = meshBumpPdf(optimal_params,H,L,n,'autumn');
        title('Final Bumps')
        %legend('Initial Points','Optimal Bump Function','Initial Bump Function')
        
        
        % Plot initial points
        figure(4)
        hold on
        %initial point 2 fit
        scatter(x, y, 'b', 'filled');  % 'b' specifies blue color, 'filled' fills the markers
        % compute and plot mesh with initial parameters
        [Xi,Yi] = meshBumpPdf(x0_n,H,L,n,'parula');
        title('Initial bumps')
        %legend('Initial Points','Optimal Bump Function','Initial Bump Function')
        
        
        % compute initial pdf
        % [Xi,Yi] = meshBumpPdf (iparams,H,L,n);
        
        
        
        % % Create a contour plot
        % contour(Xo, Yo, Zo, 'LineWidth', 2, 'LineColor', 'red');
        % colormap('autumn');
        % contour(Xi,Yi,Zi,'LineWidth', 1 , 'LineColor', 'black')
        % 
        % xlabel('x');
        % ylabel('y');
        % zlabel('z');
        % title('Optimization Bump PDF');
        % axis equal
        % % Add legend with text of the initial parameters
        % legend_text_in = sprintf('initial nx = %.1f, initial ny = %.1f, initial teta = %.1f', iparams(2), iparams(4), iparams(5));
        % legend_text_opt = sprintf('optimal nx = %.1f, optimal ny = %.1f, optimal teta = %.1f', optimal_params(2), optimal_params(4), optimal_params(5));
        % legend('Generated Data',legend_text_opt,legend_text_in);
        % 




    case 'generated'

        chunks_num = 7; % number of clusters
        n_fit = 20; % number of gmm fitting to initialize the model
        max_chunks = 12;

        % paper parameters
        L = 200; 
        H = 200;
        n = 500; % mesh definition for "contour"

        % define chunks parameters
        ux = [-5 -40 5 40 -20 20 10];
        nx = [30 35 30 20 15 12 13];
        uy = [-60 0 53 2 20 -20 13];
        ny = [8 10 10 5 8 12 5];
        teta = [7 9 13 -5 3 -6 9];
        N = [100 100 100 100 100 100 100];
        

        % generete points in [ux,uy] with [nx,ny] (vectorized)
        [x,y] =  GenerateRectangle_points (ux,nx,uy,ny,N);

       
        % rotate points (vectorized)
        [rotated_x, rotated_y] = rotate_array_vector(teta,x,y,ux,uy,0);
       

        % plot
        figure(1)
        hold on
        scatter(x,y,'blue','filled')
        scatter(rotated_x,rotated_y,'filled','o','r')
        legend('original',sprintf('rotated teta = %.1f',teta))
        title('point generation')

%% Initialization of GMM defining best chunk_num

xplot = [];
min_aic = [];
min_bic = [];
figure
for i = 2:max_chunks

    try
        [min_aic(i-1),min_bic(i-1)] = fit_gmm_and_compute_likelihood(rotated_x, rotated_y, i, n_fit,0);
    catch exception
        disp('There was an error fitting the Gaussian mixture model - try regularization...')
        error = exception.message
        [min_aic(i-1),min_bic(i-1)] = fit_gmm_and_compute_likelihood(rotated_x, rotated_y, i, n_fit,1);
    end
    xplot = [xplot i];
end


plot(xplot,min_aic,'bo')
hold on
plot(xplot,min_aic,'b')

plot(xplot,min_bic,'ro')
plot(xplot,min_bic,'r')
hold off
legend('aic','','bic','')
title('Optimization GMM')
xlabel('Number of Gaussian')
ylabel('Optimiztion Value BIC-AIC')

[val,chunks_num] = min(min_bic);
best_chunk = chunks_num + 1;
hold on
% Plot the big green circle around the best chunks num
xl = xline(best_chunk,'-.','Best Cluster Number','DisplayName','Best Cluster');

%% Initialization of GMM 

% add vargin gain to std
iparams = init_with_fitgmdist(rotated_x,rotated_y,chunks_num,n_fit);

% boost std of iparams
gain = 3;
iparams(2,:) = gain.*iparams(2,:) ;
iparams(4,:) = gain.*iparams(4,:) ;


        


% %% Define Optimization problem
% % Initial guess for the variables ux,nx,uy,ny
% 
% iparams(1) = 0; % ux
% iparams(2) = 30; % nx
% iparams(3) = 50; % uy
% iparams(4) = 10; % ny
% iparams(5) = 0;
% x0 = iparams;
% 
% % Bounds on  
%       %ux   %nx   %uy    %ny  %teta
% lb = [-inf  0    -inf    0    -60]; % ux: Lmin <= ux <= Lmax % nx: Lmin/2 <= nx <= Lmax/2
% ub = [inf   100   inf   100   60];  % uy: Hmin <= uy <= Hmax % ny: Hmin/2 <= ny <= Hmax/2
% 
% % inequality const for manual initialization
% A = [-1 -1  0  0 0;... %  ux+nx >= max(x)
%       1 -1  0  0 0;... %  ux-nx <= min(x)
%       0 -1  0  0 0;... %  nx >= 1
%       0  0 -1 -1 0;... %  uy-ny >= max(y)
%       0  0  1 -1 0;... %  ux-nx <= min(y)
%       0  0  0 -1 0];   %  ny >= 1
% 
% b = [-max(rotated_x)+0.001 min(rotated_x)-0.001  -1,...
%      -max(rotated_y)+0.001 min(rotated_y)-0.001  -1];
% 
% 
% 
% 
% %% Initialization Optimization Problem
% % check if initial point satisfy constrain
% % Inizialization
% if (any(((A*iparams') < b') == 0))
% 
%     f = zeros(size(x0)); % assumes x0 is the initial point
%     % Bounds on  
%       %ux   %nx   %uy    %ny  %teta
%     lb_linear = [-inf  0    -inf    0    -20]; % ux: Lmin <= ux <= Lmax % nx: Lmin/2 <= nx <= Lmax/2
%     ub_linear = [inf   100   inf   100   20];  % uy: Hmin <= uy <= Hmax % ny: Hmin/2 <= ny <= Hmax/2
% 
%     % need in Matlab R2024a (or Mac)
%     options = optimoptions('linprog','Algorithm','interior-point');
%     [xnew,fval,exitflag,output]= linprog(f,A,b,[],[], lb_linear ,ub_linear,options);  % Solve the linear programming problem to see if there is a feasible point
%     x0 = xnew;
%     % compute and plot mesh with initial parameters
%     figure(2)
%     hold on
%     scatter(rotated_x, rotated_y, 'b', 'filled');  % 'b' specifies blue color, 'filled' fills the markers
%     [a] = meshBumpPdf(x0,H,L,n,'black');
%     title('Optimized initial points')
%     legend('Initial Points','Optimal initial Bump Function')
% 
% else 
%     x0 = iparams;
%     figure(2)
%     hold on
%     scatter(rotated_x, rotated_y, 'b', 'filled');  % 'b' specifies blue color, 'filled' fills the markers
%     [a] = meshBumpPdf(x0,H,L,n,'black');
%     title('Optimized initial points')
%     legend('Initial Points','Optimal initial Bump Function')
% end

%% define mixture parameters
% replicate x0 for each bump function
% x0_n = repmat(x0,1,chunks_num);
% Bounds on  
  %ux   %nx   %uy    %ny  %teta
lb_linear = [-inf  0    -inf    0    -20]; % ux: Lmin <= ux <= Lmax % nx: Lmin/2 <= nx <= Lmax/2
ub_linear = [inf   100   inf   100   20];  % uy: Hmin <= uy <= Hmax % ny: Hmin/2 <= ny <= Hmax/2
lb_n = repmat(lb_linear',1,chunks_num); 
ub_n = repmat(ub_linear',1,chunks_num); 

% define Wk initial and bounds
W_k_n = (ones(chunks_num,1)./chunks_num)';
x0_n = [iparams;W_k_n];
lb_n = [lb_n; zeros(chunks_num,1)'];
ub_n = [ub_n; ones(chunks_num,1)'];

% Plot initial points
figure(4)
hold on
%initial point 2 fit
scatter(rotated_x, rotated_y, 'b', 'filled');  % 'b' specifies blue color, 'filled' fills the markers
% compute and plot mesh with initial parameters
x0_n(5,:) = -x0_n(5,:);
[Xi,Yi] = meshBumpPdf(x0_n,H,L,n,'parula');
title('Initial bumps')
%legend('Initial Points','Optimal Bump Function','Initial Bump Function')

%% solving optimization problem
% Objective function: 
obj = @(params) fun2minimize_mixture (params,rotated_x,rotated_y);

% Define the problem structure
fun = obj;

% Define the nonlinear constraint function
nonlinear_constraint_function = @(params) nonlinear_constraints_mixture(params,rotated_x,rotated_y);


% Options for the solver (optional)
options = optimoptions('fmincon', 'Display', 'iter','EnableFeasibilityMode',false,'MaxFunctionEvaluations',10^5);

% Solve the problem
[optimal_params, ~, exitflag, ~] = fmincon(fun, x0_n, [], [], [], [], lb_n, ub_n, nonlinear_constraint_function, options);

%results
optimal_ux = optimal_params(1);
optimal_nx = optimal_params(2);
optimal_uy = optimal_params(3);
optimal_ny = optimal_params(4);
optimal_teta = optimal_params(5);

%% Plot final PDF guessing
%
% Plot initial points
figure(3)
hold on
%initial point 2 fit
scatter(rotated_x, rotated_y, 'b', 'filled');  % 'b' specifies blue color, 'filled' fills the markers
% compute and plot mesh with optimal parameters
[Xo,Yo] = meshBumpPdf(optimal_params,H,L,n,'autumn');
title('Final Bumps')
%legend('Initial Points','Optimal Bump Function','Initial Bump Function')


% Plot initial points
figure(4)
hold on
%initial point 2 fit
scatter(rotated_x, rotated_y, 'b', 'filled');  % 'b' specifies blue color, 'filled' fills the markers
% compute and plot mesh with initial parameters
[Xi,Yi] = meshBumpPdf(x0_n,H,L,n,'parula');
title('Initial bumps')
%legend('Initial Points','Optimal Bump Function','Initial Bump Function')


% compute initial pdf
% [Xi,Yi] = meshBumpPdf (iparams,H,L,n);



% % Create a contour plot
% contour(Xo, Yo, Zo, 'LineWidth', 2, 'LineColor', 'red');
% colormap('autumn');
% contour(Xi,Yi,Zi,'LineWidth', 1 , 'LineColor', 'black')
% 
% xlabel('x');
% ylabel('y');
% zlabel('z');
% title('Optimization Bump PDF');
% axis equal
% % Add legend with text of the initial parameters
% legend_text_in = sprintf('initial nx = %.1f, initial ny = %.1f, initial teta = %.1f', iparams(2), iparams(4), iparams(5));
% legend_text_opt = sprintf('optimal nx = %.1f, optimal ny = %.1f, optimal teta = %.1f', optimal_params(2), optimal_params(4), optimal_params(5));
% legend('Generated Data',legend_text_opt,legend_text_in);
% 
end