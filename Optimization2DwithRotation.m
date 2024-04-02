%% This example fit a Bump Function over 2d x,y coordinates handwriting dataset
% The Bump-Function is defined in 'Bump1.m' as y = exp(-1./(nx.^2-(x-ux).^2)).*(abs(x-ux) <= nx);

% Input of Bump: 
% x: input data with bump pdf
% nx: half length of the bump
% ux: center of the bump
%--------------------------------------------------------------------------
%
clear all
close all
clc

case_variable = 'generated'; % 'real' or 'generated'

switch case_variable
    case 'real'
       
    case 'generated'
        % define parameters
        L = 200;
        H = 200;
        N = 300;
        ux = 20;
        nx = 30;
        uy = 0;
        ny = 10;
        n = 500;
        teta = 15;

        %vgenerete points in [ux,uy] with [nx,ny]
        [x,y] =  GenerateRectangle_points (ux,nx,uy,ny,N);
        
        % rotate counter clock-wise around the center [ux,uy]
        [rotated_x, rotated_y] = rotate_array_vector(teta,x,y,ux,uy,0);
       
        % plot
        figure(1)
        hold on
        scatter(x,y,'filled','o','b')
        scatter(rotated_x,rotated_y,'filled','o','r')
        legend('original',sprintf('rotated teta = %.1f',teta))
        title('point generation')
end

%% Define Optimization problem
% Initial guess for the variables ux,nx,uy,ny

iparams(1) = 0; % ux
iparams(2) = 30; % nx
iparams(3) = 50; % uy
iparams(4) = 10; % ny
iparams(5) = 0;
x0 = iparams;

% Bounds on  
      %ux   %nx   %uy    %ny  %teta
lb = [-inf  0    -inf    0    -60]; % ux: Lmin <= ux <= Lmax % nx: Lmin/2 <= nx <= Lmax/2
ub = [inf   100   inf   100   60];  % uy: Hmin <= uy <= Hmax % ny: Hmin/2 <= ny <= Hmax/2

% inequality const
A = [-1 -1  0  0 0;... %  ux+nx >= max(x)
      1 -1  0  0 0;... %  ux-nx <= min(x)
      0 -1  0  0 0;... %  nx >= 1
      0  0 -1 -1 0;... %  uy-ny >= max(y)
      0  0  1 -1 0;... %  ux-nx <= min(y)
      0  0  0 -1 0];   %  ny >= 1

b = [-max(rotated_x)+0.001 min(rotated_x)-0.001  -1,...
     -max(rotated_y)+0.001 min(rotated_y)-0.001  -1];

% Define the nonlinear constraint function
nonlinear_constraint_function = @(params) nonlinear_constraints(params,rotated_x,rotated_y);

%% solving Optimization Problem
% check if initial point satisfy constrain
% Inizialization
if (any(((A*iparams') < b') == 0))

    f = zeros(size(x0)); % assumes x0 is the initial point
    % Bounds on  
      %ux   %nx   %uy    %ny  %teta
    lb_linear = [-inf  0    -inf    0    -20]; % ux: Lmin <= ux <= Lmax % nx: Lmin/2 <= nx <= Lmax/2
    ub_linear = [inf   100   inf   100   20];  % uy: Hmin <= uy <= Hmax % ny: Hmin/2 <= ny <= Hmax/2


    xnew = linprog(f,A,b,[],[], lb_linear ,ub_linear);  % Solve the linear programming problem to see if there is a feasible point
    x0 = xnew;
    % compute and plot mesh with initial parameters
    figure(2)
    hold on
    scatter(rotated_x, rotated_y, 'b', 'filled');  % 'b' specifies blue color, 'filled' fills the markers
    [a] = meshBumpPdf(x0,H,L,n,'black');
    title('Optimized initial points')
    legend('Initial Points','Optimal initial Bump Function')
    
else 
    x0 = iparams;
    figure(2)
    hold on
    scatter(rotated_x, rotated_y, 'b', 'filled');  % 'b' specifies blue color, 'filled' fills the markers
    [a] = meshBumpPdf(x0,H,L,n,'black');
    title('Optimized initial points')
    legend('Initial Points','Optimal initial Bump Function')
end

% Objective function: 
obj = @(params) fun2minimize (params,rotated_x,rotated_y);

% Define the problem structure
fun = obj;

% Options for the solver (optional)
options = optimoptions('fmincon', 'Display', 'iter','EnableFeasibilityMode',false);

% Solve the problem
[optimal_params, ~, exitflag, ~] = fmincon(fun, x0, [], [], [], [], lb, ub, [], options);

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
[Xo,Yo] = meshBumpPdf(optimal_params,H,L,n,'red');

% compute and plot mesh with initial parameters
[Xi,Yi] = meshBumpPdf(iparams,H,L,n,'black');

title('Optimization Bump Function')
legend('Initial Points','Optimal Bump Function','Initial Bump Function')


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
