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

L = 200;
H = 200;
case_variable = 'generated'; % 'real' or 'generated'

switch case_variable
    case 'real'
        [filepath] = fullfile('C:\Users\aboschi\OneDrive - Fondazione Istituto Italiano Tecnologia\Documents\MATLAB\Handwriting\Prove\Pazienti\pz3_racconto.csv');  % Select data file
        orient = 'l';  % Variable defining tablet orientation
        % Read Data
        wd = read(wacomdata,filepath,orient); % Call the wacomdata class constructor
        test = wd.position;
        test(1300:end,:) = []; % Select 1 chunk from the text (first row)
        %data2find
        x = test(:,1);  % Data
        y = test(:,2);
        
        % plot
        figure;
        histogram2(Data(:,1),Data(:,2),'Normalization', 'pdf', 'FaceColor', 'r', 'EdgeColor', 'none');
        xlabel('x');
        ylabel('y');
        title('Chunk PDF');
        axis equal
        hold off;

    case 'generated'
        % define parameters
        ux = 10;
        nx = 30;
        uy = 0;
        ny = 5;
        n = 1000;
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
                        accepted_points(end + 1,:) = [x_rand; y_rand];
           
        end

        % plot generated data
        x = accepted_points(:,1);
        y = accepted_points(:,2);
        scatter(x, y, 'b', 'filled');  % 'b' specifies blue color, 'filled' fills the markers
        xlabel('X');
        ylabel('Y');
        title('2D Cloud Points');
        % Add a margin to the plot
        margin = 10;  % Margin size (adjust as needed)
        xlim([-100, 100]);  % Set x-axis limits with margin
        ylim([-100, 100]);  % Set y-axis limits with margin

       
end

%% Define Optimization problem
% Initial guess for the variables ux,nx,uy,ny

iparams(1) = mean(x)+10; % ux
iparams(2) = std(x); % nx
iparams(3) = mean(y); % ux
iparams(4) = std(y); % nx
x0 = iparams;

% Bounds on  
      %ux   %nx   %uy    %ny
lb = [-100  0    -100    0]; % ux: Lmin <= ux <= Lmax % nx: Lmin/2 <= nx <= Lmax/2
ub = [100   50    100    50];  % uy: Hmin <= uy <= Hmax % ny: Hmin/2 <= ny <= Hmax/2

% inequality const
A = [-1 -1  0  0;... %  ux+nx >= max(x)
      1 -1  0  0;... %  ux-nx <= min(x)
      0 -1  0  0;... %  nx >= 1
      0  0 -1 -1;... %  uy-ny >= max(y)
      0  0  1 -1;... %  ux-nx <= min(y)
      0  0  0 -1];   %  ny >= 1

b = [-max(x)+0.001 min(x)+0.001  -1,...
     -max(y)+0.001 min(y)+0.001  -1];

% check if initial point satisfy constrain
% Inizialization
% if (any(((A*iparams') < b') == 0))
% 
%     f = zeros(size(x0)); % assumes x0 is the initial point
% 
%     xnew = linprog(f,A,b,[],[],lb,ub);  % Solve the linear programming problem to see if there is a feasible point
%     x0 = xnew;
% else 
%     x0 = iparams;
% end

 % plot initial pdf
%define domain
x2plot = linspace(- 100, + 100, 1000);
y2plot = linspace(- 100, + 100, 1000);

[X, Y] = meshgrid(x2plot, y2plot);

% Compute bump function
Z = NormBump2(iparams,X,Y);

% Plot initial bump function
figure

hold on
scatter(x, y, 'b', 'filled');  % 'b' specifies blue color, 'filled' fills the markers

% Create a contour plot
contour(X, Y, Z, 'LineWidth', 2, 'LineColor', 'r');

% Set the colormap to red-based
colormap('autumn');
xlabel('x');
ylabel('y');
zlabel('z');
title('Initial Condition');

% Add legend with text of the initial parameters
legend_text = sprintf('ux=%d, nx=%d, uy=%d, ny=%d', iparams(1), iparams(2), iparams(3), iparams(4));
legend('Generated Hand writing Data',legend_text);

%% solving Optimization Problem

% Objective function: f(x) = x^2
obj = @(params) objectiveFunction2D(params,x,y);

% Define the problem structure
fun = obj;

% Options for the solver (optional)
options = optimoptions('fmincon', 'Display', 'iter','EnableFeasibilityMode',true);

% Solve the problem
[optimal_params] = fmincon(fun, x0, A, b, [], [], lb, ub, [], options);

%results
optimal_ux = optimal_params(1);
optimal_nx = optimal_params(2);
optimal_uy = optimal_params(3);
optimal_ny = optimal_params(4);

%% Plot final PDF guessing
 % plot corrispondent pdf
%define domain
x2plot = linspace(-100,100,1000);
y2plot = linspace(-100,100,1000);

[X, Y] = meshgrid(x2plot, y2plot);

% Compute bump function
Z = NormBump2([optimal_ux optimal_nx optimal_uy optimal_ny],X,Y);

% Plot initial bump function
figure

% Create a contour plot
contour(X, Y, Z, 'LineWidth', 2);

% Set the colormap to red-based
colormap('autumn');
xlabel('x');
ylabel('y');
zlabel('z');
title('Fitted Bump Function');
hold on
scatter(x, y, 'b', 'filled');  % 'b' specifies blue color, 'filled' fills the markers
% Add legend with text of the initial parameters
legend_text = sprintf('ux=%d, nx=%d, uy=%d, ny=%d', round(optimal_ux,2), round(optimal_nx,2),...
    round(optimal_ny,2),round(optimal_uy,2));
legend(legend_text,'Generated Hand writing Data');
