%% This example fit a Bump Function over 1d horizontal coordinates handwriting dataset
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
        [filepath] = fullfile('C:\Users\aboschi\OneDrive - Fondazione Istituto Italiano Tecnologia\Documents\MATLAB\Handwriting\Prove\Pazienti\pz3_racconto.csv');  % Select data file
        orient = 'l';  % Variable defining tablet orientation
        % Read Data
        wd = read(wacomdata,filepath,orient); % Call the wacomdata class constructor
        test = wd.position;
        test(1300:end,:) = []; % Select 1 chunk from the text (first row)
        %data2find
        Data = test(:,1);  % Data

        figure;
        histogram(Data, round(length(Data)./10), 'Normalization', 'pdf', 'FaceColor', 'r', 'EdgeColor', 'none');
        xlabel('x');
        ylabel('Probability Density');
        title('Chunk PDF');
        legend('X Data Points');
        hold off;

    case 'generated'
        ux = 10;
        nx = 30;
        n = 1000;
        Bnx = integral(@(x) Bump1(x, ux,nx),ux-nx,ux+nx);
        % Define the range for generating random points
        xmin = ux-nx;
        xmax = ux+nx;
        ymax = max(Bump1(linspace(xmin, xmax, n),ux,nx)./Bnx); % Upper bound for y
        
        % Initialize variables to store accepted points
        accepted_points = zeros(1, 0);
        
        % Generate random points and accept/reject based on the PDF
        while numel(accepted_points) < n
            % Generate random point coordinates
            x_rand = (xmax - xmin) * rand() + xmin;
            y_rand = ymax * rand();
            
            % Evaluate the PDF at the random point
            pdf_value = Bump1(x_rand,ux,nx)./Bnx;
            
            % Accept or reject the point
            if y_rand < pdf_value
                % Accept the point
                accepted_points(end + 1) = x_rand;
            end
        end
        
        % Plot the accepted points along with the PDF
        x_values = linspace(xmin, xmax, 1000);
        pdf_values = Bump1(x_values,ux,nx)./Bnx;
        Data = accepted_points;
        figure;
        plot(x_values, pdf_values, 'b-', 'LineWidth', 2);
        hold on;
        histogram(accepted_points,n/5, 'Normalization', 'pdf', 'FaceColor', 'r', 'EdgeColor', 'none');
        xlabel('x');
        ylabel('Probability Density');
        title('Histogram of Data Points and Bump Function PDF');
        legend('Bump Function PDF', 'Generated Data Points');
        hold off;

end
%% Optimization problem

% Initial guess for the variables ux,nx,uy,ny,teta
% iparams(1) =  -125;%ux
% iparams(2) = 125; %nx
iparams(1)= mean(Data)+3.*std(Data); % ux
iparams(2) = 3.*std(Data); % nx
 
% Bounds on   ux: Lmin <= ux <= Lmax
            % nx: Lmin/2 <= nx <= Lmax/2
            % uy: Hmin <= uy <= Hmax
            % ny: Hmin/2 <= ny <= Hmax/2
            % teta: 0 <= teta <= 90;
lb = [min(Data)-30 0 ];
ub = [max(Data)+30 (max(Data)+30)+3.*std(Data)];

% inequality const
A = [-1 -1  ;...
      1 -1  ;...
      0 -1  ];

b = [-max(Data)+0.001 min(Data)+0.001  -10];

% check if initial point satisfy constrain
% Inizialization
if (any((A*iparams') < b' == 0))

    f = zeros(size(x0)); % assumes x0 is the initial point

    xnew = linprog(f,A,b,[],[],lb,ub);  % Solve the linear programming problem to see if there is a feasible point
    x0 = xnew;
else 
    x0 = iparams;
end

%% solving Optimization Problem

% Objective function: f(x) = x^2
obj = @(params) objectiveFunction1D(params,Data);

% Define the problem structure
fun = obj;


% Options for the solver (optional)
options = optimoptions('fmincon', 'Display', 'iter');

% Solve the problem
[optimal_params] = fmincon(fun, xnew, A, b, [], [], lb, ub, [], options);

%results
optimal_ux = optimal_params(1);
optimal_nx = optimal_params(2);

%% Plot Results
% Plot the accepted points along with the PDF
x_values = linspace(optimal_ux-optimal_nx, optimal_ux+optimal_nx, 1000);

Bnx = integral( @(x) Bump1(x,optimal_ux,optimal_nx), optimal_ux-optimal_nx , optimal_ux+optimal_nx);

pdf_values = Bump1(x_values,optimal_ux,optimal_nx)./Bnx;

figure;
plot(x_values, pdf_values, 'b-', 'LineWidth', 2);
hold on;
histogram(Data,20, 'Normalization', 'pdf', 'FaceColor', 'r', 'EdgeColor', 'none');
xlabel('x');
ylabel('Probability Density');
title('Histogram of Data Points and Bump Function PDF');
legend('Fitted Bump Function PDF', 'Data Points');
hold off;
