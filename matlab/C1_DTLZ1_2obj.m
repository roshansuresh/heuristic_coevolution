%% 2-objective C1_DTLZ1
clear
close all
clc

%% Define problem parameters
k = 5;
n_obj = 2;
n_var = n_obj + k -1;
n_constr = 1;

% Variable sample points
n_samples_x = 11;
X_samples = linspace(0,1,n_samples_x);
X_des = combvec(X_samples, X_samples, X_samples, X_samples, X_samples, X_samples);

%% Compute objectives and constraints
f_X = zeros(size(X_des, 2), n_obj);
g_X = zeros(size(X_des, 2), n_constr);

for i = 1:size(X_des, 2)
    X_current = X_des(:, i);
    
    % Compute objectives
    g_xm_sum = 0;
    X_k = X_current(end-k+1:end, 1);
    for j = 1:k
        g_xm_sum = g_xm_sum + ((X_k(j) - 0.5)^2 - cos(20*pi*(X_current(j) - 0.5)));
    end
    g_xm = 100*(k + g_xm_sum);
    f1_x = (1/2)*X_current(1)*(1 + g_xm);
    f2_x = (1/2)*(1 - X_current(1))*(1 + g_xm);
    f_X(i, :) = [f1_x, f2_x];
    c_x = 1 - f2_x/0.6 - f1_x/0.5;
    g_X(i, :) = c_x; 
end

%% Plot objectives and constraints
figure
scatter(f_X(:,1), f_X(:,2),[],g_X,'filled')
colorbar
colormap jet

%% Plot only feasible designs
f_X_feas = f_X(g_X >= 0, :);
figure
scatter(f_X_feas(:, 1), f_X_feas(:, 2),'filled')