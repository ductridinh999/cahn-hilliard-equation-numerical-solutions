clear; clc; close all;

% --- Parameters ---
N = 128;
T = 300;
ep = 0.01;
dt = 1e-4; 
seed = 42;
k = 2; % Cosine initial condition

test_values = [1, 5, 10]; 

results_const = {};
results_var = {};

fprintf('Running Constant Mobility Cases (CH1D13)...\n');
for i = 1:length(test_values)
    val = test_values(i);
    [cvecs, ~] = CH1D13(N, T, ep, dt, seed, k, val);
    results_const{i} = Energy_1D(cvecs, ep);
end

fprintf('Running Variable Mobility Cases (CH1D12_VarMob_D)...\n');
for i = 1:length(test_values)
    val = test_values(i);
    [cvecs, ~] = CH1D12_Variable_Mobility(N, T, ep, dt, seed, k, val);
    results_var{i} = Energy_1D(cvecs, ep);
end

c_map = lines(length(test_values));

% FIGURE 1: Separate Analyses 
figure('Name', 'Mobility Analysis: Separate', 'Color', 'w', 'Position', [100, 100, 1000, 500]);

subplot(1,2,1); hold on;
for i = 1:length(test_values)
    plot(results_const{i}, 'LineWidth', 2, 'Color', c_map(i,:), ...
        'DisplayName', ['M = ' num2str(test_values(i))]);
end
title('Constant Mobility M');
xlabel('Time Steps'); ylabel('Free Energy');
set(gca, 'YScale', 'log'); legend('Location', 'best'); grid on; axis tight;

subplot(1,2,2); hold on;
for i = 1:length(test_values)
    plot(results_var{i}, '--', 'LineWidth', 2, 'Color', c_map(i,:), ...
        'DisplayName', ['D = ' num2str(test_values(i))]);
end
title('Variable Mobility M(c) = D(1-c^2)');
xlabel('Time Steps'); ylabel('Free Energy');
set(gca, 'YScale', 'log'); legend('Location', 'best'); grid on; axis tight;


% FIGURE 2: Combined Overlay 
figure('Name', 'Mobility Analysis: Combined Overlay', 'Color', 'w', 'Position', [150, 150, 700, 500]);
hold on;

for i = 1:length(test_values)
    val = test_values(i);
    col = c_map(i,:);
    
    % Plot Constant M (Solid Line)
    plot(results_const{i}, '-', 'LineWidth', 2, 'Color', col, ...
        'DisplayName', ['Const M = ' num2str(val)]);
        
    % Plot Variable M (Dashed Line)
    plot(results_var{i}, '--', 'LineWidth', 2, 'Color', col, ...
        'DisplayName', ['Var D = ' num2str(val)]);
end

title('Direct Comparison: Constant vs Variable Mobility');
xlabel('Time Steps'); ylabel('Free Energy');
set(gca, 'YScale', 'log'); 
grid on; axis tight;

legend('Location', 'northeast');
