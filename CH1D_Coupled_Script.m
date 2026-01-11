clear all; close all; clc;

N = 128;        
T = 500;        % Total time steps
ep = 0.01;      
dt = 1e-4;      
seed = 42;      
k = 2;          % Random noise

disp('Running Coupled Simulation...');
[cvecs, u_vecs, time] = CH1D_Coupled(N, T, ep, dt, seed, k);
disp('Done.');

indices = [1, 50, 200, T];
stage_names = {'Initial', 'Spinodal Decomposition', 'Coarsening', 'Final State'};

figure('Name', 'Coupled Evolution Stages (T=500)', 'Color', 'w', 'Position', [100, 50, 700, 1000]);

x = linspace(0, 1, N+1);

for i = 1:4
    idx = indices(i);
    
    subplot(4, 1, i);
    
    yyaxis left
    plot(x, cvecs(:, idx), 'b-', 'LineWidth', 2);
    ylabel('Conc. c', 'Color', 'b');
    set(gca, 'YColor', 'b');
    ylim([-1.1, 1.1]);
    
    yyaxis right
    % Calculate strain (gradient of displacement)
    u_prime = gradient(u_vecs(:, idx)) * N; 
    plot(x, u_prime, 'r--', 'LineWidth', 1.5);
    ylabel('Strain \epsilon', 'Color', 'r');
    set(gca, 'YColor', 'r');
    
    title([stage_names{i}, ' (Step ', num2str(idx), ')']);
    grid on;
    
    if i == 4
        xlabel('Position x');
    end
end