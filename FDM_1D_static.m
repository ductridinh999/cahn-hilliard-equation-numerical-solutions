clear; clc; close all;

N = 128;            % Grid resolution
T = 500;            % Total time steps 
ep = 0.01;          % Interface width
dt = 1e-4;          % Time step
seed = 42;          
k = 1;              % Random Noise 

disp('Running Simulation...');
[cvecs, time] = CH1D12(N, T, ep, dt, seed, k);
disp('Done.');

indices = [1, 50, 200, T]; 
titles = {'Stage 1: Initial Conditions', ...
          'Stage 2: Spinodal Decomposition', ...
          'Stage 3: Coarsening', ...
          'Stage 4: Final State'};

figure('Name', 'Evolution Stages (T=500)', 'Color', 'w', 'Position', [100, 50, 600, 900]);
x_grid = linspace(0, 1, N+1);

for i = 1:4
    idx = indices(i);
    subplot(4, 1, i);
    
    % Plot Concentration
    plot(x_grid, cvecs(:, idx), 'LineWidth', 2, 'Color', 'b');
    
    % Formatting
    ylim([-1.1, 1.1]);
    yline(0, 'k:', 'LineWidth', 1); % Zero reference
    title([titles{i}, ' (Step ', num2str(idx), ')']);
    ylabel('Conc. \phi');
    grid on;
    
    if i == 4
        xlabel('Position x');
    end
end