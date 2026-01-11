clear; clc; close all;

N_values = [16, 32, 64, 128, 256, 512]; 
T = 500;        % Time steps 
ep = 0.01;      % Interface width parameter
dt = 1e-4;      % Time step
seed = 42; 
k = 2;          % Cosine IC 

energy_results = {};
time_axis = (1:T) * dt; 

c_map = lines(length(N_values)); 

for i = 1:length(N_values)
    N = N_values(i);
    fprintf('Running Simulation for N = %d ...\n', N);
    
    [cvecs, ~] = CH1D12(N, T, ep, dt, seed, k);
    
    % Calculate Free Energy
    En = Energy_1D(cvecs, ep);
    energy_results{i} = En;
end

figure('Name', 'Mesh Convergence Study', 'Color', 'w', 'Position', [100, 100, 800, 700]);

subplot(2,1,1); hold on;
for i = 1:length(N_values)
    plot(time_axis, energy_results{i}, 'LineWidth', 1.5, ...
        'Color', c_map(i,:), ...
        'DisplayName', ['N = ' num2str(N_values(i))]);
end
title('\bf(a) Total Free Energy Evolution');
ylabel('Free Energy'); 
xlabel('Time t');
legend('Location', 'northeast'); 
grid on; axis tight;

subplot(2,1,2); hold on;

% baseline mesh (N=128)
ref_N = 128;
ref_idx = find(N_values == ref_N);

if isempty(ref_idx)
    error(['Baseline N=' num2str(ref_N) ' not found in N_values array.']);
end

E_ref = energy_results{ref_idx}; % Reference solution (N=128)

for i = 1:length(N_values)
    if i == ref_idx
        continue; % Skip plotting the reference against itself
    end
    
    E_curr = energy_results{i};
    
    % Difference: E_ref - E_curr
    diff_E = E_ref - E_curr; 
    
    plot(time_axis, diff_E, 'LineWidth', 1.5, ...
        'Color', c_map(i,:), ...
        'DisplayName', ['N = ' num2str(N_values(i))]);
end

yline(0, 'k--', 'HandleVisibility', 'off'); xw
title(['\bf(b) Energy Error vs Converged Baseline (N=' num2str(ref_N) ')']);
ylabel('Energy Difference (E_{ref} - E_{N})'); 
xlabel('Time t');
legend('Location', 'best'); 
grid on; axis tight;