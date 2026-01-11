clear; clc; close all;

k_values = [1, 2, 3];
k_names = {'Random Noise', 'Cosine Wave', 'Linear Gradient'};
N = 128; T = 5000; ep = 0.01; dt = 1e-5; seed = 42;

figure('Name', 'Initial Conditions', 'Color', 'w');

for i = 1:3
    k = k_values(i);
    [cvecs, ~] = CH1D12(N, T, ep, dt, seed, k);
    
    subplot(3,1,i);
    plot(cvecs(:,1), 'k--', 'DisplayName', 'Initial'); hold on;
    plot(cvecs(:,end), 'r-', 'LineWidth', 2, 'DisplayName', 'Final');
    title(['Condition: ' k_names{i}]);
    legend; grid on; axis tight;
end