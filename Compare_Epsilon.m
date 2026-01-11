% Script: Compare_Epsilon.m
% Purpose: Study the effect of interfacial width parameter
clear; clc; close all;

eps_values = [0.005, 0.01, 0.02, 0.04];
N = 256;        % High resolution required for smallest epsilon
T = 500;       
dt = 1e-4;
seed = 42; 
k = 2;          % Cosine IC

figure('Name', 'Parameter Sensitivity: Epsilon', 'Color', 'w'); 
hold on;
c_map = lines(length(eps_values));

for i = 1:length(eps_values)
    ep = eps_values(i);
    [cvecs, ~] = CH1D12(N, T, ep, dt, seed, k);
    
    % Plot final concentration profile
    plot(linspace(0,1,N+1), cvecs(:,end), 'LineWidth', 2, ...
        'Color', c_map(i,:), ...
        'DisplayName', ['\epsilon = ' num2str(ep)]);
end

title('Effect of Gradient Energy Coefficient \epsilon');
xlabel('Domain x'); ylabel('Concentration c');
legend('Location', 'best'); grid on;
ylim([-1.1, 1.1]);