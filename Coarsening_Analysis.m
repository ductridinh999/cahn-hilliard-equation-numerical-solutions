% Script: Coarsening_Analysis.m
clear; clc; close all;

% Parameters
N = 128; 
T = 5000; % Long time to see coarsening
ep = 0.01; 
mu = 1; 
seed = 42; 
k = 1; % Random initial condition

% Run Simulation
disp('Running Simulation...');
[cvecs, t_comp] = CH1D12(N, T, ep, mu, seed, k);

% Calculate Energy
En = Energy_1D(cvecs, ep);

% Plot Results
figure('Name', 'Coarsening Dynamics', 'Color', 'w');
subplot(2,1,1);
imagesc(cvecs);
colormap jet; colorbar;
title('Phase Evolution (Space-Time Plot)');
xlabel('Time Step'); ylabel('Space x');

subplot(2,1,2);
plot(En, 'LineWidth', 2, 'Color', 'r');
title('Total Free Energy Evolution');
xlabel('Time Step'); ylabel('Energy');
grid on;