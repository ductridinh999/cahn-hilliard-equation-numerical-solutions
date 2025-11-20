clear all; close all; clc;

% Parameters
N = 100;        % Grid points
T = 3000;       % Time steps
ep = 0.01;      % Interface width
dt = 1e-5;      % Time step
seed = 42;      % Random seed
k = 1;          % Random initial condition, 2 for smooth cosine

% Run  Solver
[cvecs, u_vecs, time] = CH1D_Coupled(N, T, ep, dt, seed, k);

figure('Name', 'Coupled Evolution', 'Color', 'w');
% Subplot 1: Concentration Evolution
subplot(2,1,1);
hold on;
plot(linspace(0,1,N+1), cvecs(:,1), 'k--', 'LineWidth', 1, 'DisplayName', 't=0');
% plot(linspace(0,1,N+1), cvecs(:,round(T/2)), 'b', 'LineWidth', 1.5, 'DisplayName', ['t=',num2str(round(T/2))]);
plot(linspace(0,1,N+1), cvecs(:,end), 'r', 'LineWidth', 2, 'DisplayName', ['t=',num2str(T)]);
ylim([-1.1, 1.1]);
ylabel('Concentration c');
legend; title('Phase Separation');
grid on;

% Subplot 2: Displacement / Strain Field at final time
subplot(2,1,2);
x = linspace(0,1,N+1);
yyaxis left
plot(x, u_vecs(:,end), 'g-', 'LineWidth', 2);
ylabel('Displacement u');
yyaxis right
u_prime = gradient(u_vecs(:,end)) * N; 
plot(x, u_prime, 'm:', 'LineWidth', 1.5);
ylabel('Strain u''');
title('Mechanical State (Final Time)');
grid on;