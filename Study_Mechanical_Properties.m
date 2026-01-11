clear; clc; close all;

% Parameters
N = 128; T = 1000; ep = 0.01; dt = 1e-4; seed = 42; k = 2;
x = linspace(0, 1, N+1);
dx = 1/N;

% 1. Homogeneous: E_matrix = 100, E_precip = 100
p1.E_alpha = 100; p1.E_beta = 100; p1.eig_strain = 0.05;

% 2. Soft Precipitate: Matrix(alpha) = 100, Precip(beta) = 20
p2.E_alpha = 100; p2.E_beta = 20;  p2.eig_strain = 0.05;

% 3. Hard Precipitate: Matrix(alpha) = 20, Precip(beta) = 100
p3.E_alpha = 20;  p3.E_beta = 100; p3.eig_strain = 0.05;

% --- Run Simulations ---
fprintf('Running Homogeneous Case...\n');
[c1, u1, ~] = CH1D_Coupled_Generic(N, T, ep, dt, seed, k, p1);

fprintf('Running Soft Precipitate Case...\n');
[c2, u2, ~] = CH1D_Coupled_Generic(N, T, ep, dt, seed, k, p2);

fprintf('Running Hard Precipitate Case...\n');
[c3, u3, ~] = CH1D_Coupled_Generic(N, T, ep, dt, seed, k, p3);

% Total Strain (du/dx)
strain1 = gradient(u1(:,end)) / dx;
strain2 = gradient(u2(:,end)) / dx;
strain3 = gradient(u3(:,end)) / dx;

figure('Name', 'Mechanical Field Analysis', 'Color', 'w', 'Position', [100, 100, 1000, 800]);

% Plot 1: Concentration 
subplot(2,1,1); hold on;
plot(x, c1(:,end), 'k--', 'LineWidth', 1.5, 'DisplayName', 'Homogeneous');
plot(x, c2(:,end), 'b-', 'LineWidth', 2, 'DisplayName', 'Soft Precipitate');
plot(x, c3(:,end), 'r-', 'LineWidth', 2, 'DisplayName', 'Hard Precipitate');
title('(a) Final Concentration Profiles (c)');
ylabel('Concentration'); legend('Location', 'best'); grid on;
ylim([-1.2, 1.2]);

% Plot 2: Total Strain 
subplot(2,1,2); hold on;
plot(x, strain1, 'k--', 'LineWidth', 1.5, 'DisplayName', 'Homogeneous');
plot(x, strain2, 'b-', 'LineWidth', 2, 'DisplayName', 'Soft Precipitate');
plot(x, strain3, 'r-', 'LineWidth', 2, 'DisplayName', 'Hard Precipitate');
title('(b) Total Strain Distribution (\epsilon_{total} = du/dx)');
ylabel('Strain \epsilon'); xlabel('Position x');
legend('Location', 'best'); grid on;


disp('Plotting Complete.');