clear all; close all; clc;

N = 80; 
T = 200; 
ep = 0.01; 
seed = 10; 
k = 1; 

h = 1/N;
dt = 0.0001; 

[cvecs, time] = CH1D12_Variable_Mobility(N,T,ep,dt,seed,k); 

disp(['Simulation time: ', num2str(time), ' seconds']);

% Plot the evolution
CH1D_Plot_Evolution(cvecs); %

% Calculate energy
En = Energy_1D(cvecs,ep); %

% Plot energy
figure(); 
plot(En, 'LineWidth', 2); 
title('Energy Evolution'); 
ylabel('Energy'); 
xlabel('Time-step');
grid on;