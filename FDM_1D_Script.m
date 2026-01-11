N = 128; 
T = 5000; 
ep = 0.01; 
dt = 1e-5; 
seed = 42; 
k = 2; % random / cosine initial condition

[cvecs, time] = CH1D12(N,T,ep,dt,seed,k);

% Plot the evolution
CH1D_Plot_Evolution(cvecs);

% Calculate energy
En = Energy_1D(cvecs,ep);

% Plot energy
figure();
plot(En, 'LineWidth', 2); 
title('Energy Evolution'); 
ylabel('Energy'); 
xlabel('Time-step');
grid on;
