N = 80; 
T = 200; 
ep = 0.01; 
mu = 1; 
seed = 10; 
k = 1; % random / cosine initial condition

[cvecs, time] = CH1D12(N,T,ep,mu,seed,k);

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
