function [cvecs] = CH_initial_1D(N, T, k, seed)
    cvecs = zeros(N+1, T);
    switch(k)
        case 1 % Random initial condition
            rng(seed);
            cvecs(:,1) = 0.1*(2*rand(N+1, 1) - 1);
        case 2 % Smooth cosine initial condition
            x = linspace(0,1,N+1)';
            cvecs(:,1) = cos(6*pi*x)*0.1; 
        case 3 % Linear initial condition 
            x = linspace(0,1,N+1)';           
            cvecs(:,1) = 1.6 * x - 0.8;
        case 4 %  Sinusoidal
            dx = 1/N;
            for i = 3:N-2
                cvecs(i,1) = 0.1 * sin(2*pi*(i-3)*dx) ...
                      + 0.01 * cos(4*pi*(i-3)*dx) ...
                      + 0.06 * sin(4*pi*(i-3)*dx) ...
                      + 0.02 * cos(10*pi*(i-3)*dx);
            end
    end
end