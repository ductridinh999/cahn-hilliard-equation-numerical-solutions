function [cvecs] = CH_initial_1D(N, T, k, seed)
    cvecs = zeros(N+1, T);
    % Select initial condition
    switch(k)
        case 1 % Random initial condition
            rng(seed);
            cvecs(:,1) = 2*rand(N+1, 1) - 1;
        case 2 % Smooth cosine initial condition
            x = linspace(0,1,N+1)';
            cvecs(:,1) = cos(2*pi*x); 
    end
end