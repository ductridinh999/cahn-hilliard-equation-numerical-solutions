function [cvecs, t] = CH1D12(N, T, ep, dt, seed, k)
    % N: Grid points
    % T: Time steps
    % ep: Epsilon (interface width)
    % dt: Time step (Fixed input)
    % seed: Random seed
    % k: Initial condition flag

    h = 1/N;
    
    [cvecs] = CH_initial_1D(N, T, k, seed);
    
    [P, Q, S, D] = Generate_1D_Matrices(N,dt);
        
    R = 2*speye(N+1) - (ep^2)*D;

    tic;
    for n = 2:T
        co = cvecs(:,n-1);
        A = [P, Q; R, S];
        b = [co; (3*co - co.^3)];
        
        x = A \ b;
        cvecs(:,n) = x(1:(N+1));
    end
    t = toc;
end