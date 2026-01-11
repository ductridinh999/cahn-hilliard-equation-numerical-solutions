function [cvecs, t] = CH1D13(N, T, ep, dt, seed, k, M)
    % Inputs:
    %   M: Mobility scalar (e.g., 1, 5, 10)
    
    h = 1/N;
    
    % Initialize
    [cvecs] = CH_initial_1D(N, T, k, seed);
    
    % Generate Base Matrices
    % Q_base = -dt * Laplacian
    [P, Q_base, S, D] = Generate_1D_Matrices(N, dt);
    
    % Apply Mobility M to the diffusion term
    Q = M * Q_base; 
        
    % R matrix (Chemical Potential Update)
    R = 2*speye(N+1) - (ep^2)*D;

    tic;
    for n = 2:T
        co = cvecs(:,n-1);
        
        % System: 
        % [ P   Q ] [ c_new ] = [ c_old ]
        % [ R   S ] [ mu_new]   [ f(c_old) ]
        
        A = [P, Q; R, S];
        b = [co; (3*co - co.^3)];
        
        x = A \ b;
        cvecs(:,n) = x(1:(N+1));
    end
    t = toc;
end