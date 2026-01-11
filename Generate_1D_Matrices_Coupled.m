function [P, Q, S, D2, D1] = Generate_1D_Matrices_Coupled(N, dt)
    h = 1/N;
    e = ones(N+1, 1);
    
    % Laplacian (Standard Central Difference)
    D2 = spdiags([e -2*e e], -1:1, N+1, N+1);
    
    % Periodic BCs 
    % D2(1,N+1) = 1;
    % D2(N+1,1) = 1;
    % Neumann BCs
    D2(1,2) = 2;         
    D2(N+1,N) = 2;

    D2 = D2 / h^2;
    
    % Gradient (Central Difference) 
    % D1 u = (u_{i+1} - u_{i-1}) / (2h)
    D1 = spdiags([-e, e], [-1,1], N+1, N+1);
    
    % Periodic BCs for Gradient
    % D1(1,N+1) = -1; 
    % D1(N+1,1) = 1;
    % Neumann BCs
    D1(1,2) = 0;
    D1(N+1,N) = 0;
    D1 = D1 / (2*h);
    
    P = speye(N+1);
    Q = -dt * D2; 
    S = -speye(N+1);
end