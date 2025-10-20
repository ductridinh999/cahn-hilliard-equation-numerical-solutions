function [P, Q, S, D] = Generate_1D_Matrices(N, mu)
    h = 1/N;
    e = ones(N+1, 1);
    
    % Second-order Laplacian D with Neumann BCs
    D = spdiags([e -2*e e], -1:1, N+1, N+1);
    % Periodic BCs
    D(1,N+1) = 1;
    D(N+1,1) = 1;
    % Newmann BCs
    %D(1,2) = 2;
    %D(N+1,N) = 2;
    D = D/h^2;
    
    % Define block matrices
    P = speye(N+1);
    Q = -mu*D; % incorrect from 2D, use dt instead
    S = -speye(N+1);
end