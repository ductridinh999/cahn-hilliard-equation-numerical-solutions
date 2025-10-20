function [cvecs, t] = CH1D12(N,T,ep,mu,seed,k)
    h = 1/N;
    dt = mu*h^2; 
    %dt = 0.01;
    [cvecs] = CH_initial_1D(N,T,k,seed);

    [P,~,S,D] = Generate_1D_Matrices(N,mu);
    
    Q = -dt*D;
    
    R = 2*speye(N+1) - (ep^2)*D;

    tic;
    for n = 2:T
        co = cvecs(:,n-1);
        A = [P,Q; R, S];
        b = [co; (3*co - co.^3)];
        x = A\b;
        cvecs(:,n) = x(1:(N+1));
    end
    t = toc;
end