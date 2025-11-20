function [cvecs, t] = CH1D12_Variable_Mobility(N,T,ep,dt,seed,k)
    h = 1/N;
    
    [cvecs] = CH_initial_1D(N,T,k,seed); 
    [P, ~, S, D] = Generate_1D_Matrices(N, 1); %
    
    R = 2*speye(N+1) - (ep^2)*D;

    % Degenerate mobility M(c) = 1 - c^2
    M_func = @(c) max(0, 1 - c.^2); 
    
    % Constant mobility M(c) = 1 ---
    % M_func = @(c) 1 + 0*c; 

    tic;
    for n = 2:T
        co = cvecs(:,n-1);
       
        M = M_func(co);

        M_plus = circshift(M, -1);  
        M_minus = circshift(M, 1);
        
        M_avg_super = (M + M_plus) / 2;  
        M_avg_sub = (M + M_minus) / 2; 
        
        band_sub = M_avg_sub / h^2;
        band_super = M_avg_super / h^2;
        band_main = -(band_sub + band_super);
        
        bands = [band_sub, band_main, band_super];
        L_M = spdiags(bands, [-1, 0, 1], N+1, N+1);
        
        L_M(1, N+1) = band_sub(1);
        L_M(N+1, 1) = band_super(N+1);

        Q = -dt * L_M; 
        A = [P, Q; R, S];
        
        b = [co; (3*co - co.^3)];
        x = A\b;
        cvecs(:,n) = x(1:(N+1));
    end
    t = toc;
end