function [cvecs, u_vecs, t] = CH1D_Coupled(N, T, ep, dt, seed, k)
    h = 1/N;
    
    E_alpha = 100;  % Stiff phase
    E_beta  = 50;   % Soft phase
    eig_strain = 0.05; % Misfit strain magnitude
    
    % theta(c): 0 at c=-1, 1 at c=1
    theta = @(c) (c + 1) / 2;
    dtheta_dc = @(c) 0.5 + 0*c;
    
    % Young's Modulus E(c)
    get_E = @(c) theta(c)*E_beta + (1-theta(c))*E_alpha;
    get_dEdc = @(c) dtheta_dc(c) * (E_beta - E_alpha);
    
    % Eigenstrain eps*(c) = alpha * dT
    get_eps_star = @(c) eig_strain * c; 
    get_deps_star_dc = @(c) eig_strain + 0*c;

    % Initialize
    [cvecs] = CH_initial_1D(N, T, k, seed);
    u_vecs = zeros(N+1, T);
    
    % Get Matrices (D2=L_h, D1=D_h)
    [P, Q, S, D2, D1] = Generate_1D_Matrices_Coupled(N, dt);
    
    % Constant parts)
    R_base = 2*speye(N+1) - (ep^2)*D2;
    
    tic;
    for n = 2:T
        c_old = cvecs(:,n-1);
        
        % 1. Solve Mechanics, calculate based on c_old
        E_vals = get_E(c_old);
        dE_dx = D1 * E_vals;
        eps_star = get_eps_star(c_old);
        % LHS
        K = spdiags(dE_dx, 0, N+1, N+1) * D1 + ...
            spdiags(E_vals, 0, N+1, N+1) * D2;
        % RHS
        stress_free_term = E_vals .* eps_star;
        b_mech = D1 * stress_free_term;
        
        % Mechanical BCs (clamped)
        K(1,:) = 0; K(1,1) = 1; b_mech(1) = 0;
        K(end,:) = 0; K(end,end) = 1; b_mech(end) = 0;
        % Solve for u
        u = K \ b_mech;
        u_vecs(:,n) = u;
        
        % 2. Compute mu_el
        u_prime = D1 * u;
        elastic_strain = u_prime - eps_star;
        dE_dc = get_dEdc(c_old);
        deps_dc = get_deps_star_dc(c_old);
        
        mu_el = 0.5 * dE_dc .* (elastic_strain.^2) + ...
                E_vals .* (elastic_strain) .* (-deps_dc);
            
        % 3. Solve the Coupled CH        
        rhs_chem = (3*c_old - c_old.^3); 
        rhs_total = rhs_chem - mu_el; % Add elastic potential 
        
        A = [P, Q; R_base, S];
        b = [c_old; rhs_total];
        
        sol = A \ b;
        cvecs(:,n) = sol(1:N+1);
    end
    t = toc;
end