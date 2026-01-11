function [cvecs, u_vecs, t] = CH1D_Coupled_Generic(N, T, ep, dt, seed, k, mech_params)
    % CH1D_Coupled_Generic: Solves Coupled CH-Elasticity with variable props
    % Inputs:
    %   mech_params: struct with fields .E_alpha, .E_beta, .eig_strain
    
    h = 1/N;
    
    % Unpack Mechanical Parameters
    E_alpha = mech_params.E_alpha;     % Stiffness of Phase -1
    E_beta  = mech_params.E_beta;      % Stiffness of Phase +1
    eig_strain = mech_params.eig_strain; % Magnitude of misfit
    
    % Interpolation Functions
    theta = @(c) (c + 1) / 2;
    dtheta_dc = @(c) 0.5 + 0*c;
    
    % Stiffness E(c)
    get_E = @(c) theta(c)*E_beta + (1-theta(c))*E_alpha;
    get_dEdc = @(c) dtheta_dc(c) * (E_beta - E_alpha);
    
    % Eigenstrain eps*(c)
    get_eps_star = @(c) eig_strain * c; 
    get_deps_star_dc = @(c) eig_strain + 0*c;

    % Initialize
    [cvecs] = CH_initial_1D(N, T, k, seed);
    u_vecs = zeros(N+1, T);
    
    % Get Matrices (Requires Generate_1D_Matrices_Coupled.m from your upload)
    [P, Q, S, D2, D1] = Generate_1D_Matrices_Coupled(N, dt);
    
    R_base = 2*speye(N+1) - (ep^2)*D2;
    
    tic;
    for n = 2:T
        c_old = cvecs(:,n-1);
        
        % --- 1. Solve Mechanics (Equilibrium: div(sigma) = 0) ---
        E_vals = get_E(c_old);
        dE_dx = D1 * E_vals;
        eps_star = get_eps_star(c_old);
        
        % Stiffness Matrix K = d/dx ( E(c) d/dx )
        K = spdiags(dE_dx, 0, N+1, N+1) * D1 + ...
            spdiags(E_vals, 0, N+1, N+1) * D2;
            
        % Force Vector from Eigenstrain (RHS)
        stress_free_term = E_vals .* eps_star;
        b_mech = D1 * stress_free_term;
        
        % BCs (Clamped)
        K(1,:) = 0; K(1,1) = 1; b_mech(1) = 0;
        K(end,:) = 0; K(end,end) = 1; b_mech(end) = 0;
        
        u = K \ b_mech;
        u_vecs(:,n) = u;
        
        % --- 2. Compute Elastic Potential ---
        u_prime = D1 * u; % Total Strain
        elastic_strain = u_prime - eps_star;
        dE_dc = get_dEdc(c_old);
        deps_dc = get_deps_star_dc(c_old);
        
        % Derivative of elastic energy density
        mu_el = 0.5 * dE_dc .* (elastic_strain.^2) + ...
                E_vals .* (elastic_strain) .* (-deps_dc);
            
        % --- 3. Solve Cahn-Hilliard ---        
        rhs_chem = (3*c_old - c_old.^3); 
        rhs_total = rhs_chem - mu_el; 
        
        A = [P, Q; R_base, S];
        b = [c_old; rhs_total];
        
        sol = A \ b;
        cvecs(:,n) = sol(1:N+1);
    end
    t = toc;
end