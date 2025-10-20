function En = Energy_1D(C,ep)

% Get dimensions
[n_pts,T] = size(C); 
h = 1/(n_pts-1);

% Energy vector
En = zeros(T,1);

% Simpson weights
s = zeros(n_pts,1);
s(1) = 1; s(n_pts) = 1; s(2:2:end-1) = 4; s(3:2:end-2) = 2;
S = s*(h/3);

% First-order differentiation matrix (for gradient term)
e = ones(n_pts,1);
D1 = spdiags([-e, e], [-1,1], n_pts, n_pts)/(2*h);
%D1(1,2) = 0; D1(n_pts,n_pts-1) = 0; % Neumann BCs
D1(1,n_pts) = 0;
D1(n_pts,1) = 0;
g = @(c) (1 - c.^2).^2/4 + ep^2/2 * (D1*c).^2;
    for tt = 1:T
        En(tt) = S'*g(C(:,tt));
    end
    
end