function phinew = solve_elliptic_nosource(d,gbc,eps2,p,t,e,M) 
% solve the elliptic equation for phi with given coefficients and director
% field d, right hand side only involves boundary condtions (no free charges)
% M mass matrix
% d director field
% eps2 coefficient in equation
% gbc boundary condition interpolated on the mesh
% p vertices
% t triangles
% e boundary edges
    
    % Assemble system
    A = assemble_coeff_stiffness_matrix_homDirichlet(p, t, e, d,eps2,M);
    A_bc = assemble_coeff_stiffness_matrix_bc2(p, t, e,d,eps2,M);
    int_nodes = extract_interior_nodes(p,e);
    rhs = -A_bc*gbc;
    N = size(p,1);
    phinew = zeros(N,1);
    
    % Solve system
    u0 = A\rhs;

    phinew(int_nodes) = u0;
    phinew = phinew + gbc; % phi = u_0 + \tilde{g}
end

