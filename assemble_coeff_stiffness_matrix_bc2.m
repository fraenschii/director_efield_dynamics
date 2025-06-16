function A_bc = assemble_coeff_stiffness_matrix_bc2(p, t,e, d,eps2,M)
    % Assemble the stiffness matrix for right hand side with boundary
    % conditions
    % p: nodal points coordinates
    % t: triangle indices
    % e: boundary edges
    % d: director field approximations at node values p
    % eps2: coefficient \eps_2 in equations
    % M is mass matrix

    N = size(p, 1);        % Number of nodes
    NT = size(t, 1);       % Number of triangles
    int_nodes = extract_interior_nodes(p, e);
    % Initialize matrices
    A = spalloc(N,N,6*N);
    
    % Assemble stiffness matrix and load vector
    for i = 1:NT
        % Get vertices of triangle
        vertices = t(i, :);
        x1 = p(vertices(1), 1); y1 = p(vertices(1), 2);
        x2 = p(vertices(2), 1); y2 = p(vertices(2), 2);
        x3 = p(vertices(3), 1); y3 = p(vertices(3), 2);
        
        % Calculate area of triangle
        area = 0.5 * abs((x2-x1)*(y3-y1) - (x3-x1)*(y2-y1));
        
        % get values of d
        d_loc = d(vertices,1:2);
        % get local mass matrix
        M_loc = M(vertices,vertices);
        % Calculate gradients of basis functions
        grad_b = zeros(3, 2);
        grad_b(1, :) = [y2-y3, x3-x2] / (2*area);
        grad_b(2, :) = [y3-y1, x1-x3] / (2*area);
        grad_b(3, :) = [y1-y2, x2-x1] / (2*area);
       
        % Evaluate matrix-valued coefficient 
        % anisotropic part
        c_matrix = eps2*d_loc'*M_loc*d_loc;

        % isotropic part
        iso = area*eye(2);
        
        % Calculate local stiffness matrix
        A_local = zeros(3, 3);
        for j = 1:3
            for k = 1:3
                % Compute grad_phi_j^T * c_matrix * grad_phi_k
                A_local(j, k) = grad_b(j, :) * (iso+c_matrix) * grad_b(k, :)';
            end
        end
        
        % Add to global matrices
        A(vertices,vertices) = A(vertices,vertices) + A_local;
       
    end
    
    A_bc = A(int_nodes,:);

end
