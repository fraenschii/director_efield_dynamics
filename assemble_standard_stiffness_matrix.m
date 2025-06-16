function A = assemble_standard_stiffness_matrix(p, t)
    % Assemble the stiffness matrix for \int \Grad u\cdot \Grad v dx (no
    % coefficient), Neumann bc
    % p: nodal points coordinates
    % t: triangle indices
    
    
    N = size(p, 1);        % Number of nodes
    NT = size(t, 1);       % Number of triangles
    
    % Initialize matrices
    %A = sparse(N, N);
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
        
        % Calculate gradients of basis functions
        grad_phi = zeros(3, 2);
        grad_phi(1, :) = [y2-y3, x3-x2] / (2*area);
        grad_phi(2, :) = [y3-y1, x1-x3] / (2*area);
        grad_phi(3, :) = [y1-y2, x2-x1] / (2*area);
        
      
        % Calculate local stiffness matrix
        A_local = zeros(3, 3);
        
        for j = 1:3
            for k = 1:3
                % Compute grad_phi_j^T * grad_phi_k
                A_local(j, k) = grad_phi(j, :) * grad_phi(k, :)' * area;
            end
        end
        
       
        
        % Add to global matrices
        A(vertices,vertices) = A(vertices,vertices) + A_local;
        

    end
    
  
     
end
