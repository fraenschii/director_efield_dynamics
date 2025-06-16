function M = assemble_mass_matrix(p, t)
    % Assemble the mass matrix for a FEM problem
    % M_ij = \int b_i b_j dx
    % 
    % Parameters:
    % p: nodal points coordinates
    % t: triangle indices
    %
    % Returns:
    % M: mass matrix (sparse)
    
    N = size(p, 1);        % Number of nodes
    NT = size(t, 1);       % Number of triangles
    
    % Initialize mass matrix
    
    M = spalloc(N,N,6*N);
    % Standard mass matrix for linear elements on a triangle
    % For local basis functions on reference triangle
    M_local = [2 1 1; 1 2 1; 1 1 2] / 12;
    
    % Assemble mass matrix
    for i = 1:NT
        % Get vertices of triangle
        vertices = t(i, :);
        x1 = p(vertices(1), 1); y1 = p(vertices(1), 2);
        x2 = p(vertices(2), 1); y2 = p(vertices(2), 2);
        x3 = p(vertices(3), 1); y3 = p(vertices(3), 2);
        
        % Calculate area of triangle
        area = 0.5 * abs((x2-x1)*(y3-y1) - (x3-x1)*(y2-y1));
        
        % Scale local mass matrix by the area
        M_elem = M_local * area;
        
        % Add to global mass matrix
        M(vertices,vertices) = M(vertices,vertices)+ M_elem;
       
    end
end

