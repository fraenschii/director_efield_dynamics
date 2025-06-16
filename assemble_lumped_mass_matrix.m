function M = assemble_lumped_mass_matrix(p, t)
    % Assemble the lumped mass matrix for a FEM problem
    % M_ij = \int b_i \delta_{ij} dx
    % 
    % Parameters:
    % p: nodal points coordinates
    % t: triangle indices
    %
    % Returns:
    % M: lumped mass matrix (sparse)
    
    N = size(p, 1);        % Number of nodes
    NT = size(t, 1);       % Number of triangles
    
    % Initialize lumped mass matrix
    
    M = spalloc(N,N,N);
   

    Mlumped_local =1/6;
    % Assemble mass matrix
    for i = 1:NT
        % Get vertices of triangle
        vertices = t(i, :);
        x1 = p(vertices(1), 1); y1 = p(vertices(1), 2);
        x2 = p(vertices(2), 1); y2 = p(vertices(2), 2);
        x3 = p(vertices(3), 1); y3 = p(vertices(3), 2);
        
        % Calculate area of triangle
        detJk =   abs((x2-x1)*(y3-y1) - (x3-x1)*(y2-y1));
        
        % Add to global mass matrix
        for j = 1:3
            %for k = 1:3
                M(vertices(j), vertices(j)) = M(vertices(j), vertices(j)) + Mlumped_local*detJk;
            %end
        end
    end
end

