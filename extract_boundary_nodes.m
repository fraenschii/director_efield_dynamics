function boundary_nodes = extract_boundary_nodes(e)
    % Extract unique boundary nodes from the edge array
    % e: boundary edges array from generate_mesh
    
    % Extract all boundary nodes from the first two columns of e
    boundary_nodes = unique(e(:,1:2));
end