function interior_nodes = extract_interior_nodes(p, e)
    % Extract interior nodes
    % p: all node coordinates
    % e: boundary edges array
    
    % Get all node indices
    all_nodes = 1:size(p, 1);
    
    % Get boundary nodes
    boundary_nodes = unique(e(:,1:2));
    
    % Interior nodes are those not on the boundary
    interior_nodes = setdiff(all_nodes, boundary_nodes);
end