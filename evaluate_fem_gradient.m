function [grad_u] = evaluate_fem_gradient(p, t, u)
    % Evaluate the gradient of a finite element approximation at grid points
    %
    % Parameters:
    % p: Node coordinates (Nx2 matrix)
    % t: Triangle connectivity (Mx3 matrix)
    % u: Solution values at nodes (N x 1 vector)
    %
    % Returns:
    % grad_u: Gradient of the solution at each node (Nx2 matrix)
    %         Each row is [du/dx, du/dy] at the corresponding node
    
    % Number of nodes
    N = size(p, 1);
    
    % Initialize gradient storage
    grad_u = zeros(N, 2);
    
    % Keep track of the number of triangles contributing to each node
    node_counts = zeros(N, 1);
    
    % Iterate through triangles
    for i = 1:size(t, 1)
        % Get vertices of current triangle
        vertices = t(i, :);
        
        % Extract coordinates of triangle vertices
        x1 = p(vertices(1), 1); y1 = p(vertices(1), 2);
        x2 = p(vertices(2), 1); y2 = p(vertices(2), 2);
        x3 = p(vertices(3), 1); y3 = p(vertices(3), 2);
        
        % Calculate area of triangle
        area = 0.5 * abs((x2-x1)*(y3-y1) - (x3-x1)*(y2-y1));
        
        % Calculate gradients of basis functions
        % These are constant over the triangle
        b1 = [y2-y3, x3-x2] / (2*area);
        b2 = [y3-y1, x1-x3] / (2*area);
        b3 = [y1-y2, x2-x1] / (2*area);
        
        % Get solution values at triangle vertices
        u1 = u(vertices(1));
        u2 = u(vertices(2));
        u3 = u(vertices(3));
        
        % Compute gradient using linear combination of basis function gradients
        grad_triangle = u1 * b1 + u2 * b2 + u3 * b3;
        
        % Accumulate gradient for each node of the triangle
        for j = 1:3
            node = vertices(j);
            grad_u(node, :) = grad_u(node, :) + grad_triangle;
            node_counts(node) = node_counts(node) + 1;
        end
    end
    
    % Average the gradients (to handle nodes shared by multiple triangles)
    for i = 1:N
        grad_u(i, :) = grad_u(i, :) / node_counts(i);
    end
end