function [p, t, e] = generate_mesh(n, domain)
    % Generate a structured triangular mesh on a square domain
    % p: nodal points coordinates
    % t: triangles (elements)
    % e: boundary edges
    
    xmin = domain(1,1); xmax = domain(1,2);
    ymin = domain(2,1); ymax = domain(2,2);
    
    h = (xmax - xmin) / n;  % Grid spacing
    
    % Generate nodal points
    [X, Y] = meshgrid(xmin:h:xmax, ymin:h:ymax);
    p = [X(:), Y(:)];
    
    % Number nodes
    nx = n + 1;
    ny = n + 1;
    N = nx * ny;
    
    % Generate triangles
    t = zeros(2*n*n, 3);
    for i = 1:n
        for j = 1:n
            % Node indices
            nd1 = j + (i-1)*nx;
            nd2 = j + i*nx;
            nd3 = j+1 + (i-1)*nx;
            nd4 = j+1 + i*nx;
            
            % Triangle 1
            idx = 2*((i-1)*n + j) - 1;
            t(idx, :) = [nd1, nd2, nd3];
            
            % Triangle 2
            t(idx+1, :) = [nd2, nd4, nd3];
        end
    end
    
    % Generate boundary edges
    e = [];
    
    % Bottom boundary
    for i = 1:n
        e = [e; i, i+1, 1];
    end
    
    % Right boundary
    for i = 1:n
        nd1 = (i+1)*nx;
        nd2 = (i)*nx;
        e = [e; nd1, nd2, 1];
    end
    
    % Top boundary
    for i = n:-1:1
        nd1 = i+1 + n*nx;
        nd2 = i + n*nx;
        e = [e; nd1, nd2, 1];
    end
    
    % Left boundary
    for i = n:-1:1
        nd1 = (i-1)*nx + 1;
        nd2 = i*nx + 1;
        e = [e; nd1, nd2, 1];
    end
end