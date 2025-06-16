function S = compute_source(p, t, d,eps1,phiold,phinew,M)
    % compute source term in w equation
    % Input:
    % p: nodal points coordinates
    % t: triangle indices
    % d: director field approximations at node values p
    % eps1: coefficient \eps_1 in equations
    % phiold, phinew : phi at previous time step, latest update of phi
    % M is mass matrix

    % Output: S: source term in w equation at nodes
    
    N = size(p, 1);        % Number of nodes
    NT = size(t, 1);       % Number of triangles
    SA = zeros(N,3);
    SB = zeros(N,3);
    
    % Compute contribution of source on each triangle
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
        % get values of phiold
        phiold_loc = phiold(vertices);
        % get values of phinew locally
        phinew_loc = phinew(vertices);
        % compute cross products (since we assume that \partial_{x_3} \phi
        % = 0, the last component is not needed)
        cross1 = [zeros(3,1),d(vertices,3)]; 
        cross2 = [-d(vertices,3),zeros(3,1)]; 
        cross3 = [d(vertices,2),-d(vertices,1)];
        % Calculate gradients of basis functions
        grad_b = zeros(3, 2);
        grad_b(1, :) = [y2-y3, x3-x2] / (2*area);
        grad_b(2, :) = [y3-y1, x1-x3] / (2*area);
        grad_b(3, :) = [y1-y2, x2-x1] / (2*area);
        

        % compute local source term
        cross1A = phinew_loc'*(grad_b*cross1');
        cross1B = phiold_loc'*(grad_b*cross1');
        cross2A = phinew_loc'*(grad_b*cross2');
        cross2B = phiold_loc'*(grad_b*cross2');
        cross3A = phinew_loc'*(grad_b*cross3');
        cross3B = phiold_loc'*(grad_b*cross3');

        tempA = phiold_loc'*(grad_b*d_loc')*M_loc;
        tempB = phinew_loc'*(grad_b*d_loc')*M_loc;

        for j=1:3
            SA_loc(j,1)=cross1A(j)*tempA(j);
            SA_loc(j,2)=cross2A(j)*tempA(j);
            SA_loc(j,3)=cross3A(j)*tempA(j);
            SB_loc(j,1)=cross1B(j)*tempB(j);
            SB_loc(j,2)=cross2B(j)*tempB(j);
            SB_loc(j,3)=cross3B(j)*tempB(j);
        end
       
        % add to global source terms
       SA(vertices,:) = SA(vertices,:) +SA_loc;
       SB(vertices,:) = SB(vertices,:) +SB_loc;       
       
    end
    
    S = eps1*(SA + SB)/2;
end
