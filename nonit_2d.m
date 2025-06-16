%% fixed point iteration
function [dnew,wnew,phinew,nit]=nonit_2d(dold,wold,phiold,dt,toll,disp,dk,eps1,gbc,eps2,oneconst,p,t,e,M,A,Ml) 
% fixed point iteration to solve nonlinear algebraic system in 2D, 1 time
% step of length dt
% dold, wold, phiold: d,w,phi at previous time step
% dt : time step size
% dx : mesh size (needed??)
% toll : tolerance until which iteration is run
% disp : sigma constant in front of w_t term
% dk : damping constant (alpha in paper)
% eps1, eps2 : coefficients in the source term and the elliptic equation
% (\eps_1, \eps_2)
% gbc : boundary condition interplated at desired time on mesh
% oneconst : elastic constant k
% p : nodes of mesh
% t : triangles of mesh
% e : boundary nodes
% M : mass matrix (maybe more efficient to generate it each time instead of
% handing it around?)
% A : standard stiffness matrix

% Output: 
% nit = number of iterations used

%boundary_nodes = extract_boundary_nodes(e);
%interior_nodes = extract_interior_nodes(p, e);

nit=1;
% first step of iteration (forward Euler kind of)
dtemp = compute_d_2d(dold,wold,dt); 
dnew = dtemp;
% dnew = interpolate_boundary_values(p, interior_nodes,boundary_nodes, dtemp); % neumann boundary conditions
phinew = solve_elliptic_nosource(dnew,gbc,eps2,p,t,e,M);
wnew=compute_w_2d(dnew,dold,wold,phiold,phinew,dt,dk,disp,eps1,oneconst,M,A,p,t,Ml);

% compute the initial error
% square of error
err = (wnew(:,1)-wold(:,1))'*Ml*(wnew(:,1)-wold(:,1)) ...
    + (dnew(:,1)-dold(:,1))'*A*(dnew(:,1)-dold(:,1)) ...
    + (wnew(:,2)-wold(:,2))'*Ml*(wnew(:,2)-wold(:,2)) ...
    + (dnew(:,2)-dold(:,2))'*A*(dnew(:,2)-dold(:,2)) ...
    + (wnew(:,3)-wold(:,3))'*Ml*(wnew(:,3)-wold(:,3)) ...
    + (dnew(:,3)-dold(:,3))'*A*(dnew(:,3)-dold(:,3)) ...
    + (phinew-phiold)'*A*(phinew-phiold);



% err= sum(sum(sum(dx^2*(wnew-wold).^2+...
%       (diff(dnew(1:end-1,2:end-1,:)-dold(1:end-1,2:end-1,:),1,1).^2)...
%       +(diff(dnew(2:end-1,1:end-1,:)-dold(2:end-1,1:end-1,:),1,2).^2) ))) ...
%       +sum(sum(diff(phinew-phiold,1,1).^2))+sum(sum(diff(phinew-phiold,1,2).^2)); 

% do the fixed point iteration     
    while(nit<200 && err>toll^2)  
        dmid=dnew;
        wmid=wnew;
        phimid = phinew;
        dtemp = compute_d_2d(dold,(wold+wmid)/2,dt); % compute d^{m,s+1}
        %dnew = interpolate_boundary_values(p, interior_nodes,boundary_nodes, dtemp); % neumann boundary conditions
        dnew = dtemp;
        phinew = solve_elliptic_nosource(dnew,gbc,eps2,p,t,e,M); % compute phi^{m,s+1}
        wnew=compute_w_2d(dnew,dold,wold,phiold,phinew,dt,dk,disp,eps1,oneconst,M,A,p,t,Ml); % compute w^{m,s+1}
        nit=nit+1;
    
        % compute the error in the nth iteration
        err = (wnew(:,1)-wmid(:,1))'*Ml*(wnew(:,1)-wmid(:,1)) ...
            + (dnew(:,1)-dmid(:,1))'*A*(dnew(:,1)-dmid(:,1)) ...
            + (wnew(:,2)-wmid(:,2))'*Ml*(wnew(:,2)-wmid(:,2)) ...
            + (dnew(:,2)-dmid(:,2))'*A*(dnew(:,2)-dmid(:,2)) ...
            + (wnew(:,3)-wmid(:,3))'*Ml*(wnew(:,3)-wmid(:,3)) ...
            + (dnew(:,3)-dmid(:,3))'*A*(dnew(:,3)-dmid(:,3)) ...
            + (phinew-phimid)'*A*(phinew-phimid);

        %err= sum(sum(sum(dx^2*(wnew-wmid).^2+...
        %    (diff(dnew(1:end-1,2:end-1,:)-dmid(1:end-1,2:end-1,:),1,1).^2)...
        %    +(diff(dnew(2:end-1,1:end-1,:)-dmid(2:end-1,1:end-1,:),1,2).^2) ))) ...
        %    +sum(sum(diff(phinew-phimid,1,1).^2))+sum(sum(diff(phinew-phimid,1,2).^2));  
    end

resterr = (wnew(:,1)-wmid(:,1))'*M*(wnew(:,1)-wmid(:,1)) ...
            + (dnew(:,1)-dmid(:,1))'*A*(dnew(:,1)-dmid(:,1)) ...
            + (wnew(:,2)-wmid(:,2))'*M*(wnew(:,2)-wmid(:,2)) ...
            + (dnew(:,2)-dmid(:,2))'*A*(dnew(:,2)-dmid(:,2)) ...
            + (wnew(:,3)-wmid(:,3))'*M*(wnew(:,3)-wmid(:,3)) ...
            + (dnew(:,3)-dmid(:,3))'*A*(dnew(:,3)-dmid(:,3));
phierr =  (phinew-phimid)'*A*(phinew-phimid);


% if more than 200 iterations are needed, print out the error and stop
% iteration
    if nit>=200
        fprintf('No convergence\n');
        fprintf(['error=',num2str(err),'\n']);
        fprintf(['phierror=',num2str(phierr),'\n']);
        fprintf(['resterror=',num2str(resterr),'\n']);       
    end   
    %fprintf(['Number of iterations',num2str(nit),'\n'])
end