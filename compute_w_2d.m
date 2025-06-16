%% w evolution
function wn=compute_w_2d(dn,dp,w,phiold,phinew,dt,dk,disp,eps1,oneconst,M,A,p,t,Ml)
% update for w equation 
% Input:
% dp: old nodal values for d
% dn: new nodal values for d
% w: old nodal values for w
% phiold: old nodal values for phi
% phinew: new nodal values for phi
% dt: time step
% dk: damping constant (alpha in paper)
% disp: inertia constant (sigma in paper)
% eps1: constant in front of source term involving electric field
% oneconst: elastic constant (k in paper)
% M : standard mass matrix
% A : standard stiffness matrix of problem (Neumann bc)
% p : coordinates of nodes of mesh
% t : triangles of mesh with nodes
% Ml : mass-lumped matrix

% Output:
% wn: nodal values for w at next iteration step

fracm=1/(disp+dk*dt/2);
fracp=(disp-dk*dt/2);
da=(dn+dp)/2; % average of d at previous time step and current iteration step

% compute discrete Laplacian % 
lpd(:,1) = -oneconst*((A*da(:,2)).*da(:,3)-(A*da(:,3)).*da(:,2));
lpd(:,2) = -oneconst*((A*da(:,3)).*da(:,1)-(A*da(:,1)).*da(:,3));
lpd(:,3) = -oneconst*((A*da(:,1)).*da(:,2)-(A*da(:,2)).*da(:,1));

% compute source term with electric field
source_t = compute_source(p, t, da,eps1,phiold,phinew,M);

wn = fracm*(fracp*w+dt*(Ml\(lpd+source_t)));
end