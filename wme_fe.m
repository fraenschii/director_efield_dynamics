function [d,dp,df, w,wp,mphi,phip,phif,Ep,Em,Ef,damping,ed,ew,tt,tf,p,t,e,max_it,av_it]=wme_fe(inputdata,expname)
% Finite element method for the coupled wave map electric field model for 2D
% spatial domain, but director field can vary in S^2
% 2025-03-23
% piecewise linears for d and w and phi, regular grid on a square
% Equations: 
% d_t = w\times d
% disp* w_t=k\Delta d\times d- dk*w + \eps1 (d . E) E\times d.
% E = \Grad \phi, -\Div(\Grad \phi +\eps2 (d . \Grad\phi) d) = 0
% 
% Input variables: inputdata - struct containing input parameters for
% simulation, 
% for an example, see file inputdata1.m, expname is a string to name the
% experiment and the output plots

% Output variables: 
% d       - director field approximation at nodes
% dp      - d at different time steps at nodes
% df      - this is dp interpolated on a grid  
% w       - angular momentum d_t x d
% wp      - w at different time steps at nodes
% mphi    - approximation of electric potential at final time  
% phip    - approximation of electric potential at many time steps 
% phif    - phip interpolated on a grid
% Ep      - electric field (\nabla \phi) at different time steps  
% Ef      - approximation of the electric field, interpolation on a regular grid
% Em      - approximation of reduced energy, evolution over time
% damping - approximation of the damping term of the energy
% ed      - elastic part of the energy
% ew      - kinetic part of the energy
% tt      - vector with all time steps at which approximation was computed
% tf      - vector with time steps coresponding to .f and .p variables
% p       - coordinates of nodes of mesh
% t       - triangles of the mesh with their nodes
% e       - boundary nodes of mesh
% max_it  - variable to store the maximum number of fixed point iterations
% used
% av_it   - variable to store the average number of fixed point iterations neeeded over all time steps 





if  nargin < 1 % some sample initial data if initdata is not specified
    expname = 'unspecified_name';
    disp= 1; % constant in front of d_{tt}
    dk = 0.0;   % constant in front of damping term d_t
    oneconst = 1; % constant k, in front of elastic term |\Grad d|^2 in paper
    T = 2.8; % final time of simulation
    
    dx = 1/32;  % grid size
    toll = 0.1*dx^2; % tolerance in nonlinear iteration (stopping tolerance)
    plotornot = 1; % whether to plot the director field while simulating or not
    % initial data for director field
    af=@(r)((1-2*r).^4);
    rf=@(x,y)(sqrt(x.^2+y.^2));
    denominatorf=@(x,y)((rf(x,y).^2+af(rf(x,y)).^2));
    d1_0=@(x,y)(2*(rf(x,y)<=0.5).*(x).*af(rf(x,y))./(denominatorf(x,y)));
    d2_0=@(x,y)(2*(rf(x,y)<=0.5).*(y).*af(rf(x,y))./(denominatorf(x,y)));
    d3_0=@(x,y)((rf(x,y)<=0.5).*(af(rf(x,y)).^2-rf(x,y).^2)./(denominatorf(x,y))-(rf(x,y)>0.5));

    % initial data for angular momentum
    w1_0=@(x,y)(zeros(size(x)));
    w2_0=@(x,y)(zeros(size(x)));
    w3_0=@(x,y)(zeros(size(x)));
    gbcfcn = @(t,x,y)(3*(x+0.5)*sin(2*pi*t+0.2)); % boundary condition for elliptic equation (potential)

    D=[-0.5,0.5;-0.5,0.5]; % computational domain
    eps1 = 5; % constant \eps_1 in paper, which appears in front of source term for w equation
    eps2 = 0.2; % constant in paper which appears in elliptic equation, nonhomogeneous part 
else
    eps1 = inputdata.eps1;
    eps2 = inputdata.eps2;
    disp = inputdata.disp;
    oneconst = inputdata.oneconst;
    dk = inputdata.dk;
    T = inputdata.T;
    dx = inputdata.dx;
    toll = inputdata.toll;
    plotornot = inputdata.plotornot;
    d1_0 = inputdata.d1_0;
    d2_0 = inputdata.d2_0;
    d3_0 = inputdata.d3_0;
    w1_0 = inputdata.w1_0;
    w2_0 = inputdata.w2_0;
    w3_0 = inputdata.w3_0;
    D = inputdata.D;
    gbcfcn = inputdata.gbcfcn;
    
   
end

%% specify time step, CFL condition, grid etc.

Nx = ceil((D(1,2)-D(1,1))/dx);

dx = (D(1,2)-D(1,1))/Nx;
cfl=0.4;
dt=cfl*dx*sqrt(dk*dx^2+disp);
nt=ceil(T/dt);
dt = T/nt;
[p, t, e] = generate_mesh(Nx,D); % generates a uniform mesh on a square
boundary_nodes = extract_boundary_nodes(e);
boundary_node_coords = p(boundary_nodes, :);
Np = size(p,1); % number of nodes

av_it = 0; % variable for the average number of iterations used
max_it = 0; % variable for the maximum number of iterations used

% initialize energy and damping variables
Em=zeros(1,nt+1);
ew=zeros(1,nt);
ed=zeros(1,nt);
damping = zeros(1,nt+1);
 
% initialize d and w variables and set initial data
d=zeros(Np,3);
gbc = zeros(Np,1);

w=zeros(Np,3);
sizep = size(p);
d(:,1)=d1_0(p(:,1),p(:,2));
d(:,2)=d2_0(p(:,1),p(:,2));
d(:,3)=d3_0(p(:,1),p(:,2));
w(:,1)=w1_0(p(:,1),p(:,2));
w(:,2)=w2_0(p(:,1),p(:,2));
w(:,3)=w3_0(p(:,1),p(:,2));

% initialize boundary data for phi
gbc(boundary_nodes,1) = gbcfcn(0,boundary_node_coords(:,1),boundary_node_coords(:,2));
M = assemble_mass_matrix(p, t); % assemble mass matrix
A = assemble_standard_stiffness_matrix(p,t); % assemble stiffness matrix
Ml = assemble_lumped_mass_matrix(p, t); % mass lumped matrix

% compute initial energies
Em(1) = disp* (w(:,1)'*Ml*w(:,1)+w(:,2)'*Ml*w(:,2)+w(:,3)'*Ml*w(:,3)) ...
    +oneconst*(d(:,1)'*A*d(:,1)+d(:,2)'*A*d(:,2)+d(:,3)'*A*d(:,3)); 


cumt=0; % current time in simulation
h_plot =1/40;
[X,Y] = meshgrid(D(1,1):h_plot:D(1,2)); % grid for plots
    
if plotornot % plot director field in quiver plot live
    figure;
    d1 = griddata(p(:,1), p(:,2), d(:,1), X, Y);
    d2 = griddata(p(:,1), p(:,2), d(:,2), X, Y);
    hq = quiver(X,Y,d1,d2);
    axis image;
    htitle=get(gca,'Title');
    drawnow;
    plotrate=90;
    
end

% space time data collection for movies
movierate = 10;
[NX1,NX2]= size(X);
NT=ceil(nt/movierate);
 
% variables at nodes of mesh
dp=zeros(sizep(1),3,NT); % director field
wp=zeros(sizep(1),3,NT); % angular momentum
phip=zeros(sizep(1),NT); % electric potential
Ep=zeros(sizep(1),2,NT); % electric field
dp(:,:,1) = d;
wp(:,:,1)= w;

% to interpolate variables in 2d regular grid
df=zeros(NX1,NX2,3,NT); % director field
wf=zeros(NX1,NX2,3,NT); % angular momentum
phif=zeros(NX1,NX2,NT); % electric potential
Ef=zeros(NX1,NX2,2,NT); % electric field

tf = zeros(1,NT); % variable for time steps
% set initial data, interpolate on grid X,Y
df(:,:,1,1)=griddata(p(:,1), p(:,2), d(:,1), X, Y);
df(:,:,2,1)=griddata(p(:,1), p(:,2), d(:,2), X, Y);
df(:,:,3,1)=griddata(p(:,1), p(:,2), d(:,3), X, Y);
wf(:,:,1,1)=griddata(p(:,1), p(:,2), w(:,1), X, Y);
wf(:,:,2,1)=griddata(p(:,1), p(:,2), w(:,2), X, Y);
wf(:,:,3,1)=griddata(p(:,1), p(:,2), w(:,3), X, Y);

% find initial approximation of phi
phiold = solve_elliptic_nosource(d,gbc,eps2,p,t,e,M); % initial approximation of phi
phif(:,:,1)=griddata(p(:,1), p(:,2), phiold, X, Y);
phip(:,1) = phiold;
grad_phif = evaluate_fem_gradient(p, t, phiold); % the approximation of the electric field E =\nabla \phi
Ep(:,:,1) = grad_phif;
Ef(:,:,1,1) = griddata(p(:,1), p(:,2), grad_phif(:,1), X, Y); 
Ef(:,:,2,1) = griddata(p(:,1), p(:,2), grad_phif(:,2), X, Y); 

countt=2; % number of time steps executed
for n=1:nt
    cumt=cumt+dt; % in crease current time
    % set boundary data
    gbc(boundary_nodes,1) = gbcfcn(cumt,boundary_node_coords(:,1),boundary_node_coords(:,2));
    
    % solve for next time step d^{m+1},w^{m+1}, phi^{m+1} using fixed pt
    % iteration
    [dn,wn,phin,no_it]=nonit_2d(d,w,phiold,dt,toll,disp,dk,eps1,gbc,eps2,oneconst,p,t,e,M,A,Ml);
    av_it = av_it + no_it; % add number of iterations needed in this time step
    max_it = max(max_it,no_it); % maximum number of iterations needed is updated
    w_av = (w+wn)/2;  
    
    % update energy dissipation
    damping(n+1)=damping(n)+2*dk*dt*(w_av(:,1)'*Ml*w_av(:,1) ...
        +w_av(:,2)'*Ml*w_av(:,2)+w_av(:,3)'*Ml*w_av(:,3)); 
    
    d=dn; %update d
    w=wn; % update w
    phiold = phin;
    
    % compute energy 
    e_w=disp*(w(:,1)'*Ml*w(:,1)+w(:,2)'*Ml*w(:,2)+w(:,3)'*Ml*w(:,3)) ; % contribution from sigma |w|^2 term to energy
    e_d=oneconst*(d(:,1)'*A*d(:,1)+d(:,2)'*A*d(:,2)+d(:,3)'*A*d(:,3)); % contribution of elastic term |\Grad d|^2 to energy
    ew(n) = e_w;
    ed(n) = e_d;
    Em(n+1) = e_w+e_d; % total energy \int \disp |w|^2+ k |\Grad d|^2 dx at t^{n}
   
    if (plotornot && mod(n,plotrate)==0) % update quiver plot for director field
        d1=griddata(p(:,1), p(:,2), d(:,1), X, Y);
        d2=griddata(p(:,1), p(:,2), d(:,2), X, Y);
        set(hq,'Udata',d1,'Vdata',d2);
            
        set(htitle,'String', sprintf(' timestep no. %d, T = %7.3f ',n,cumt));
        drawnow;
    end
    
    % every 'movierate'th time step, save approximation for d,w and phi in
    % df, wf,phif and Ef (electric field)
    if mod(n,movierate)==0
        
        dp(:,:,countt)=d;
        wp(:,:,countt)=w;
        phip(:,countt)=phin;
        df(:,:,1,countt)=griddata(p(:,1), p(:,2), d(:,1), X, Y);
        df(:,:,2,countt)=griddata(p(:,1), p(:,2), d(:,2), X, Y);
        df(:,:,3,countt)=griddata(p(:,1), p(:,2), d(:,3), X, Y);
        wf(:,:,1,countt)=griddata(p(:,1), p(:,2), w(:,1), X, Y);
        wf(:,:,2,countt)=griddata(p(:,1), p(:,2), w(:,2), X, Y);
        wf(:,:,3,countt)=griddata(p(:,1), p(:,2), w(:,3), X, Y);
        phif(:,:,countt)=griddata(p(:,1), p(:,2), phin, X, Y);
        grad_phif = evaluate_fem_gradient(p, t, phin);
        Ep(:,:,countt) = grad_phif;
        Ef(:,:,1,countt) = griddata(p(:,1), p(:,2), grad_phif(:,1), X, Y); 
        Ef(:,:,2,countt) = griddata(p(:,1), p(:,2), grad_phif(:,2), X, Y);
        
        % save time step
        countt=countt+1;
        tf(countt)=cumt;
    end
end
av_it = av_it/nt;

mphi = phin; 
Em = Em/2;
damping = damping/2;
% plot director field
figure;quiver(X,Y,df(:,:,1,end),df(:,:,2,end)); 
axis image;
ylabel('$y$','rotation',0,'FontSize',16,'Interpreter','LateX')
xlabel('$x$','FontSize',16,'Interpreter','LateX')
ax = gca;
ax.FontSize = 12; 
title(['director field,T=',num2str(T)])
filename1 = [expname,num2str(T),'_d.fig'];
savefig(filename1);
close(gcf);

% plot electric field
figure;quiver(X,Y,Ef(:,:,1,end), Ef(:,:,2,end));
axis image;
ylabel('$y$','rotation',0,'FontSize',16,'Interpreter','LateX')
xlabel('$x$','FontSize',16,'Interpreter','LateX')
ax = gca;
ax.FontSize = 12; 
title(['electric field,T=',num2str(T)])
filename2 = [expname,num2str(T),'_E.fig'];
savefig(filename2);
close(gcf);

% plot energy and damping
tt=0:dt:T;
figure; plot(tt,Em,'-',tt,damping,'--','LineWidth',2);
legend('$E_m$', 'damping','FontSize',12,'Interpreter','LateX');
xlabel('$t$','FontSize',16,'Interpreter','LateX');
ax = gca;
ax.FontSize = 12; 
filename3 = [expname,'_energy.fig'];
savefig(filename3);
close(gcf);

end