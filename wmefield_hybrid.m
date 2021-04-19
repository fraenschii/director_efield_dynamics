function [d,df,w,wf,mphi,phif,Ef,Em,damping,maxgrad,tt,tf]=wmefield_hybrid(inputdata)
% Numerical method for the coupled wave map electric field model for 2D
% spatial domain, but director field can vary in S^2
% 2021-02-25
% This should be a hybrid scheme, pw constants for d and w and pw linears
% for phi, regular grid on a square, triangles for phi with nodes the same
% as corners of cells for d and w
% d_t = w\times d
% disp* w_t=k\Delta d\times d- dk*w + \eps1 (d . E) E\times d.
% E = \Grad \phi, -\Div(\Grad \phi +\eps2 (d . \Grad\phi) d) = 0
% 
% Input variables: inputdata - struct containing input parameters for
% simulation, for an example, see file inputdata1.m

% Output variables:
% d       - solution of wave map eq
% df      - d for movies, with time dependence
% w       - angular momentum \pm d x d_t
% wf      - w for movies with time dependence
% phi     - approximation of electric potential at final time  
% phif    - approximation of electric potential at many time steps (for movie)  
% e_w     - "w"-part of energy
% e_d     - "\Grad d" part of energy
% Em      - e_w+e_d, this part is conserved for dk=0
% et      - 
% Hm      - discrete energy as in paper (without w but instead D_t d)
% maxgrad - maximum of discrete gradient at each time step

if  nargin < 1 % some sample initial data if initdata is not specified
    
    disp= .1; % constant in front of d_{tt}
    dk = 0.15; % constant in front of damping term d_t
    oneconst = 1; % constant k, in front of elastic term |\Grad d|^2 in paper
    T = .2; % final time of simulation
    bc = 'N'; % boundary conditions for liquid crystal director, choose either periodic or Neumann ('p' or 'N')
    dx = 1/64;  % grid size
    toll = 0.1*dx^2; % tolerance in nonlinear iteration (stopping tolerance)
    plotornot = 0; % whether to plot the director field while simulating or not
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
    gbcfcn = @(t,x,y)(5*(x+0.5)*sin(2*pi*t+0.2)); % boundary condition for elliptic equation (potential)
    
    D=[-0.5,0.5;-0.5,0.5]; % computational domain
    eps1 = 0.0; % constant \eps_1 in paper, which appears in front of source term for w equation
    eps2 = 0.2; % constant in paper which appears in elliptic equation, nonhomogeneous part 
else
    eps1 = inputdata.eps1;
    eps2 = inputdata.eps2;
    disp = inputdata.disp;
    oneconst = inputdata.oneconst;
    dk = inputdata.dk;
    T = inputdata.T;
    bc = inputdata.bc;
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

rat=(D(1,2)-D(1,1))/dx/32;

% specify time step, CFL condition, grid etc.
cfl=0.1;
dt=cfl*dx*sqrt(dk*dx^2+disp);
nt=ceil(T/dt);
dt = T/nt;
x=D(1,1)+dx/2:dx:D(1,2)-dx/2;
xphi = D(1,1):dx:D(1,2); % grid for phi interpolation including boundary
[Xphi,Yphi]=meshgrid(xphi);
if bc=='p' %periodic boundary conditions for d
    x_=[x(end),x,x(1)];
    [X_,Y_]=meshgrid(x_);
else
    if bc=='N' % Neumann boundary conditions for d
        x_=[x(1),x,x(end)];
        [X_,Y_]=meshgrid(x_,x_);
    else
        if strcmp(bc,'mixed') % mixed boundary conditions for d
            x_=[x(end),x,x(1)];
            y_=[x(1),x,x(end)]; 
            [X_,Y_]=meshgrid(x_,y_);
        end     
    end
end

% create grid for d and w
nx=length(x);
nxphi = length(xphi);
[X,Y]=meshgrid(x);
[Xrat,Yrat]=meshgrid(x(rat:rat:end));

% initialize energy and damping variables
Em=zeros(1,nt+1);
damping = zeros(1,nt+1);
maxgrad = zeros(1,nt+1);

% initialize d and w variables and set initial data
d=zeros(nx+2,nx+2,3);
w=zeros(nx,nx,3);

d(:,:,1)=d1_0(X_,Y_);
d(:,:,2)=d2_0(X_,Y_);
d(:,:,3)=d3_0(X_,Y_);
w(:,:,1)=w1_0(X,Y);
w(:,:,2)=w2_0(X,Y);
w(:,:,3)=w3_0(X,Y);
% initialize boundary data for phi
gbc = gbcfcn(0,Xphi,Yphi);

% compute initial energies
Em(1)=disp*dx^2*sum(sum(sum(w.^2)))+oneconst*sum(sum(sum((diff(d(1:end-1,2:end-1,:),1,1).^2)+(diff(d(2:end-1,1:end-1,:),1,2).^2))));
d_deltax=diff(d(1:end-1,2:end-1,1:3),1,1)/dx;
d_deltay=diff(d(2:end-1,1:end-1,1:3),1,2)/dx;
maxgrad(1)=max(max(sqrt(abs(d_deltax(:,:,1)).^2+abs(d_deltax(:,:,2)).^2+...
            abs(d_deltay(:,:,1)).^2+abs(d_deltay(:,:,2)).^2+...
            abs(d_deltax(:,:,3)).^2+abs(d_deltay(:,:,3)).^2)));

cumt=0; % current time in simulation
if plotornot % plot director field in quiver plot live
    figure;
    hq=quiver(X(rat:rat:end,rat:rat:end),Y(rat:rat:end,rat:rat:end),...
        d(rat+1:rat:end-1,1+rat:rat:end-1,1),d(1+rat:rat:end-1,1+rat:rat:end-1,2));
    axis image;
    htitle=get(gca,'Title');
    drawnow;
    plotrate=90;
end

% space time data collection for movies
movierate = 10;
NT=ceil(nt/movierate);
ncx=nx/rat;
df=zeros(ncx,ncx,3,NT); % director field
wf=zeros(ncx,ncx,3,NT); % angular momentum
phif=zeros(nxphi,nxphi,NT); % electric potential
Ef=zeros(nxphi-1,nxphi-1,2,NT); % electric field

tf = zeros(1,NT); % variable for time steps
% set initial data
df(:,:,:,1)=d(rat+1:rat:end-1,1+rat:rat:end-1,:);
wf(:,:,:,1)=w(rat:rat:end,rat:rat:end,:);
phiold = solve_elliptic_nosource(d(:,:,1),d(:,:,2),gbc,eps2); % initial approximation of phi
phif(:,:,1)=phiold;
Ef(:,:,1,1) = diff(phiold(:,2:end),1,1)/dx;
Ef(:,:,2,1) = diff(phiold(2:end,:),1,2)/dx;
       
countt=2; % number of time steps executed
for n=1:nt
    cumt=cumt+dt; % in crease current time
    gbc = gbcfcn(cumt,Xphi,Yphi); % set boundary condition
    % solve for next time step d^{m+1},w^{m+1}, phi^{m+1} using fixed pt
    % iteration
    [dn,wn,phin]=nonit_2d(d,w,phiold,dx,dt,toll,disp,dk,eps1,bc,eps2,oneconst,gbc);
  
    damping(n+1)=damping(n) ...
        +2*dk*dt*dx^2*sum(sum(sum((w+wn).^2/4))); % total damping at time t^n
      
    d=dn; %update d
    w=wn; % update w
    phiold = phin;
    
    % compute energy
    e_w=disp*dx^2*sum(sum(sum(w.^2))); % contribution from sigma |w|^2 term to energy
    e_d=oneconst*sum(sum(sum((diff(d(1:end-1,2:end-1,:),1,1).^2) ...
        +(diff(d(2:end-1,1:end-1,:),1,2).^2)))); % contribution of elastic term |\Grad d|^2 to energy
    Em(n+1) = e_w+e_d; % total energy \int \disp |w|^2+ k |\Grad d|^2 dx at t^{n}
    % compute maximum discrete gradient
    d_deltax=diff(d(1:end-1,2:end-1,1:3),1,1)/dx;
    d_deltay=diff(d(2:end-1,1:end-1,1:3),1,2)/dx;
    maxgrad(n+1)=max(max(sqrt(abs(d_deltax(:,:,1)).^2+abs(d_deltax(:,:,2)).^2 ...
        +abs(d_deltay(:,:,1)).^2+abs(d_deltay(:,:,2)).^2 ...
        +abs(d_deltax(:,:,3)).^2+abs(d_deltay(:,:,3)).^2)));  % maximum gradient (not checked)      
      
    if (plotornot && mod(n,plotrate)==0) % update quiver plot for director field
        set(hq,'Udata',d(rat+1:rat:end-1,1+rat:rat:end-1,1),...
            'Vdata',d(1+rat:rat:end-1,1+rat:rat:end-1,2));
             %quiver(X(rat:rat:end,rat:rat:end),Y(rat:rat:end,rat:rat:end),...
             %            d_(rat+1:rat:end-1,1+rat:rat:end-1,1),d_(1+rat:rat:end-1,1+rat:rat:end-1,2),'.');
        set(htitle,'String', sprintf(' timestep no. %d, T = %7.3f ',n,cumt));
        drawnow;
    end
    
    % every 'movierate'th time step, save approximation for d,w and phi in
    % df, wf,phif and Ef (electric field)
    if mod(n,movierate)==0
        df(:,:,:,countt)=d(rat+1:rat:end-1,1+rat:rat:end-1,:);
        wf(:,:,:,countt)=w(rat:rat:end,rat:rat:end,:);
        phif(:,:,countt)=phin; 
        Ef(:,:,1,countt) = diff(phin(:,2:end),1,1)/dx;
        Ef(:,:,2,countt) = diff(phin(2:end,:),1,2)/dx;
       
        % save time step
        countt=countt+1;
        tf(countt)=cumt;
    end
end

mphi = phin; 
Em = Em/2;
damping = damping/2;
% plot director field
figure;quiver(Xrat,Yrat,d(rat+1:rat:end-1,rat+1:rat:end-1,1), ...
    d(rat+1:rat:end-1,rat+1:rat:end-1,2)); 
axis image;
ylabel('$y$','rotation',0,'FontSize',16,'Interpreter','LateX')
xlabel('$x$','FontSize',16,'Interpreter','LateX')
ax = gca;
ax.FontSize = 12; 
title(['T=',num2str(T)])
% plot electric field
figure;quiver(Xrat,Yrat,Ef(rat:rat:end,rat:rat:end,1,end), ...
    Ef(rat:rat:end,rat:rat:end,2,end));
axis image;
ylabel('$y$','rotation',0,'FontSize',16,'Interpreter','LateX')
xlabel('$x$','FontSize',16,'Interpreter','LateX')
ax = gca;
ax.FontSize = 12; 
title(['T=',num2str(T)])

% plot energy and damping
tt=0:dt:T;
figure; plot(tt,Em,'-',tt,damping,'--','LineWidth',2);
legend('$E_m$', 'damping','FontSize',12,'Interpreter','LateX');
% legend('$E_m$','H_m', 'damping');
xlabel('$t$','FontSize',16,'Interpreter','LateX');
ax = gca;
ax.FontSize = 12; 
end

%% fixed point iteration
function [dnew,wnew,phinew]=nonit_2d(dold,wold,phiold,dx,dt,toll,disp,dk,eps1,bc,eps2,oneconst,bcdata) 
% fixed point iteration to solve nonlinear algebraic system in 2D, 1 time
% step of length dt

nit=1;
% first step of iteration (forward Euler kind of)
dnew = compute_d_2d(dold(2:end-1,2:end-1,:),wold,dt,bc); 
phinew = solve_elliptic_nosource(dnew(:,:,1),dnew(:,:,2),bcdata,eps2);
wnew=compute_w_2d(dnew,dold,wold,phiold,phinew,dt,dx,dk,disp,eps1,oneconst);

% compute the initial error
% square of error
err= sum(sum(sum(dx^2*(wnew-wold).^2+...
      (diff(dnew(1:end-1,2:end-1,:)-dold(1:end-1,2:end-1,:),1,1).^2)...
      +(diff(dnew(2:end-1,1:end-1,:)-dold(2:end-1,1:end-1,:),1,2).^2) ))) ...
      +sum(sum(diff(phinew-phiold,1,1).^2))+sum(sum(diff(phinew-phiold,1,2).^2)); 

% do the fixed point iteration     
while(nit<400 && err>toll^2)  
    dmid=dnew;
    wmid=wnew;
    phimid = phinew;
    dnew =compute_d_2d(dold(2:end-1,2:end-1,:),(wold+wmid)/2,dt,bc); % compute d^{m,s+1}
    phinew = solve_elliptic_nosource(dnew(:,:,1),dnew(:,:,2),bcdata,eps2); % compute phi^{m,s+1}
    wnew=compute_w_2d(dnew,dold,wold,phiold,phinew,dt,dx,dk,disp,eps1,oneconst); % compute w^{m,s+1}
    nit=nit+1;
    
    % compute the error in the nth iteration
    err= sum(sum(sum(dx^2*(wnew-wmid).^2+...
      (diff(dnew(1:end-1,2:end-1,:)-dmid(1:end-1,2:end-1,:),1,1).^2)...
      +(diff(dnew(2:end-1,1:end-1,:)-dmid(2:end-1,1:end-1,:),1,2).^2) ))) ...
      +sum(sum(diff(phinew-phimid,1,1).^2))+sum(sum(diff(phinew-phimid,1,2).^2));  
end

resterr=sum(sum(sum(dx^2*(wnew-wmid).^2+...
      (diff(dnew(1:end-1,2:end-1,:)-dmid(1:end-1,2:end-1,:),1,1).^2)...
      +(diff(dnew(2:end-1,1:end-1,:)-dmid(2:end-1,1:end-1,:),1,2).^2) )));
phierr =  sum(sum(diff(phinew-phimid,1,1).^2))+sum(sum(diff(phinew-phimid,1,2).^2));

% if more than 399 iterations are needed, print out the error and stop
% iteration
if nit>=399
        fprintf('No convergence\n');
        fprintf(['error=',num2str(err),'\n']);
        fprintf(['phierror=',num2str(phierr),'\n']);
        fprintf(['resterror=',num2str(resterr),'\n']);       
end     
end

function v=bc_2d(v, bc) 
% 2D boundary conditions for d 
    switch bc
        case 'N'
            v([1 end],2:end-1,:)=v([2 end-1],2:end-1,:);
            v(2:end-1,[1 end],:)=v(2:end-1,[2 end-1],:);
        case 'p'
            v([1 end],2:end-1,:)=v([end-1 2],2:end-1,:);
            v(2:end-1,[1 end],:)=v(2:end-1,[end-1 2],:);
    end
end

%%  evolution of d
function dn=compute_d_2d(d,w,dt,bc) 
% one step evolution of d 
dn=zeros(size(d,1)+2,size(d,2)+2,size(d,3));
dtm =dt/2;
dtm2 = dtm^2;
c1=d(:,:,1)+dtm* (d(:,:,2).*w(:,:,3)-d(:,:,3).*w(:,:,2));
c2=d(:,:,2)+dtm* (d(:,:,3).*w(:,:,1)-d(:,:,1).*w(:,:,3));
c3=d(:,:,3)+dtm* (d(:,:,1).*w(:,:,2)-d(:,:,2).*w(:,:,1));

a11=1+dtm2*(w(:,:,1).^2);
a12=dtm2*w(:,:,2).*w(:,:,1)+dtm*w(:,:,3);
a13=dtm2*w(:,:,1).*w(:,:,3)-dtm*w(:,:,2);
a21=dtm2*w(:,:,1).*w(:,:,2)-dtm*w(:,:,3);
a22=1+dtm2*(w(:,:,2).^2);
a23=dtm2*w(:,:,2).*w(:,:,3)+dtm*w(:,:,1);
a31=dtm2*w(:,:,1).*w(:,:,3)+dtm*w(:,:,2);
a32=dtm2*w(:,:,2).*w(:,:,3)-dtm*w(:,:,1);
a33=1+dtm2*(w(:,:,3).^2);
detn=1+dtm2*((w(:,:,1).^2)+(w(:,:,2).^2)+(w(:,:,3).^2));
dn(2:end-1,2:end-1,1)=(a11.*c1+a12.*c2+a13.*c3)./detn;
dn(2:end-1,2:end-1,2)=(a21.*c1+a22.*c2+a23.*c3)./detn;
dn(2:end-1,2:end-1,3)=(a31.*c1+a32.*c2+a33.*c3)./detn;
dn = bc_2d(dn,bc);
end

%% w evolution
function wn=compute_w_2d(dn,dp,w,phiold,phinew,dt,dx,dk,disp,eps1,oneconst)
%dn, dp with boundary data
fracm=1/(disp+dk*dt/2);
fracp=(disp-dk*dt/2);
wn=zeros(size(w));
da=(dn+dp)/2;
dd=da(2:end-1,2:end-1,:);
lpd=oneconst*(da(3:end,2:end-1,:)+da(1:end-2,2:end-1,:)+da(2:end-1,3:end,:)+...
    da(2:end-1,1:end-2,:)-4*dd)/dx^2; % discrete 5pt stencil laplacian
% compute source term with electric field
source_t = compute_source_weq_hybrid(da(:,:,1),da(:,:,2),da(:,:,3),phiold,phinew,dx,eps1);
% update w
wn(:,:,1)=fracm*(fracp*w(:,:,1)+dt*(lpd(:,:,2).*dd(:,:,3) ... 
    - lpd(:,:,3).*dd(:,:,2) + source_t(:,:,1)));
wn(:,:,2)=fracm*(fracp*w(:,:,2)+dt*(lpd(:,:,3).*dd(:,:,1) ...
    - lpd(:,:,1).*dd(:,:,3)  + source_t(:,:,2)));
wn(:,:,3)=fracm*(fracp*w(:,:,3)+dt*(lpd(:,:,1).*dd(:,:,2) ...
    - lpd(:,:,2).*dd(:,:,1)  + source_t(:,:,3)));
end

function source = compute_source_weq_hybrid(dav1,dav2,dav3,phiold,phinew,dx,eps1)
% compute source term with electric field for w equation
% dav1, dav2 1st and 2nd component of averaged d, phiold old values of phi,
% phinew next time step values of phi, dx grid size, eps1 parameter from
% paper

[Ntot,~] = size(dav1);
N = Ntot-2;
source = zeros(N,N,3);

comp1 = eps1*(dav1(2:end-1,2:end-1).*( ...
    (phiold(2:end,2:end)-phiold(1:end-1,2:end)).* ...
    (phinew(2:end,2:end)-phinew(1:end-1,2:end)) ...
    +(phiold(2:end,1:end-1)-phiold(1:end-1,1:end-1)).* ...
    (phinew(2:end,1:end-1)-phinew(1:end-1,1:end-1))) + ...
    dav2(2:end-1,2:end-1).* ...
    ((phiold(1:end-1,2:end)-phiold(1:end-1,1:end-1)).* ...
    (phinew(2:end,2:end)-phinew(1:end-1,2:end)) + ... 
    (phiold(2:end,2:end)-phiold(2:end,1:end-1)).* ...
    (phinew(2:end,1:end-1)-phinew(1:end-1,1:end-1))))/4/dx^2;

comp2 = eps1*( dav1(2:end-1,2:end-1).* ( ...
    (phiold(2:end,2:end)-phiold(1:end-1,2:end)).* ...
    (phinew(1:end-1,2:end)-phinew(1:end-1,1:end-1)) + ...
    (phiold(2:end,1:end-1)-phiold(1:end-1,1:end-1)).* ...
    (phinew(2:end,2:end)-phinew(2:end,1:end-1))) ...
    +dav2(2:end-1,2:end-1).* ( ...
    (phiold(1:end-1,2:end)-phiold(1:end-1,1:end-1)).* ...
    (phinew(1:end-1,2:end)-phinew(1:end-1,1:end-1)) + ...
    (phiold(2:end,2:end)-phiold(2:end,1:end-1)).* ...
    (phinew(2:end,2:end)-phinew(2:end,1:end-1))))/4/dx^2;

source(:,:,1) = comp2.*dav3(2:end-1,2:end-1);
source(:,:,2) = -comp1.*dav3(2:end-1,2:end-1);
source(:,:,3) = comp1.*dav2(2:end-1,2:end-1)-comp2.*dav1(2:end-1,2:end-1);
end

%% solve the elliptic equation
function A = stiffness_matrix(d1,d2,eps2)
% create the stiffness matrix for the equation for phi
% d1, d2 are first and second component of director field, eps2 the 
% constant \epsi_2

[N,~]= size(d1);
Nint = N-3;

% case 1: central diagonal
diag0_square = 4+ eps2/2*( d1(2:end-2,2:end-2).^2+d2(2:end-2,2:end-2).^2 ...
    +d1(3:end-1,3:end-1).^2+d2(3:end-1,3:end-1).^2 + ...
    (d1(2:end-2,3:end-1)-d2(2:end-2,3:end-1)).^2 + ...
    (d1(3:end-1,2:end-2)-d2(3:end-1,2:end-2)).^2);
diag0 = reshape(diag0_square,[],1);

% case 2a: first upper diagonal
diagp1_square = -1 +eps2/2*(d1(3:end-2,2:end-2).*(d2(3:end-2,2:end-2)-d1(3:end-2,2:end-2)) ...
    + d1(3:end-2,3:end-1).*(d2(3:end-2,3:end-1)-d1(3:end-2,3:end-1)));
diagp1_square = [diagp1_square;zeros(1,Nint)];
diagp1 = reshape(diagp1_square,[],1);

% case 2b: first lower diagonal
diagm1_square = -1 +eps2/2*(d1(3:end-2,3:end-1).*(d2(3:end-2,3:end-1)-d1(3:end-2,3:end-1)) ...
    + d1(3:end-2,2:end-2).*(d2(3:end-2,2:end-2)-d1(3:end-2,2:end-2)));
diagm1_square = [zeros(1,Nint);diagm1_square];
diagm1 = reshape(diagm1_square,[],1);

% Case 3a: Nth upper diagonal
diagpN_square = -1+ eps2/2*(d2(2:end-2,3:end-2).*(d1(2:end-2,3:end-2)-d2(2:end-2,3:end-2)) ...
    + d2(3:end-1,3:end-2).*(d1(3:end-1,3:end-2)-d2(3:end-1,3:end-2)));
diagpN_square = [diagpN_square,zeros(Nint,1)];
diagpN = reshape(diagpN_square,[],1);

% case 3b: Nth lower diagonal
diagmN_square = -1 +eps2/2*(d2(2:end-2,3:end-2).*(d1(2:end-2,3:end-2)-d2(2:end-2,3:end-2)) ...
    + d2(3:end-1,3:end-2).*(d1(3:end-1,3:end-2)-d2(3:end-1,3:end-2)));
diagmN_square = [zeros(Nint,1),diagmN_square];
diagmN = reshape(diagmN_square,[],1);

% Case 4:N+1th upper diagonal
diagpN1_square = -eps2*d1(3:end-2,3:end-2).*d2(3:end-2,3:end-2);
diagpN1_square = [diagpN1_square,zeros(Nint-1,1);zeros(1,Nint)];
diagpN1 = reshape(diagpN1_square,[],1);

% Case 5: N+1th lower diagonal
diagmN1_square =  -eps2*d1(3:end-2,3:end-2).*d2(3:end-2,3:end-2);
diagmN1_square = [zeros(1,Nint);zeros(Nint-1,1),diagmN1_square];
diagmN1 = reshape(diagmN1_square,[],1);

coeffs = [diagpN1,diagpN,diagp1,diag0,diagm1,diagmN,diagmN1];

A = spdiags(coeffs,[-Nint-1,-Nint,-1,0,1,Nint,Nint+1],Nint^2,Nint^2);
A = A'; % because of the way Matlab creates sparse matrices...
end


function F = right_hand_side_nosource(gbc,d1,d2,eps2)
% compute right hand side for elliptic equation without source term f
% gbc boundary data, assume it is a N+1 x N+1 matrix, can be zeros in
% interior points
% d1, d2 coefficients from 1st and 2nd component of d, eps2 the
% coefficient \eps_2 from the paper

% case 1: central diagonal
diag0_square = 4+ eps2/2*( d1(2:end-2,2:end-2).^2+d2(2:end-2,2:end-2).^2 ...
    +d1(3:end-1,3:end-1).^2+d2(3:end-1,3:end-1).^2 + ...
    (d1(2:end-2,3:end-1)-d2(2:end-2,3:end-1)).^2 + ...
    (d1(3:end-1,2:end-2)-d2(3:end-1,2:end-2)).^2); 
centercoeff = diag0_square.*gbc(2:end-1,2:end-1); 

% case 2a: 1st upper diagonal
diagp1_square = -1 +eps2/2*(d1(3:end-1,2:end-2).*(d2(3:end-1,2:end-2)-d1(3:end-1,2:end-2)) ...
    + d1(3:end-1,3:end-1).*(d2(3:end-1,3:end-1)-d1(3:end-1,3:end-1)));
coeffp1 = diagp1_square.*gbc(3:end,2:end-1);

% case 2b: 1st lower diagonal
diagm1_square = -1 +eps2/2*(d1(2:end-2,3:end-1).*(d2(2:end-2,3:end-1)-d1(2:end-2,3:end-1)) ...
    + d1(2:end-2,2:end-2).*(d2(2:end-2,2:end-2)-d1(2:end-2,2:end-2)));
coeffm1 = diagm1_square.*gbc(1:end-2,2:end-1);

% Case 3a: Nth upper diagonal
diagpN_square = -1+ eps2/2*(d2(2:end-2,3:end-1).*(d1(2:end-2,3:end-1)-d2(2:end-2,3:end-1)) ...
    + d2(3:end-1,3:end-1).*(d1(3:end-1,3:end-1)-d2(3:end-1,3:end-1)));
coeffpN = diagpN_square.*gbc(2:end-1,3:end);

% Case 3b: Nth lower diagonal
diagmN_square = -1 +eps2/2*(d2(2:end-2,2:end-2).*(d1(2:end-2,2:end-2)-d2(2:end-2,2:end-2)) ...
    + d2(3:end-1,2:end-2).*(d1(3:end-1,2:end-2)-d2(3:end-1,2:end-2)));
coeffmN = diagmN_square.*gbc(2:end-1,1:end-2);

% Case 4: N+1th upper diagonal
diagpN1_square = -eps2*d1(3:end-1,3:end-1).*d2(3:end-1,3:end-1);
coeffpN1 = diagpN1_square.*gbc(3:end,3:end);

% Case 5: N+1th lower diagonal
diagmN1_square =  -eps2*d1(2:end-2,2:end-2).*d2(2:end-2,2:end-2);
coeffmN1 = diagmN1_square.*gbc(1:end-2,1:end-2);

rhssquare = centercoeff + coeffp1 + coeffm1 + coeffpN + coeffmN + coeffmN1 ...
    + coeffpN1;
F = -reshape(rhssquare,[],1);
end

function phinew = solve_elliptic_nosource(d1,d2,gbc,eps2)
% solve the elliptic equation for phi with given coefficients and director
% field d, right hand side only involves boundary condtions (no free charges)
    [N,~]= size(gbc);
    Nint = N-2;
    F = right_hand_side_nosource(gbc,d1,d2,eps2);
    A = stiffness_matrix(d1,d2,eps2);
    
    phivec = A\F;
    phinew = gbc;
    phinew(2:end-1,2:end-1) =gbc(2:end-1,2:end-1)+ reshape(phivec,Nint,Nint);
end