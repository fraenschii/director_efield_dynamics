% Initial data numerical experiments director field coupled to electric
% field
% 2021-04-18
% Experiment 1
% swap sign of eps1 and eps2

 disp= 0.5; % constant in front of d_{tt}
 dk = 2; % constant in front of damping term d_t
 oneconst = 1; % constant k, in front of elastic term |\Grad d|^2 in paper
 T = 2; % final time of simulation
 bc = 'N'; % boundary conditions for liquid crystal director, choose either periodic or Neumann ('p' or 'N')

 dx = 1/64; % grid size
 toll = 0.05*dx^2; % tolerance in nonlinear iteration (stopping tolerance)
 plotornot = 0; % whether to plot the director field while simulating or not


 % initial data for director field d
 d1_0 = @(x,y)(sqrt(1/2)*ones(size(y)));
 d2_0 = @(x,y)(sqrt(1/2)*ones(size(y)));
 d3_0 = @(x,y)(zeros(size(y)));
 % initial data for angular momentum w
 w1_0=@(x,y)(zeros(size(x)));
 w2_0=@(x,y)(zeros(size(x)));
 w3_0=@(x,y)(zeros(size(x)));

 gbcfcn = @(t,x,y)(10*sin(2*pi*t+0.2)*(x+0.5).*sin(pi*y)); % boundary condition for elliptic equation (potential)
    
 D=[-0.5,0.5;-0.5,0.5]; % computational domain
 
 eps1 = 5; % constant \eps_1 in paper, which appears in front of source term for w equation
 eps2 = 0.5; % constant in paper which appears in elliptic equation, nonhomogeneous part 

 % create a struct with input parameters and initial data   
 inputdata = struct('eps1',eps1,'eps2',eps2,'disp',disp,'oneconst', ...
        oneconst,'dk',dk,'T',T,'bc',bc,'dx',dx,'toll',toll,'plotornot', ...
        plotornot,'d1_0',d1_0,'d2_0',d2_0,'d3_0',d3_0,'w1_0',w1_0, ...
        'w2_0',w2_0,'w3_0',w3_0,'D',D,'gbcfcn',gbcfcn);
    