% Experiment 5
% Initial data numerical experiments director field coupled to electric
% field
% 2025-06-14
% Experiment 5
% positive signs for eps_1 and eps_2, constant initial vector for d
% different mesh refinements, figure out how many iterations are used on
% average and maximal

 disp = 0.5; % constant in front of d_{tt}
 dk = 0.5; % constant in front of damping term d_t
 oneconst = 1; % constant k, in front of elastic term |\Grad d|^2 in paper
 T = 0.4; % final time of simulation
 
 dx = [1/16,1/32,1/64,1/128]; % grid size
 toll = 0.05*dx.^2; % tolerance in nonlinear iteration (stopping tolerance)
 plotornot = 0; % whether to plot the director field while simulating or not
av_it = zeros(size(dx));
max_it = zeros(size(dx));

 % initial data for director field d
 d1_0 = @(x,y)(sqrt(1/2)*ones(size(y)));
 d2_0 = @(x,y)(sqrt(1/2)*ones(size(y)));
 d3_0 = @(x,y)(zeros(size(y)));
 % initial data for angular momentum w
 w1_0=@(x,y)(zeros(size(x)));
 w2_0=@(x,y)(zeros(size(x)));
 w3_0=@(x,y)(zeros(size(x)));

 gbcfcn = @(t,x,y)(3*sin(2*pi*t+0.2)*(x+0.5).*sin(pi*y)); % boundary condition for elliptic equation (potential)
    
 D=[-0.5,0.5;-0.5,0.5]; % computational domain
 
 eps1 = 5; % constant \eps_1 in paper, which appears in front of source term for w equation
 eps2 = 0.5; % constant in paper which appears in elliptic equation, nonhomogeneous part 

 % create a struct with input parameters and initial data   
 
expname = 'exp5_positive_epsilons_noit';

for j=1:length(dx)
    inputdata = struct('eps1',eps1,'eps2',eps2,'disp',disp,'oneconst', ...
        oneconst,'dk',dk,'T',T,'dx',dx(j),'toll',toll(j),'plotornot', ...
        plotornot,'d1_0',d1_0,'d2_0',d2_0,'d3_0',d3_0,'w1_0',w1_0, ...
        'w2_0',w2_0,'w3_0',w3_0,'D',D,'gbcfcn',gbcfcn);
    [d,dp,df, w,wp,mphi,phip,phif,Ep,Em,Ef,damping,ed,ew,tt,tf,p,t,e,max_it(j),av_it(j)]=wme_fe(inputdata,expname);
end
