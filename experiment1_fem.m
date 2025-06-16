% Experiment 1
% Initial data numerical experiments director field coupled to electric
% field
% 2025-06-14
% Experiment 1
% positive signs for eps_1 and eps_2, constant initial vector for d

 disp = 0.5; % constant in front of d_{tt}
 dk = 2; % constant in front of damping term d_t
 oneconst = 1; % constant k, in front of elastic term |\Grad d|^2 in paper
 T = 2; % final time of simulation
 
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

 gbcfcn = @(t,x,y)(3*sin(2*pi*t+0.2)*(x+0.5).*sin(pi*y)); % boundary condition for elliptic equation (potential)
    
 D=[-0.5,0.5;-0.5,0.5]; % computational domain
 
 eps1 = 5; % constant \eps_1 in paper, which appears in front of source term for w equation
 eps2 = 0.5; % constant in paper which appears in elliptic equation, nonhomogeneous part 

 % create a struct with input parameters and initial data   
 inputdata = struct('eps1',eps1,'eps2',eps2,'disp',disp,'oneconst', ...
        oneconst,'dk',dk,'T',T,'dx',dx,'toll',toll,'plotornot', ...
        plotornot,'d1_0',d1_0,'d2_0',d2_0,'d3_0',d3_0,'w1_0',w1_0, ...
        'w2_0',w2_0,'w3_0',w3_0,'D',D,'gbcfcn',gbcfcn);
expname = 'exp1_positive_epsilons';
 [d,dp,df, w,wp,mphi,phip,phif,Ep,Em,Ef,damping,ed,ew,tt,tf,p,t,e,max_it,av_it]=wme_fe(inputdata,expname);
 
 h_plot =1/40;
[X,Y] = meshgrid(D(1,1):h_plot:D(1,2));
 no_ts = length(tt);
 no_dp = size(phip,2);
 quarter_pt = round(no_dp/8);
 half_pt = round(no_dp/4);
 quarter = round(no_ts/8);
 half = round(no_ts/4);
 quartertime = tt(quarter);
 halftime = tt(half);
 d_quarter = df(:,:,:,quarter_pt);
 d_half = df(:,:,:,half_pt);
 E_quarter = Ef(:,:,:,quarter_pt);
 E_half = Ef(:,:,:,half_pt);

% plot d at intermediate times
figure;quiver(X,Y,d_half(:,:,1),d_half(:,:,2)); 
axis image;
ylabel('$y$','rotation',0,'FontSize',16,'Interpreter','LateX')
xlabel('$x$','FontSize',16,'Interpreter','LateX')
ax = gca;
ax.FontSize = 12; 
title(['director field,T=',num2str(halftime)])
filename1 = [expname,num2str(halftime),'_d.fig'];
savefig(filename1);
close(gcf);


figure;quiver(X,Y,d_quarter(:,:,1),d_quarter(:,:,2)); 
axis image;
ylabel('$y$','rotation',0,'FontSize',16,'Interpreter','LateX')
xlabel('$x$','FontSize',16,'Interpreter','LateX')
ax = gca;
ax.FontSize = 12; 
title(['director field,T=',num2str(quartertime)])
filename2 = [expname,num2str(quartertime),'_d.fig'];
savefig(filename2);
close(gcf);


% plot electric field
figure;quiver(X,Y,E_half(:,:,1), E_half(:,:,2));
axis image;
ylabel('$y$','rotation',0,'FontSize',16,'Interpreter','LateX')
xlabel('$x$','FontSize',16,'Interpreter','LateX')
ax = gca;
ax.FontSize = 12; 
title(['electric field,T=',num2str(halftime)])
filename3 = [expname,num2str(halftime),'_E.fig'];
savefig(filename3);
close(gcf);

figure;quiver(X,Y,E_quarter(:,:,1), E_quarter(:,:,2));
axis image;
ylabel('$y$','rotation',0,'FontSize',16,'Interpreter','LateX')
xlabel('$x$','FontSize',16,'Interpreter','LateX')
ax = gca;
ax.FontSize = 12; 
title(['electric field,T=',num2str(quartertime)])
filename4 = [expname,num2str(quartertime),'_E.fig'];
savefig(filename4);
close(gcf);


