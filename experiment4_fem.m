% Experiment 4
% Initial data numerical experiments director field coupled to electric
% field
% 2025-06-14
% Experiment 1
% blow up example, damping dk=3, 

 disp = 0.5; % constant in front of d_{tt}
 dk = 3; % constant in front of damping term d_t
 oneconst = 1; % constant k, in front of elastic term |\Grad d|^2 in paper
 T = 1; % final time of simulation
 
 dx = 1/64; % grid size
 toll = 0.05*dx^2; % tolerance in nonlinear iteration (stopping tolerance)
 plotornot = 0; % whether to plot the director field while simulating or not


  % initial data for director field
     af=@(r)((1-2*r).^4);
     rf=@(x,y)(sqrt(x.^2+y.^2));
     denominatorf=@(x,y)((rf(x,y).^2+af(rf(x,y)).^2));
     d1_0=@(x,y)(2*(rf(x,y)<=0.5).*(x).*af(rf(x,y))./(denominatorf(x,y)));
     d2_0=@(x,y)(2*(rf(x,y)<=0.5).*(y).*af(rf(x,y))./(denominatorf(x,y)));
     d3_0=@(x,y)((rf(x,y)<=0.5).*(af(rf(x,y)).^2-rf(x,y).^2)./(denominatorf(x,y))-(rf(x,y)>0.5));
  
  % initial data for angular momentum w   
  w1_0=@(x,y)(zeros(size(x)));
  w2_0=@(x,y)(zeros(size(x)));
  w3_0=@(x,y)(zeros(size(x)));

 gbcfcn = @(t,x,y)(3*sin(2*pi*t+0.2)*(x+0.5).*sin(pi*y)); % boundary condition for elliptic equation (potential)
    
 D=[-0.5,0.5;-0.5,0.5]; % computational domain
 
 eps1 = -5; % constant \eps_1 in paper, which appears in front of source term for w equation
 eps2 = -0.2; % constant in paper which appears in elliptic equation, nonhomogeneous part 

 % create a struct with input parameters and initial data   
 inputdata = struct('eps1',eps1,'eps2',eps2,'disp',disp,'oneconst', ...
        oneconst,'dk',dk,'T',T,'dx',dx,'toll',toll,'plotornot', ...
        plotornot,'d1_0',d1_0,'d2_0',d2_0,'d3_0',d3_0,'w1_0',w1_0, ...
        'w2_0',w2_0,'w3_0',w3_0,'D',D,'gbcfcn',gbcfcn);
expname = 'exp4_blowup_largedamping';
 [d,dp,df, w,wp,mphi,phip,phif,Ep,Em,Ef,damping,ed,ew,tt,tf,p,t,e,max_it,av_it]=wme_fe(inputdata,expname);
 h_plot =1/40;
[X,Y] = meshgrid(D(1,1):h_plot:D(1,2));
 no_ts = length(tt);
 no_dp = size(phip,2);
 quarter_pt = round(no_dp/4);
 threequarter_pt = round(3*no_dp/4);
 half_pt = round(no_dp/2);
 quarter = round(no_ts/4);
 half = round(no_ts/2);
 threequarter = round(3*no_ts/4);
 quartertime = tt(quarter);
 halftime = tt(half);
 threequartertime = tt(threequarter);
 d_quarter = df(:,:,:,quarter_pt);
 d_3quarter = df(:,:,:,threequarter_pt);
 d_half = df(:,:,:,half_pt);
 E_quarter = Ef(:,:,:,quarter_pt);
 E_3quarter = Ef(:,:,:,threequarter_pt);
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

figure;quiver(X,Y,d_quarter(:,:,1),d_quarter(:,:,2)); 
axis image;
ylabel('$y$','rotation',0,'FontSize',16,'Interpreter','LateX')
xlabel('$x$','FontSize',16,'Interpreter','LateX')
ax = gca;
ax.FontSize = 12; 
title(['director field,T=',num2str(quartertime)])
filename2 = [expname,num2str(quartertime),'_d.fig'];
savefig(filename2);

figure;quiver(X,Y,d_3quarter(:,:,1),d_3quarter(:,:,2)); 
axis image;
ylabel('$y$','rotation',0,'FontSize',16,'Interpreter','LateX')
xlabel('$x$','FontSize',16,'Interpreter','LateX')
ax = gca;
ax.FontSize = 12; 
title(['director field,T=',num2str(threequartertime)])
filename5 = [expname,num2str(threequartertime),'_d.fig'];
savefig(filename5);

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

figure;quiver(X,Y,E_quarter(:,:,1), E_quarter(:,:,2));
axis image;
ylabel('$y$','rotation',0,'FontSize',16,'Interpreter','LateX')
xlabel('$x$','FontSize',16,'Interpreter','LateX')
ax = gca;
ax.FontSize = 12; 
title(['electric field,T=',num2str(quartertime)])
filename4 = [expname,num2str(quartertime),'_E.fig'];
savefig(filename4);


figure;quiver(X,Y,E_3quarter(:,:,1), E_3quarter(:,:,2));
axis image;
ylabel('$y$','rotation',0,'FontSize',16,'Interpreter','LateX')
xlabel('$x$','FontSize',16,'Interpreter','LateX')
ax = gca;
ax.FontSize = 12; 
title(['electric field,T=',num2str(threequartertime)])
filename6 = [expname,num2str(threequartertime),'_E.fig'];
savefig(filename6);


