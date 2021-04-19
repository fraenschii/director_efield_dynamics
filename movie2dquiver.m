function F=movie2dquiver(xx,mm,tt)
% a function to create movies using quiver plot (for director and electric
% field)
% Input:        xx -    vector with grid points in spatial dimension
%               mm -    vector with function values mm=mm(x,t)
%               tt -    vector with time step values, discrete times where
%                       mm was evaluated

% Output:       F  -    struct, containing cdata and colormap

% 2021-04-19
[X,Y]=meshgrid(xx);
figure('Renderer','zbuffer')
quiver(X,Y,mm(:,:,1,1),mm(:,:,2,1));
   axis tight manual
set(gca,'NextPlot','replaceChildren');
[~,~,~,nt]=size(mm);

F(nt)=struct('cdata',[],'colormap',[]);
for n=1:nt
     quiver(X,Y,mm(:,:,1,n),mm(:,:,2,n));
    axis image;
    title(['t=', num2str(tt(n))])
    F(n)=getframe;
end