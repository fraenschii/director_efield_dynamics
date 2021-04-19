function F=movie2d(mm,tt)
% a function to create movies for 2D variables like the electric potential
% Input:        
%               mm -    vector with function values mm=mm(x,t)
%               tt -    vector with time step values, discrete times where
%                       mm was evaluated
%
% Output:       F  -    struct, containing cdata and colormap

% 2021-04-19

figure('Renderer','zbuffer')
imagesc(mm(:,:,1)); colorbar
axis tight manual
set(gca,'NextPlot','replaceChildren');
[~,~,nt]=size(mm);

F(nt)=struct('cdata',[],'colormap',[]);
for n=1:nt
    imagesc(mm(:,:,n)); colorbar
    axis tight manual 
    title(['t=', num2str(tt(n))])
    F(n)=getframe;
end
