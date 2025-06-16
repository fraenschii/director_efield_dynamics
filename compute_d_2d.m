function dn=compute_d_2d(d,w,dt) 
% evolution/update for equation d_t = d x w
% Input: d,w previous nodal values and dt time step
% Output: dn new nodal values, dn =d + dt*0.5* (d+dn) x w
% one step evolution of d 
dn=zeros(size(d));
dtm =dt/2;
dtm2 = dtm^2;
c1=d(:,1)+dtm* (d(:,2).*w(:,3)-d(:,3).*w(:,2));
c2=d(:,2)+dtm* (d(:,3).*w(:,1)-d(:,1).*w(:,3));
c3=d(:,3)+dtm* (d(:,1).*w(:,2)-d(:,2).*w(:,1));

a11=1+dtm2*(w(:,1).^2);
a12=dtm2*w(:,2).*w(:,1)+dtm*w(:,3);
a13=dtm2*w(:,1).*w(:,3)-dtm*w(:,2);
a21=dtm2*w(:,1).*w(:,2)-dtm*w(:,3);
a22=1+dtm2*(w(:,2).^2);
a23=dtm2*w(:,2).*w(:,3)+dtm*w(:,1);
a31=dtm2*w(:,1).*w(:,3)+dtm*w(:,2);
a32=dtm2*w(:,2).*w(:,3)-dtm*w(:,1);
a33=1+dtm2*(w(:,3).^2);
detn=1+dtm2*((w(:,1).^2)+(w(:,2).^2)+(w(:,3).^2));
dn(:,1)=(a11.*c1+a12.*c2+a13.*c3)./detn;
dn(:,2)=(a21.*c1+a22.*c2+a23.*c3)./detn;
dn(:,3)=(a31.*c1+a32.*c2+a33.*c3)./detn;

end