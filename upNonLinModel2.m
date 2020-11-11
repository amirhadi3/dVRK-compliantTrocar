function y = upNonLinModel2(x,Ls,Nn,ft,m,maxati)

% x(1:6*m) = A(6*m)
% x(6*m+1) = L;
% x(6*m+2) = cx;
% x(6*m+3) = cy;
% x(6*m+4) = Lp;
% x(6*m+5) = offset in Ls

A = reshape(x(1:6*m),6,m);
csf = Nn*A'; %%dim(yy) = n*m*m*6 = n*6

L = x(6*m+1);
cx = x(6*m+2);
cy = x(6*m+3);
Lp = x(6*m+4);
offset_LS = x(6*m+5);

Ls = offset_LS-Ls;

gx = Ls.^2./(cx+2*Ls.^3);
gy = Ls.^2./(cy+2*Ls.^3);

F11 = 1-(3*L-Ls).*gy;
F15 = -3*gy;
F22 = 1-(3*L-Ls).*gx;
F24 = 3*gx;
F42 = (3*L-Ls).*(Ls-Lp).*gx-(L-Lp);
F44 = 1-3*(Ls-Lp).*gx;
F51 = (L-Lp)-(3*L-Ls).*(Ls-Lp).*gy;
F55 = 1-3*(Ls-Lp).*gy;

det1 = F11.*F55-F15.*F51;
det2 = F22.*F44-F24.*F42;

F(:,1) = (1./det1).*(F55.*csf(:,1)-F15.*csf(:,5))/maxati(1);
F(:,5) = (1./det1).*(-F51.*csf(:,1)+F11.*csf(:,5))/maxati(5);
F(:,3) = csf(:,3)/maxati(3);
F(:,2) = (1./det2).*(F44.*csf(:,2)-F24.*csf(:,4))/maxati(2);
F(:,4) = (1./det2).*(-F42.*csf(:,2)+F22.*csf(:,4))/maxati(4);
F(:,6) = csf(:,6)/maxati(6);
% parfor i = 1:size(ft,1)
% %{
%        B is a 6 by 6 matrix mapping forces at the tip to
%        the forces at a cross section located at Lp
%                 
%        F(6x1) = B(6x6) * ft(6x1)
%        F(nx6) = ft(nx6) * B(6x6)'
% %}
%     H = [F11(i) 0 0 0 F15(i) 0;
%          0 F22(i) 0 F24(i) 0 0;
%          0 0 1 0 0 0;
%          0 F42(i) 0 F44(i) 0 0;
%          F51(i) 0 0 0 F55(i) 0;
%          0 0 0 0 0 1];
%     
%     F(i,:) = (pinv(H)*csf(i,:)')'*mat;
% end
y = (F-ft)*diag([1,1,3,10,10,1]);