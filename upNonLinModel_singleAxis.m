function y = upNonLinModel_singleAxis(x,Ls,Nn,ft,m,maxati)

A = reshape(x(1:6*m),6,m);
csf = Nn*A';

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
F(:,2) = (1./det2).*(F22.*csf(:,2)-F24.*csf(:,4))/maxati(2);
F(:,4) = (1./det2).*(-F42.*csf(:,2)+F44.*csf(:,4))/maxati(4);
F(:,6) = csf(:,6)/maxati(6);

y = F(:,[1,2])-ft(:,[1,2]);