function y = upNonLinModel2(x,Ls,Nn,ft,m,maxati)

% x(1:6*m) = A(6*m)
% x(6*m+1) = L;
% x(6*m+2) = cx;
% x(6*m+3) = cy;
% x(6*m+4) = Lp;
% x(6*m+5) = offset in Ls
% DOFs = an vector of 6x1 which indicates which DOFs to include in the
% calibration example (DOFS = [1,2,6])

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
F22 = 1 - (3*L-Ls).*gx;
F24 = 3*gx;
F42 = (3*L-Ls).*(Ls-Lp).*gx-(L-Lp);
F44 = 1-3*(Ls-Lp).*gx;
F51 = (L-Lp)-(3*L-Ls).*(Ls-Lp).*gy;
F55 = 1-3*(Ls-Lp).*gy;

F(:,1) = F11.*ft(:,1)+F15.*ft(:,5);
F(:,2) = F22.*ft(:,2)+F24.*ft(:,4);
F(:,3) = ft(:,3);
F(:,4) = F42.*ft(:,2)+F44.*ft(:,4);
F(:,5) = F51.*ft(:,1)+F55.*ft(:,5);
F(:,6) = ft(:,6);
%{
for counter = 1:size(ft,1)
    ls = Ls(counter);
    g = ls^2/(c+2*ls^3);
    F11 = 1 - (3*L-ls)*g;
    F15 = -3*g;
    F42 = (3*L-ls)*(ls-Lp)*g-(L-Lp);
    F44 = 1-3*(ls-Lp)*g;
    omega = [0 1;-1 0];
%{
       B is a 6 by 6 matrix mapping forces at the tip to
       the forces at a cross section located at Lp
                
       F(6x1) = B(6x6) * ft(6x1)
       F(nx6) = ft(nx6) * B(6x6)'
%}
    B = [F11*eye(2) zeros(2,1) F15*omega zeros(2,1);
        zeros(1,2) 1 zeros(1,3);
        F42*omega zeros(2,1) F44*eye(2) zeros(2,1);
        zeros(1,5) 1];
    
    F(counter,:) = ft(counter,:)*B.';
end
%}
y = csf-F;