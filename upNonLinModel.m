function y = upNonLinModel(x,Ls,Nn,ft)
% x(1:36) = A(6x6)
% x(37) = L;
% x(38) = c;
% x(39) = Lp;
% x(40) = offset in Ls
%
A = reshape(x(1:36),6,6);
yy = Nn*A';

L = x(37);
cx = x(38);
cy = x(39);
Lp = x(40);
offset_LS = x(41);
Ls = offset_LS+Ls;

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
y = yy-F;