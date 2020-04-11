function y = upNonLinModel(x,Ls,Nn,ft)
% x(1:36) = A(6x6)
% x(37) = L;
% x(38) = c;
% x(39) = Lp;
%
A = reshape(x(1:36),6,6);
yy = Nn*A';

L = x(37);
c = x(38);
Lp = x(39);

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

y = yy-F;