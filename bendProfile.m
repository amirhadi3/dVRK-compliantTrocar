clc;clear all;close all;
%%
figure; set(gcf,'position',[500   390   700   300])
ks = 10000;
Fxt = 0;
Myt = 0.2;
for i=1:numel(ks)
    [z,deltax] = bendProf(ks(i),Myt,Fxt);
%     subplot(2,1,1);
    plot(z,deltax,'linewidth',2);
    hold on;
end
line([0.2,0.35],[0.00015,max(0.00015,deltax(find(z>0.000155,1)))],'LineStyle','--','color','red','linewidth',2)
line([0.2,0.35],[-0.00015,-2*0.00015+max(0.00015,deltax(find(z>0.35,1)))],'LineStyle','--','color','red','linewidth',2)
grid on;
xlim([0,0.45])
v = max(abs(deltax));
ylim([-1.5*v,1.5*v]);
set(gca,'linewidth',2)
xtickangle(90)
ytickangle(90)
set(gca,'FontSize',15);
print(sprintf('bend_%d_%d',Fxt,10*Myt),'-djpeg','-r300')
% Fxt = 0:2:10;
% for i=1:numel(Fxt)
%     subplot(2,1,2)
%     [z,deltax] = bendProf(ks,0,Fxt(i));
%     plot(z,deltax);
%     hold on;
% end
% grid on;
% ylim([-.001 .001])


function [z,deltax] = bendProf(ks,Myt,Fxt)
%% Shaft Properties
E = 45E9;               % modulus of elasticity (N/m2)
G = 3E9;                % shear modulus of elasticity (N/m2)

ri = 5.99/2/1000;       % inner radius of the shaft (m)
ro = 8.36/2/1000;       % outer radius of the shaft (m)
L = (400)/1000;         % length of the shaft (m)
Lp = 50e-3;             % Lp is the distance from the base of the instrument
                        % at which the sensor is mounted in (m)

Iyy = 1/4*pi*(ro^4-ri^4);
%% ks is the spring stiffness in the trocar
cy = 6*E*Iyy/ks;
Ls = 0.35;
z = (0:0.01:L)';
deltas = Ls.^2./(6*E*Iyy+2*ks*Ls.^3)*(3*Myt+(3*L-Ls).*Fxt);
Fs = ks*deltas;
deltax = Myt.*z.^2/2/E/Iyy + Fxt.*z.^2.*(3*L-z)/6/E/Iyy - ...
    Fs.*z.^2/6/E/Iyy.*(3*Ls-z);
end