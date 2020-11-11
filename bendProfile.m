clc;clear all;close all;
%%
ks = 10000:10000:50000;
Fxt = 0;
Myt = 0.5;
for i=1:numel(ks)
    [z,deltax] = bendProf(ks(i),Myt,Fxt);
    subplot(2,1,1);
    plot(z,deltax);
    hold on;
end
grid on;
ylim([-.001 .001])

Fxt = 0:2:10;
for i=1:numel(Fxt)
    subplot(2,1,2)
    [z,deltax] = bendProf(ks,0,Fxt(i));
    plot(z,deltax);
    hold on;
end
grid on;
ylim([-.001 .001])


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