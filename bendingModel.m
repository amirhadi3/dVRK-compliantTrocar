clc;clear all;close all;
%% Initialize the shaft properties
E = 100E9;              % modulus of elasticity (N/m2)
G = 3E9;                % shear modulus of elasticity (N/m2)

ri = 5.99/2/1000;       % inner radius of the shaft (m)
ro = 8.36/2/1000;       % outer radius of the shaft (m)
L = (430)/1000;         % length of the shaft (m)
Lp = 50e-3;             % Lp is the distance from the base of the instrument
% at which the sensor is mounted in (m)
%% ks is the spring stiffness in the trocar
ks = 10000;             % ks is in N/m
%% construct a shaft object
shaftObj = shaft(L,ri,ro,E,G,ks);
%% read the force data applied at the tip of the instrument
ft = readForceData();
%% construct the instrument penetration into the trocar
Ls_min = 0.160;
Ls_max = L - 0.01;
numofCycles = 5;
Ls_period = (size(ft,1))/numofCycles;
Ls = Ls_min+(Ls_max-Ls_min)*(sin((1:size(ft,1))'/2/Ls_period*2*pi)).^2;
%% calculate the sensor normalized output
shaftObj = shaftObj.ft_to_F(ft,Ls,Lp);
shaftObj = shaftObj.F_to_Nn();
%%
c = 6*shaftObj.E*shaftObj.Ixx/shaftObj.ks;

lb = [-inf*ones(36,1);0.8*L;0.2*c;0.8*Lp];
ub = [ inf*ones(36,1);1.2*L;  5*c;1.2*Lp];

f = @(x) upNonLinModel(x,Ls,shaftObj.Nn,shaftObj.tpF);

options = optimoptions(@lsqnonlin,'MaxIterations',1000,'Algorithm','trust-region-reflective','Display','iter');
[optimal_x,resnorm] = lsqnonlin(f,rand(39,1),lb,ub,options);
optimal_x(39)
inv(reshape(optimal_x(1:36),6,6))
%%
function ati = readForceData()
ati = table2array(readtable('atiData.txt'));
ati = ati-mean(ati(1:500,:));
ati(:,1:2) = -ati(:,1:2);
ati(:,3)= -ati(:,3);
ati(:,6) = ati(:,6);
ati = ati(1000:2000,:);
end