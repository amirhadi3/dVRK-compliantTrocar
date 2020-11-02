clc;clear all;close all
%% add the location of the code
warning('off','all');
rmpath(genpath(sprintf('%s%s',pwd)));
warning('on','all');
addpath(sprintf('%s%s',pwd,'\20201030-xy'));
fit_active = true;
print_on = false;
plot_suffix = '_linear';
%% read data
ati = readForceData();
ati(:,4:end) = ati(:,4:end)/1000;  % convert moment units from N.mm to N.m
maxati = max(abs(ati));
atiN = ati./maxati;
diffsig = table2array(readtable('diffData'));
sumsig = table2array(readtable('sumData'));
possig = table2array(readtable('posData'));
%% data sampling parameters
numPoint = size(ati,1);
Fs = 1500;
Ts = 1/Fs;
t = (1:numPoint)*Ts;
%% Training Data and Test Data
Nsig = diffsig./sumsig;
Nsig = Nsig - mean(Nsig(1:500,:));
Nsig = [Nsig,Nsig.^2];
%% Initialize the shaft properties
E = 45E9;            % modulus of elasticity (N/m2)
G = 3E9;              % shear modulus of elasticity (N/m2)

ri = 5.99/2/1000;     % inner radius of the shaft (m)
ro = 8.36/2/1000;     % outer radius of the shaft (m)
L = (400)/1000;       % length of the shaft (m)
Lp = 60e-3;           % Lp is the distance from the base of the instrument
                      % at which the sensor is mounted in (m)

Ixx = 1/4*pi*(ro^4-ri^4);
Iyy = 1/4*pi*(ro^4-ri^4);
Izz = 1/2*pi*(ro^4-ri^4);
%% ks is the spring stiffness in the trocar
ks = 10000;             % ks is the trocar stiffness in N/m
%% solve the optimization function with specified bounds
cx = 6*E*Ixx/ks;
m = size(Nsig,2);

Ls_min = 0;
Ls_max = L;

Ls_offset_min = 0.8*0.400;
Ls_offset_max = 1.2*0.400;

if fit_active
lb = [-inf*ones(6*m,1);0.8*L;0.2*cx;0.2*cx;Lp;Ls_offset_min];
ub = [ inf*ones(6*m,1);1.2*L;10*cx;10*cx;Lp;Ls_offset_max];

f = @(x) upNonLinModel_singleAxis(x,possig,Nsig,atiN,m,maxati);
ip = [rand(6*m,1);L;cx;cx;Lp;Ls_offset_max];

options = optimoptions(@lsqnonlin,'MaxIterations',1000,'Algorithm',...
'trust-region-reflective','Display','iter','MaxFunctionEvaluations',inf);
% [optimal_x,resnorm] = lsqnonlin(f,rand(6*m+5,1),lb,ub,options);
% pswarmopt = optimoptions('particleswarm','SwarmSize',200,'HybridFcn',@fmincon,'MaxIterations',10000,'Display','iter','UseParallel',true,'UseVectorized', false);%,'MaxFunctionEvaluations',inf);
% optimal_x = particleswarm(f,6*m+5,lb,ub,pswarmopt);

problem = createOptimProblem('lsqnonlin','x0',ip,'objective',f,...
    'lb',lb,'ub',ub,'options',options);
ms = MultiStart('PlotFcns',@gsplotbestf);
[optimal_x,errormulti] = run(ms,problem,50);

save('linear_fit.mat','optimal_x')
else
    load('linear_fit.mat')
end
id_A = reshape(optimal_x(1:6*m),6,m);
%% display the results
disp("----------------------------------------------------------------")
fprintf("%-20s%-20s%-20s\n",'Variable','True Value','Identified Value');
fprintf("%-20s%-20f%-20f\n",'L',L,optimal_x(6*m+1));
fprintf("%-20s%-20f%-20f\n",'cx',cx,optimal_x(6*m+2));
fprintf("%-20s%-20f%-20f\n",'cx',cx,optimal_x(6*m+3));
fprintf("%-20s%-20f%-20f\n",'Lp',Lp,optimal_x(6*m+4));
fprintf("%-20s%-20f%-20f\n",'Ls_offset',0,optimal_x(6*m+5));
disp("----------------------------------------------------------------")
disp("Identified A - Sensor calibration Matrix {F = A*Nn}")
disp(id_A)
%%
csf = Nsig*id_A';                                   % forces at the cross-section

Ls = optimal_x(6*m+5)-possig;
ft = ati;
cx = optimal_x(6*m+2);
cy = optimal_x(6*m+3);
Lp = optimal_x(6*m+4);
L = optimal_x(6*m+1);

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

frec = zeros(size(ati));
det1 = F11.*F55-F15.*F51;
det2 = F22.*F44-F24.*F42;

frec(:,1) = (1./det1).*(F55.*csf(:,1)-F15.*csf(:,5));
frec(:,5) = (1./det1).*(-F51.*csf(:,1)+F11.*csf(:,5));
frec(:,3) = csf(:,3);
frec(:,2) = (1./det2).*(F22.*csf(:,2)-F24.*csf(:,4));
frec(:,4) = (1./det2).*(-F42.*csf(:,2)+F44.*csf(:,4));
frec(:,6) = csf(:,6);
%%
os = 1*Fs;
atiplot_data = ft(os:end-os,:).*[1,1,1,1000,1000,1000];
Dplot_data = frec(os:end-os,:).*[1,1,1,1000,1000,1000];
tplot = (1:size(atiplot_data,1))*Ts;
eplot_data = atiplot_data-Dplot_data; 
posplot_data = possig(os:end-os,:);
%% NRMSD calculation
NRMSD = std(eplot_data)./(max(atiplot_data)-min(atiplot_data))*100;
Rsq = 1 - sum((eplot_data).^2)./sum((atiplot_data - mean(atiplot_data)).^2);
rms = std(eplot_data);
pv = max(eplot_data)-min(eplot_data);
%% plotData
plot_Data(tplot,posplot_data,atiplot_data,Dplot_data,print_on,Rsq,plot_suffix)