clc;clear all;close all
%% add the location of the code
warning('off','all');
rmpath(genpath(sprintf('%s%s',pwd)));
warning('on','all');
%% read lateral0, axial 0, lateral1, axial 1
loc = {'data_shaft/Lateral0','data_shaft/Axial0','data_shaft/Lateral1','data_shaft/Axial1'};
ati = [];
diffsig = [];
sumsig = [];
pos = [];
index = [];
for i = 1:length(loc)
    ati = [ati;table2array(readtable(sprintf('%s%s',loc{i},'/atiData')))];
    diffsig = [diffsig;table2array(readtable(sprintf('%s%s',loc{i},'/diffData')))];
    sumsig = [sumsig;table2array(readtable(sprintf('%s%s',loc{i},'/sumData')))];
    pos = [pos;table2array(readtable(sprintf('%s%s',loc{i},'/posData')))];
    index = [index;length(ati)];
end
%% preprocessing
ati = ati-mean(ati(1:500,:));
ati(:,4:end) = ati(:,4:end)/1000;               % convert moment units from N.mm to N.m
diffsig = diffsig-mean(diffsig(1:500,:));
possig = pos(:,3);
maxati = max(abs(ati));
atiN = ati./maxati;
Nsig = diffsig./sumsig;
Nsig = [Nsig,Nsig.^2];
%% Fit and save settings
fit_active = true;
print_on = false;
plot_suffix = '_linear6';
%% data sampling parameters
numPoint = size(ati,1);
Fs = 1000;
Ts = 1/Fs;
t = (1:numPoint)*Ts;
%% Initialize the shaft properties
E = 45E9;            % modulus of elasticity (N/m2)
G = 3E9;              % shear modulus of elasticity (N/m2)

ri = 5.99/2/1000;     % inner radius of the shaft (m)
ro = 8.36/2/1000;     % outer radius of the shaft (m)
L = (400)/1000;       % length of the shaft (m)
Lp = 50e-3;           % Lp is the distance from the base of the instrument
% at which the sensor is mounted in (m)

Ixx = 1/4*pi*(ro^4-ri^4);
Iyy = 1/4*pi*(ro^4-ri^4);
Izz = 1/2*pi*(ro^4-ri^4);
%% ks is the spring stiffness in the trocar
ks = 10000;             % ks is the trocar stiffness in N/m
%% solve the optimization function with specified bounds
cx = 6*E*Ixx/ks;
m = size(Nsig,2);

Ls_offset_min =0;
Ls_offset_max =0.400;

if fit_active
    lb = [-inf*ones(6*m,1);0.5*L;0.2*cx;0.2*cx;Lp;Ls_offset_min];
    ub = [ inf*ones(6*m,1);1.5*L;10*cx;10*cx;Lp;Ls_offset_max];
    dofselector = [1,1,1,1,1,1];
    f = @(x) upNonLinModel(instrument(x,m),possig(1:50:end),Nsig(1:50:end,:),atiN(1:50:end,:),maxati,dofselector);
    ip = [rand(6*m+5,1)];%,1);L;cx;cx;Lp;2/3*Ls_offset_max];
    options = optimoptions(@lsqnonlin,'MaxIterations',1000,'Algorithm',...
        'trust-region-reflective','Display','iter');
    problem = createOptimProblem('lsqnonlin','x0',ip,'objective',f,...
        'lb',lb,'ub',ub,'options',options);
    ms = MultiStart('PlotFcns',@gsplotbestf);
    [optimal_x,errormulti] = run(ms,problem,50);
    save('linear_fit.mat','optimal_x')
else
    load('linear_fit.mat')
end

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
% disp(id_A)
%%
instObj = instrument(optimal_x,m);
instObj = instObj.Nn_to_wt(Nsig,possig);
frec = instObj.wt;
%%
os = 1*Fs;
atiplot_data = ati(os:end-os,:).*[1,1,1,1000,1000,1000];
Dplot_data = frec(os:end-os,:).*[1,1,1,1000,1000,1000];
tplot = (1:size(atiplot_data,1))*Ts;
eplot_data = atiplot_data-Dplot_data;
posplot_data = pos(os:end-os,[3,4]);
%% NRMSD calculation
NRMSD = std(eplot_data)./(max(atiplot_data)-min(atiplot_data))*100;
Rsq = 1 - sum((eplot_data).^2)./sum((atiplot_data - mean(atiplot_data)).^2);
rms = std(eplot_data);
pv = max(eplot_data)-min(eplot_data);
%% plotData
plot_Data(tplot,posplot_data,atiplot_data,Dplot_data,print_on,Rsq,plot_suffix)