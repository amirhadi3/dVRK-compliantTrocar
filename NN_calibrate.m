clc;clear all;close all
%% add the location of the code
warning('off','all');
rmpath(genpath(sprintf('%s%s',pwd)));
warning('on','all');
%%
print_on = true;
do_train = true;
plot_suffix = '_wGrip_5';
%% read lateral0, axial 0, lateral1, axial 1
loc = {'data_shaft/10_cycle_wGrip1','data_shaft/10_cycle_wGrip0'};
% loc = {'data_shaft/10_cycle_NoGrip'};
ati = [];
diffsig = [];
sumsig = [];
pos = [];
index = [];
eff = [];
for i = 1:length(loc)
    ati = [ati;table2array(readtable(sprintf('%s%s',loc{i},'/atiData')))];
    diffsig = [diffsig;table2array(readtable(sprintf('%s%s',loc{i},'/diffData')))];
    sumsig = [sumsig;table2array(readtable(sprintf('%s%s',loc{i},'/sumData')))];
    pos = [pos;table2array(readtable(sprintf('%s%s',loc{i},'/posData')))];
    eff = [eff;table2array(readtable(sprintf('%s%s',loc{i},'/effData')))];
    index = [index;length(ati)];
end
%% data sampling parameters
numPoint = size(ati,1);
Fs = 1000;
Ts = 1/Fs;
t = (1:numPoint)*Ts;
atiLables = {'F_{Nx}','F_{Ny}','F_{Nz}','M_{Nx}','M_{Ny}','M_{Nz}'};
ofsLabels = {'V_{N1}','V_{N2}','V_{N3}','V_{N4}','V_{N5}','V_{N6}'};
%%
ati = ati-mean(ati(1:500,:));
% ati(:,4:end) = ati(:,4:end)/1000;               % convert moment units from N.mm to N.m
ati = ati(:,[1:2,4:6]);
diffsig = diffsig-mean(diffsig(1:500,:));
possig = pos(:,[3,4]);
effsig = eff(:,end);
Nsig = diffsig./sumsig;
Nsig = [Nsig,Nsig.^2];
q3vel = [[0,0];diff(pos(:,[3,4]))]/Ts;
x = [possig,effsig,Nsig];
% x = [possig,Nsig];
%%
%%
if do_train
    %Training Data and Test Data
    mean_x = mean(x);
    stddev_x = std(x);
    x = (x-mean_x)./stddev_x;
    xt = x(1:25:end,:)';
    x = x';
    max_ati = max(abs(ati));
    ati_train = (ati(1:25:end,:)./max_ati)';
    
    
    % Create a Fitting Network
    hiddenLayerSize = [5];
    net = fitnet(hiddenLayerSize);
    
    net.divideFcn = 'divideblock';
    
    net.divideParam.trainRatio = 0.6;
    net.divideParam.valRatio = 0.2;
    net.divideParam.testRatio = 0.2;
    net.trainParam.epochs = 1000;
    
    % Train the Network
    tic 
    [net,tr] = train(net,xt,ati_train);%,'useGPU','yes','useParallel','yes','showResources','yes');
    toc
    
    % Test the Network
    y = net(x);
    pred = y'.*max_ati;
    save(sprintf('network%s.mat',plot_suffix),'net','tr','mean_x','stddev_x','max_ati')
    
    % View the Network
    %view(net)
else
    load(sprintf('network%s.mat',plot_suffix),'net','tr','mean_x','stddev_x','max_ati')
    net = network(net);
    x = (x-mean_x)./stddev_x;
    x = x';
    % Test the Network
    y = net(x);
    tic
    pred = y'.*max_ati;
    toc
end
%% plot Results
os = floor(0.8*numPoint);
atiplot_data = ati(os:end,:);
Dplot_data = pred(os:end,:);
tplot = (1:size(atiplot_data,1))*Ts;
eplot_data = atiplot_data-Dplot_data;
posplot_data = possig(os:end,:);
%% NRMSD calculation
NRMSD = std(eplot_data)./(max(atiplot_data)-min(atiplot_data))*100;
Rsq = 1 - sum((eplot_data).^2)./sum((atiplot_data - mean(atiplot_data)).^2);
rms = std(eplot_data);
pv = max(eplot_data)-min(eplot_data);
%% plotData
plot_Data5(tplot,posplot_data,atiplot_data,Dplot_data,print_on,Rsq,plot_suffix)