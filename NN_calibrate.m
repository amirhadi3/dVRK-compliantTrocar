clc;clear all;close all
%% add the location of the code
warning('off','all');
rmpath(genpath(sprintf('%s%s',pwd)));
warning('on','all');
addpath(sprintf('%s%s',pwd,'/20200818-5'));
print_on = false;
do_train = true;
plot_suffix = '_new';
%% read data
ati = readForceData();
diffsig = table2array(readtable('diffData'));
sumsig = table2array(readtable('sumData'));
possig = table2array(readtable('posData'));
%% bias correction
diffsig = diffsig - mean(diffsig(1:500,:));
Nsig = diffsig./sumsig;
Nsig = [Nsig,Nsig.^2];
possig = [possig(1);possig];
x = [possig,Nsig];
%% data sampling parameters
numPoint = size(ati,1);
Fs = 1500;
Ts = 1/Fs;
t = (1:numPoint)*Ts;
atiLables = {'F_{Nx}','F_{Ny}','F_{Nz}','M_{Nx}','M_{Ny}','M_{Nz}'};
ofsLabels = {'V_{N1}','V_{N2}','V_{N3}','V_{N4}','V_{N5}','V_{N6}'};
%%
if do_train
    %Training Data and Test Data
    mean_x = mean(x);
    stddev_x = std(x);
    x = (x-mean_x)./stddev_x;
    x = x';
    max_ati = max(ati);
    ati_train = (ati./max_ati)';
    
    % Create a Fitting Network
    hiddenLayerSize = [10 10];
    net = fitnet(hiddenLayerSize);
    
    % Set up Division of Data for Training, Validation, Testing
    net.divideParam.trainRatio = 70/100;
    net.divideParam.valRatio = 15/100;
    net.divideParam.testRatio = 15/100;
    
    % Train the Network
    [net,tr] = train(net,x,ati_train,'useGPU','yes');%'useParallel','yes','showResources','yes');
    
    % Test the Network
    y = net(x);
    pred = y'.*max_ati;
    save(sprintf('network%s.mat',plot_suffix),'net','tr','mean_x','stddev_x','max_ati')
    % View the Network
    view(net)
else
    load('network_new.mat','net','tr','mean_x','stddev_x','max_ati')
    x = (x-mean_x)./stddev_x;
    x = x';
    % Test the Network
    y = net(x);
    tic
    pred = y'.*max_ati;
    toc
end
%% plot Results
os = 1*Fs;
atiplot_data = ati(os:end-os,:);
Dplot_data = pred(os:end-os,:);
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
%%
function ati = readForceData()
ati = table2array(readtable('atiData.txt'));
ati = ati-mean(ati(1:500,:));
ati(:,1) = -ati(:,1);
ati(:,3) = -ati(:,3);
ati(:,6) = -ati(:,6);
ati(:,4) = -ati(:,4);
end