clc;clear all;close all
%% random number generator control for consistency
rng.seed = 2020;
load 'data_pe.mat'
atiTrain = [atiTrain(:,1:2),atiTrain(:,4:6)];
print_on = true;
do_train = false;
plot_suffix = '_10_5_s_pe';
%% data sampling parameters
numPoint = size(atiTrain,1);
Fs = 1000;
Ts = 1/Fs;
t = (1:numPoint)*Ts;
%% network training
if do_train
    hiddenLayerSize = [10 5];
end
if do_train
    %Training Data and Test Data
    mean_x = mean(xTrain);
    stddev_x = std(xTrain);
    xmn = meanNormalize(xTrain,mean_x,stddev_x);
    
    ts = 1 ; te = floor(0.2*length(xmn));numelt = te-ts+1;
    vs = floor(0.2*length(xmn))+1; ve = floor(0.3*length(xmn));numelv = ve-vs+1;
    tes = floor(0.9*length(xmn))+1; tee = floor(1*length(xmn));numelte = tee-tes+1;
    
    m=50;
    xt = [xmn(ts:m:te,:);xmn(tes:m:tee,:);xmn(vs:m:ve,:)];
    xTest = xmn(vs:ve,:);

    
    max_ati = max(abs(atiTrain));
    
    ys = (atiTrain./max_ati);
    yt = [ys(ts:m:te,:);ys(tes:m:tee,:);ys(vs:m:ve,:)];
    yTest = ys(vs:ve,:);
    
    % Create a Fitting Network
    net = fitnet(hiddenLayerSize,'trainlm');
   
    % Set up Division of Data for Training, Validation, Testing
    net.divideFcn = 'divideblock';
    
    net.divideParam.trainRatio = numelt/(numelt+numelte+numelv);
    net.divideParam.valRatio = numelv/(numelt+numelte+numelv);
    net.divideParam.testRatio = numelte/(numelt+numelte+numelv);
    net.trainParam.epochs = 1000;  
    
    % Train the Network
    tic
    [net,tr] = train(net,xt',yt');%,'useGPU','yes','useParallel','yes','showResources','yes');
    toc
    % Test the Network
    %     xmnValid = meanNormalize(xValid,mean_x,stddev_x);
    predTrain = net(xmn')'.*max_ati;
    predValid = net(xTest')'.*max_ati;
    
    save(sprintf('inGripperNetwork%s.mat',plot_suffix),'net','tr','max_ati','mean_x','stddev_x');
else
    load(sprintf('inGripperNetwork%s.mat',plot_suffix),'net','tr','max_ati','mean_x','stddev_x')
    xmn = meanNormalize(xTrain,mean_x,stddev_x);
    tes = floor(0.9*length(xmn))+1; tee = floor(1*length(xmn));numelte = tee-tes+1;
    xTest = xmn(tes:tee,:);
    predTest = net(xTest')'.*max_ati;
    yTest = atiTrain(tes:tee,:);
end
%% plot Results
atiplot_data = yTest;
Dplot_data = predTest;
tplot = (1:size(atiplot_data,1))*Ts;
eplot_data = atiplot_data-Dplot_data;
posplot_data = [];
%% NRMSD calculation
NRMSD = std(eplot_data)./(max(atiplot_data)-min(atiplot_data))*100;
Rsq = 1 - sum((eplot_data).^2)./sum((atiplot_data - mean(atiplot_data)).^2);
rms = std(eplot_data);
pv = max(eplot_data)-min(eplot_data);
%% plotData
plot_Data5(tplot,posplot_data,atiplot_data,Dplot_data,print_on,Rsq,plot_suffix)
%% functions
function [x,y,pos] = readInputs(loc)
%% add the location of the code
warning('off','all');
rmpath(genpath(sprintf('%s%s',pwd)));
warning('on','all');
addpath(sprintf('%s%s',pwd,loc));

%% read data
% if exist('x','var') ~= 1
ati = table2array(readtable('atiData'));
diffsig = table2array(readtable('diffData'));
sumsig = table2array(readtable('sumData'));
pos = table2array(readtable('posData'));
vel = table2array(readtable('velData'));
eff = table2array(readtable('effData'));
temp = table2array(readtable('tempData'));
%% bias correction
diffsig = diffsig - mean(diffsig(1:500,:));
ati = ati-mean(ati(1:500,:));
Nsig = diffsig./sumsig;
Nsig = [Nsig,Nsig.^2];
x = [Nsig,pos,vel,eff,temp];
y = ati;
end

function plotLocal(x,y)
index = [1,3,5,2,4,6];
for i=1:6
    subplot(3,2,index(i))
    plot(x(:,i));hold on; plot(y(:,i))
end
end

function xm = meanNormalize(x,mean_x,stddev_x)
xm = (x-mean_x)./stddev_x;
end