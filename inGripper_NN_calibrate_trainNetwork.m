clc;clear all;close all
%% random number generator control for consistency
% rng.seed = 2020;
trainpath = '/ingripper calibration/fast - 5 cycles - pi_3';
validpath = '/ingripper calibration/fast - single cycle - pi_3';

[xTrain,atiTrain,posTrain] = readInputs(trainpath);
[xValid,atiValid,posValid] = readInputs(validpath);
xTrain = [xTrain;xValid];
atiTrain = [atiTrain;atiValid];

print_on = false;
do_train = true;
plot_suffix = '_30_10';
%% data sampling parameters
numPoint = size(atiTrain,1);
Fs = 1000;
Ts = 1/Fs;
t = (1:numPoint)*Ts;
atiLables = {'F_{Nx}','F_{Ny}','F_{Nz}','M_{Nx}','M_{Ny}','M_{Nz}'};
ofsLabels = {'V_{N1}','V_{N2}','V_{N3}','V_{N4}','V_{N5}','V_{N6}'};
%%
if do_train
    %Training Data and Test Data
    xt = xTrain(1:100:end,:);
    max_ati = max(abs(atiTrain));
    max_ati(3) = 1;
    yt = (atiTrain./max_ati);
    yt = yt(1:100:end,:);
    
    
    yValid = atiValid./max_ati;
    xValids = xValid(1:50:end,:);
    yValids = yValid(1:50:end,:);
    
    % Create a Fitting Network
    hiddenLayerSize = [30 10];
    %%net = fitnet(hiddenLayerSize);
    layers = [...
        featureInputLayer(33)
        fullyConnectedLayer(30)
        sigmoidLayer
        fullyConnectedLayer(10)
        sigmoidLayer
        fullyConnectedLayer(6)
        regressionLayer];
    options = trainingOptions('adam', ...
        'LearnRateSchedule','piecewise', ...
        'LearnRateDropFactor',0.2, ...
        'LearnRateDropPeriod',50, ...
        'MaxEpochs',1000,...
        'InitialLearnRate',1e-3, ...
        'Verbose',false, ...
        'Plots','training-progress', ...
        'ExecutionEnvironment','gpu', ...
        'ValidationData',{xValids,yValids});
    
    % Train the Network
    [net, tr] = trainNetwork(xt,yt,layers,options);
    % net = feedforwardnet(hiddenLayerSize, 'trainlm');
    
    % Test the Network
    predTrain = predict(net,xTrain).*max_ati;
    predValid = predict(net,xValid).*max_ati;
    
    save(sprintf('inGripperNetwork%s.mat',plot_suffix),'net','tr','max_ati');
else
    load(sprintf('inGripperNetwork%s.mat',plot_suffix),'net','tr','max_ati')
    predTrain = predict(net,xTrain).*max_ati;
    predValid = predict(net,xValid).*max_ati;
end
%% plot Results
figure;plotLocal(atiTrain,predTrain)
figure;plotLocal(atiValid,predValid)
% os = 1*Fs;
% atiplot_data = ati(os:end-os,:);
% Dplot_data = pred(os:end-os,:);
% tplot = (1:size(atiplot_data,1))*Ts;
% eplot_data = atiplot_data-Dplot_data;
% posplot_data = pos(os:end-os,:);
% %% NRMSD calculation
% NRMSD = std(eplot_data)./(max(atiplot_data)-min(atiplot_data))*100;
% Rsq = 1 - sum((eplot_data).^2)./sum((atiplot_data - mean(atiplot_data)).^2);
% rms = std(eplot_data);
% pv = max(eplot_data)-min(eplot_data);
% %% plotData
% plot_Data(tplot,posplot_data,atiplot_data,Dplot_data,print_on,Rsq,plot_suffix)
% 
% %% valid:
% atiplot_datavalid = yvalid(os:end-os,:);
% Dplot_datavalid = pred(os:end-os,:);
% tplotvalid = (1:size(atiplot_data,1))*Ts;
% eplot_datavalid = atiplot_data-Dplot_data;
% posplot_datavalid = pos(os:end-os,:);

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
% temp = table2array(readtable('tempData'));
%% bias correction
diffsig = diffsig - mean(diffsig(1:500,:));
ati = ati-mean(ati(1:500,:));
Nsig = diffsig./sumsig;
Nsig = [Nsig,Nsig.^2];
x = [Nsig,pos,vel,eff];
% ati(:,3) = 0;
y = ati;
end

function plotLocal(x,y)
index = [1,3,5,2,4,6];
for i=1:6
    subplot(3,2,index(i))
    plot(x(:,i));hold on; plot(y(:,i))
end
end