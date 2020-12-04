clc;clear all;close all;
%%
xTrain = [];
atiTrain = [];
for i = 0:49
    clearvars -except xTrain atiTrain i
    %% add the location of the code
    loc = sprintf('/ingripper calibration/Data_1/%d',i)
    warning('off','all');
    rmpath(genpath(sprintf('%s%s',pwd)));
    warning('on','all');
    addpath(sprintf('%s%s',pwd,loc));
    
    %% read data
    ati = table2array(readtable('atiData.txt'));
    diffsig = table2array(readtable('diffData.txt'));
    sumsig = table2array(readtable('sumData.txt'));
    pos = table2array(readtable('posData.txt'));
    eff = table2array(readtable('effData.txt'));
    temp = table2array(readtable('tempData'));
    %% bias correction
    diffsig = diffsig - mean(diffsig(1:500,:));
    ati = ati-mean(ati(1:500,:));
    Nsig = diffsig./sumsig;
    Nsig = [Nsig,Nsig.^2];
    x = [pos,eff];
    xTrain = [xTrain;x];
    atiTrain = [atiTrain;ati];
end
save('data_pe.mat','xTrain','atiTrain')

