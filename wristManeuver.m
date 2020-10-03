clc;clear all;close all
%% add the location of the code
warning('off','all');
rmpath(genpath(sprintf('%s%s',pwd)));
warning('on','all');
addpath(sprintf('%s%s',pwd,'\20200821-wristManeuver'));
print_on = false;
do_train = false;
network_suffix = '_10_10';
%% read data
diffsig = table2array(readtable('diffData'));
sumsig = table2array(readtable('sumData'));
pos2sig = table2array(readtable('pos2Data'));
pos3sig = table2array(readtable('pos3Data'));
pos4sig = table2array(readtable('pos4Data'));
pos5sig = table2array(readtable('pos5Data'));
jawPossig = table2array(readtable('jawPosData'));
jawEffortsig = table2array(readtable('jawEffortData'));
proGraspMomentArm = 0.017; %m
%% bias correction
diffsig = diffsig - mean(diffsig(1:500,:));
Nsig = diffsig./sumsig;
Nsig = [Nsig,Nsig.^2];
%% data sampling parameters
numPoint = size(diffsig,1);
Fs = 1500;
Ts = 1/Fs;
t = (1:numPoint)*Ts;
atiLables = {'F_{Nx}','F_{Ny}','F_{Nz}','M_{Nx}','M_{Ny}','M_{Nz}'};
ofsLabels = {'V_{N1}','V_{N2}','V_{N3}','V_{N4}','V_{N5}','V_{N6}'};
%% Training Data and Test Data
x = [pos2sig,Nsig];
load(sprintf('network%s.mat',network_suffix),'net','tr','mean_x','stddev_x','max_ati');
x = (x-mean_x)./stddev_x;
x = x';
%% Find calibration results
net = network(net);
y = net(x);
pred = y'.*max_ati;
%%
% [ha, pos] = tight_subplot(4, 2, 0.065, 0.1, 0.05);
% set(gcf,'units','normalized','OuterPosition',[0 0 1 1]);
% index = reshape(1:6,2,3).';
figure;
counter = 1;
set(gcf,'units','normalized','OuterPosition',[0 0 1 1]);
for axisNum=1:8
    switch axisNum
        case 1
            subplot(4,2,1);
            plot(t,pos4sig*180/pi,'r','linewidth',2);
            hold on
            plot(t,pos5sig*180/pi,'k','linewidth',2);
            plot(t,jawPossig*180/pi,'m','linewidth',2);
            legend('q4','q5','q6');
            xlabel('time(s');
            ylabel('deg');
        case 2
            subplot(4,2,2);
            plot(t,jawEffortsig/proGraspMomentArm,'b','linewidth',2);
            grid on;
            xlabel('time(s');
            ylabel('Jaw Effort (N)');
            ax = gca;
            ax.FontSize = 20;
            set(gca,'linewidth',2)
        case 3
            subplot(4,2,3);
            plot(t,pred(:,1));
        case 4
            subplot(4,2,4);
            plot(t,pred(:,4));
        case 5
            subplot(4,2,5);
            plot(t,pred(:,2));
        case 6
            subplot(4,2,6);
            plot(t,pred(:,5));
        case 7
            subplot(4,2,7);
            plot(t,pred(:,3));
        case 8
            subplot(4,2,8);
            plot(t,pred(:,6));
    end
    grid on;
    ax = gca;
    ax.FontSize = 20;
    set(gca,'linewidth',2)
end

%%
% subplot(3,2,1);
% plot(t,pred(:,1));
% subplot(3,2,2);
% plot(t,pred(:,2));
% subplot(3,2,3);
% plot(t,pred(:,3));
% subplot(3,2,4);
% plot(t,pred(:,4));
% subplot(3,2,5);
% plot(t,pred(:,5));
% subplot(3,2,6);
% plot(t,pred(:,6));
%%
function ati = readForceData()
ati = table2array(readtable('atiData.txt'));
ati = ati-mean(ati(1:500,:));
ati(:,1) = -ati(:,1);
ati(:,3) = -ati(:,3);
ati(:,6) = -ati(:,6);
ati(:,4) = -ati(:,4);
end