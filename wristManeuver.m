clc;clear all;close all
%% add the location of the code
warning('off','all');
rmpath(genpath(sprintf('%s%s',pwd)));
warning('on','all');
addpath(sprintf('%s%s',pwd,'\20200821-wristManeuver'));
print_on = true;
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
close all;
figure(1);
set(gcf,'units','normalized','position',[0.4 0.2 0.3 0.6]);
for axisNum=1:8
    switch axisNum
        case 1
            subplot(4,2,1);
            plot(t,pos4sig*180/pi,'b','linewidth',2);
            hold on
            plot(t,pos5sig*180/pi,'r','linewidth',2);
            plot(t,jawPossig*180/pi,'k','linewidth',2);
            legend('q4','q5','q6','NumColumns',3);
            ylabel('deg.','fontweight','b');
        case 2
            subplot(4,2,3);
            plot(t,pred(:,1),'b','linewidth',2);
            ylabel('F_x (N)','fontweight','b');
        case 3
            subplot(4,2,5);
            plot(t,pred(:,2),'b','linewidth',2);
            ylabel('F_y (N)','fontweight','b');
        case 4
            subplot(4,2,7);
            plot(t,pred(:,3),'b','linewidth',2);
            ylabel('F_z (N)','fontweight','b');
        case 5
            subplot(4,2,2);
            plot(t,jawEffortsig/proGraspMomentArm,'b','linewidth',2);
            ylabel('F_{jaw} (N)','fontweight','b');
        case 6
            subplot(4,2,4);
            plot(t,pred(:,4)/1000,'b','linewidth',2);
            ylabel('M_x (N.m)','fontweight','b');
        case 7
            subplot(4,2,6);
            plot(t,pred(:,5)/1000,'b','linewidth',2);
            ylabel('M_y (N.m)','fontweight','b');
        case 8
            subplot(4,2,8);
            plot(t,pred(:,6)/1000,'b','linewidth',2);
            ylabel('M_z (N.m)','fontweight','b');
    end
    
    grid on;
    ax = gca;
    ax.FontSize = 18;
    set(gca,'linewidth',2)
    
    if axisNum == 4 || axisNum == 8
        set(gca,'XTickLabel',xticklabel);
        xlabel('Time(s)','fontweight','b');
    else
        xticklabel = get(gca,'XTickLabel');
        set(gca,'XTickLabel',[]);
    end
    
end

ha=get(gcf,'children');
set(ha(1),'position',[0.6014 0.1100 0.35 0.17])
set(ha(2),'position',[0.6014 0.3300 0.35 0.17])
set(ha(3),'position',[0.6014 0.5500 0.35 0.17])
set(ha(4),'position',[0.6014 0.7700 0.35 0.17])
set(ha(5),'position',[0.1300 0.1100 0.35 0.17])
set(ha(6),'position',[0.1300 0.3300 0.35 0.17])
set(ha(7),'position',[0.1300 0.5500 0.35 0.17])
set(ha(8),'position',[0.159 0.945 0.2969 0.0359]);
set(ha(9),'position',[0.1300 0.7700 0.35 0.17])
set(gcf,'color','w')
set(gcf, 'InvertHardCopy', 'off');
if print_on
    print('wristManuever','-djpeg','-r600');
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