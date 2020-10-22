clc;clear all;close all
%% add the location of the code
print_on = true;
%% read data
[t0,jawPossig0,jawEffortsig0,pred0] = estimate('20201021-pinch0');
[t1,jawPossig1,jawEffortsig1,pred1] = estimate('20201021-pinch1');
[t2,jawPossig2,jawEffortsig2,pred2] = estimate('20201021-pinch2');
[t4,jawPossig4,jawEffortsig4,pred4] = estimate('20201021-pinch4');

index0 = find(abs(jawPossig0-jawPossig0(1))>1/180*3.1415,1);
index1 = find(abs(jawPossig1-jawPossig1(1))>1/180*3.1415,1);
index2 = find(abs(jawPossig2-jawPossig2(1))>1/180*3.1415,1);
index4 = find(abs(jawPossig4-jawPossig4(1))>1/180*3.1415,1);

offset = 2000;
t0 = t0(index0-offset:end)-t0(index0);
jawPossig0 = jawPossig0(index0-offset:end);
jawEffortsig0 = jawEffortsig0(index0-offset:end);
pred0 = pred0(index0-offset:end,:);

t1 = t1(index1-offset:end)-t1(index1);
jawPossig1 = jawPossig1(index1-offset:end);
jawEffortsig1 = jawEffortsig1(index1-offset:end);
pred1 = pred1(index1-offset:end,:);

t2 = t2(index2-offset:end)-t2(index2);
jawPossig2 = jawPossig2(index2-offset:end);
jawEffortsig2 = jawEffortsig2(index2-offset:end);
pred2 = pred2(index2-offset:end,:);

t4 = t4(index4-offset:end)-t4(index4);
jawPossig4 = jawPossig4(index4-offset:end);
jawEffortsig4 = jawEffortsig4(index4-offset:end);
pred4 = pred4(index4-offset:end,:);

proGraspMomentArm = 0.017; %m
%%
close all;
figure(1);
set(gcf,'units','normalized','position',[0.4 0.2 0.3 0.6]);
for axisNum=1:8
    switch axisNum
        case 1
            subplot(4,2,1);
            plot(t0,jawPossig0*180/pi,'linewidth',2);hold on;
            plot(t1,jawPossig1*180/pi,'b','linewidth',2);
            plot(t2,jawPossig2*180/pi,'r','linewidth',2);
            plot(t4,jawPossig4*180/pi,'k','linewidth',2);
            ylabel('q_6(deg.)','fontweight','b');
            legend('no foam','1-layer','2-layers','4-layers','NumColumns',4)
        case 2
            subplot(4,2,3);
            plot(t0,pred0(:,1),'linewidth',1);hold on;
            plot(t1,pred1(:,1),'b','linewidth',1);
            plot(t2,pred2(:,1),'r','linewidth',1);
            plot(t4,pred4(:,1),'k','linewidth',2);
            ylabel('F_x (N)','fontweight','b');
        case 3
            subplot(4,2,5);
            plot(t0,pred0(:,2),'linewidth',2);hold on;
            plot(t1,pred1(:,2),'b','linewidth',2);
            plot(t2,pred2(:,2),'r','linewidth',2);
            plot(t4,pred4(:,2),'k','linewidth',2);
            ylabel('F_y (N)','fontweight','b');
        case 4
            subplot(4,2,7);
            plot(t0,pred0(:,3),'linewidth',2);hold on;
            plot(t1,pred1(:,3),'b','linewidth',2);
            plot(t2,pred2(:,3),'r','linewidth',2);
            plot(t4,pred4(:,3),'k','linewidth',2);
            ylabel('F_z (N)','fontweight','b');
        case 5
            subplot(4,2,2);
            plot(t0,jawEffortsig0/proGraspMomentArm,'linewidth',2);hold on;
            plot(t1,jawEffortsig1/proGraspMomentArm,'b','linewidth',2);
            plot(t2,jawEffortsig2/proGraspMomentArm,'r','linewidth',2);
            plot(t4,jawEffortsig4/proGraspMomentArm,'k','linewidth',2);
            ylabel('F_{jaw} (N)','fontweight','b');
        case 6
            subplot(4,2,4);
            plot(t0,pred0(:,4)/1000,'linewidth',2);hold on;
            plot(t1,pred1(:,4)/1000,'b','linewidth',2);
            plot(t2,pred2(:,4)/1000,'r','linewidth',2);
            plot(t4,pred4(:,4)/1000,'k','linewidth',2);
            ylabel('M_x (N.m)','fontweight','b');
        case 7
            subplot(4,2,6);
            plot(t0,pred0(:,5)/1000,'linewidth',2);hold on;
            plot(t1,pred1(:,5)/1000,'b','linewidth',2);
            plot(t2,pred2(:,5)/1000,'r','linewidth',2);
            plot(t4,pred4(:,5)/1000,'k','linewidth',2);
            ylabel('M_y (N.m)','fontweight','b');
        case 8
            subplot(4,2,8);
            plot(t0,pred0(:,6)/1000,'linewidth',2);hold on;
            plot(t1,pred1(:,6)/1000,'b','linewidth',2);
            plot(t2,pred2(:,6)/1000,'r','linewidth',2);
            plot(t4,pred4(:,6)/1000,'k','linewidth',2);
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
set(ha(1),'position',[0.6014 0.1100 0.32 0.17])
set(ha(2),'position',[0.6014 0.3100 0.32 0.17])
set(ha(3),'position',[0.6014 0.5100 0.32 0.17])
set(ha(4),'position',[0.6014 0.7100 0.32 0.17])
set(ha(5),'position',[0.1300 0.1100 0.32 0.17])
set(ha(6),'position',[0.1300 0.3100 0.32 0.17])
set(ha(7),'position',[0.1300 0.5100 0.32 0.17])
set(ha(8),'position',[0.1300 0.9 0.6014+0.32-0.13 0.0359]);
set(ha(9),'position',[0.1300 0.7100 0.32 0.17])
set(gcf,'color','w')
set(gcf, 'InvertHardCopy', 'off');
if print_on
    print('pinch','-djpeg','-r600');
end

%%
function [t,jawPossig,jawEffortsig,pred] = estimate(file)
warning('off','all');
rmpath(genpath(sprintf('%s%s',pwd)));
warning('on','all');
addpath(sprintf('%s\\%s',pwd,file));
diffsig = table2array(readtable('diffData'));
sumsig = table2array(readtable('sumData'));
pos2sig = table2array(readtable('posData'));
pos3sig = table2array(readtable('pos3Data'));
pos4sig = table2array(readtable('pos4Data'));
pos5sig = table2array(readtable('pos5Data'));
jawPossig = table2array(readtable('jawPosData'));
jawEffortsig = table2array(readtable('jawEffortData'));
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
network_suffix = '_10_10';
x = [pos2sig,Nsig];
load(sprintf('network%s.mat',network_suffix),'net','tr','mean_x','stddev_x','max_ati');
x = (x-mean_x)./stddev_x;
x = x';
%% Find calibration results
net = network(net);
y = net(x);
pred = y'.*max_ati;
end