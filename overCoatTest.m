clc;clear all;close all
%% add the location of the code
warning('off','all');
rmpath(genpath(sprintf('%s%s',pwd)));
warning('on','all');
addpath(sprintf('%s%s',pwd,'\20180819-overcoatTest'));
print_on = false;
plot_suffix = '_10_10_overcoat';
%% read data
ati = readForceData();
diffsig = table2array(readtable('diffData'));
sumsig = table2array(readtable('sumData'));
possig = table2array(readtable('posData'));
%% bias correction
diffsig = diffsig - mean(diffsig(1:500,:));
Nsig = diffsig./sumsig;
Nsig = [Nsig,Nsig.^2];
x = [possig(1:size(Nsig,1)),Nsig];
%% data sampling parameters
numPoint = size(ati,1);
Fs = 1500;
Ts = 1/Fs;
t = (1:numPoint)*Ts;
atiLables = {'F_{Nx}','F_{Ny}','F_{Nz}','M_{Nx}','M_{Ny}','M_{Nz}'};
ofsLabels = {'V_{N1}','V_{N2}','V_{N3}','V_{N4}','V_{N5}','V_{N6}'};
%%
load('network_10_10.mat','net','tr','mean_x','stddev_x','max_ati')
net = network(net);
x = (x-mean_x)./stddev_x;
x = x';
% Test the Network
y = net(x);

pred = y'.*max_ati;
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
close all;
figure(1);
% set(gcf,'units','normalized','OuterPosition',[0 0 1 1]);
set(gcf,'units','normalized','position',[0.4 0.2 0.3 0.6]);
counter = 1;
for axisNum=1:6
    switch axisNum
        case 1
            subplot(3,2,1);
            plot(t,ati(:,1),'b','linewidth',2);
            hold on;
            plot(t,pred(:,1),'r','linewidth',2);
            ylabel('F_x (N)','fontweight','b');
            legend('Cannula Outer Tube','Instrument','NumColumns',2);
        case 2
                        subplot(3,2,3);
            plot(t,ati(:,2),'b','linewidth',2);
            hold on;
            plot(t,pred(:,2),'r','linewidth',2);
            ylabel('F_y (N)','fontweight','b');
        case 3
                        subplot(3,2,5);
            plot(t,ati(:,3),'b','linewidth',2);
            hold on;
            plot(t,pred(:,3),'r','linewidth',2);
            ylabel('F_z (N)','fontweight','b');
        case 4
                        subplot(3,2,2);
            plot(t,ati(:,4)/1000,'b','linewidth',2);
            hold on;
            plot(t,pred(:,4)/1000,'r','linewidth',2);
            ylabel('M_x (N.m)','fontweight','b');
            
        case 5
                        subplot(3,2,4);
            plot(t,ati(:,5)/1000,'b','linewidth',2);
            hold on;
            plot(t,pred(:,5)/1000,'r','linewidth',2);
            ylabel('M_y (N.m)','fontweight','b');
            
        case 6
                        subplot(3,2,6);
            plot(t,ati(:,6)/1000,'b','linewidth',2);
            hold on;
            plot(t,pred(:,6)/1000,'r','linewidth',2);
            ylabel('M_z (N.m)','fontweight','b');
    end
    
    grid on;
    ax = gca;
    ax.FontSize = 18;
    set(gca,'linewidth',2)
    
    if axisNum == 3 || axisNum == 6
        set(gca,'XTickLabel',xticklabel);
        xlabel('Time(s)','fontweight','b');
    else
        xticklabel = get(gca,'XTickLabel');
        set(gca,'XTickLabel',[]);
    end
    
end
set(gcf,'color','w')
set(gcf, 'InvertHardCopy', 'off');

ha=get(gcf,'children');
set(ha(1),'position',[0.6105 0.16 0.2945 0.2157])
set(ha(2),'position',[0.6105 0.40 0.2945 0.2157])
set(ha(3),'position',[0.6105 0.64 0.2945 0.2157])
set(ha(4),'position',[0.1300 0.16 0.3347 0.2157])
set(ha(5),'position',[0.1300 0.40 0.3347 0.2157])
set(ha(6),'position',[0.1300 0.88 0.6105+0.2945-0.13 0.0359])
set(ha(7),'position',[0.1300 0.64 0.3347 0.2157])
if print_on
    print(sprintf('calibrationPlot%s',plot_suffix),'-djpeg','-r600');
end