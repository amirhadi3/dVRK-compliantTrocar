clc;clear all;close all
%% add the location of the code
warning('off','all');
rmpath(genpath(sprintf('%s%s',pwd)));
warning('on','all');
print_on = false;
%% read lateral0, axial 0, lateral1, axial 1
loc = {'overcoatTest'};
ati = [];
diffsig = [];
sumsig = [];
pos = [];
index = [];
for i = 1:length(loc)
    ati = [ati;table2array(readtable(sprintf('%s%s',loc{i},'/atiData.txt')))];
    diffsig = [diffsig;table2array(readtable(sprintf('%s%s',loc{i},'/diffData')))];
    sumsig = [sumsig;table2array(readtable(sprintf('%s%s',loc{i},'/sumData')))];
    pos = [pos;table2array(readtable(sprintf('%s%s',loc{i},'/posData')))];
    index = [index;length(ati)];
end
%% bias correction
ati = ati(:,[1:2,4:6]);
ati =ati-mean(ati(1:500,:));
diffsig = diffsig - mean(diffsig(1:500,:));
Nsig = diffsig./sumsig;
Nsig = [Nsig,Nsig.^2];
possig =ones(length(Nsig),2)*0.08;
possig(:,2) = 0;
x = [possig,Nsig];
%% data sampling parameters
numPoint = size(ati,1);
Fs = 1500;
Ts = 1/Fs;
t = (1:numPoint)*Ts;
%%
load('network_noGrip_5.mat','net','tr','mean_x','stddev_x','max_ati')
net = network(net);
x = (x-mean_x)./stddev_x;
x = x';
% Test the Network
y = net(x);

pred = y'.*max_ati;
pred = pred-mean(pred(1:500,:));
%% plot Results
atiplot_data = ati;
Dplot_data = pred;
tplot = (1:size(atiplot_data,1))*Ts;
eplot_data = atiplot_data-Dplot_data;
posplot_data = possig;
%% plotData
close all;
figure(1);
% set(gcf,'units','normalized','OuterPosition',[0 0 1 1]);
set(gcf,'units','normalized','position',[0.1 0.1 0.4 0.8]);
counter = 1;
for axisNum=1:6
    switch axisNum
        case 1
            subplot(5,1,1);
            plot(t,ati(:,1),'b','linewidth',2);
            hold on;
            plot(t,pred(:,1),'r','linewidth',2);
            ylabel('f_x(N)','fontweight','b');
            legend('Cannula Outer Tube','Instrument','NumColumns',2);
        case 2
            subplot(5,1,2);
            plot(t,ati(:,2),'b','linewidth',2);
            hold on;
            plot(t,pred(:,2),'r','linewidth',2);
            ylabel('f_y(N)','fontweight','b');
        case 3
%             handle = subplot(3,2,5);
%             plot(t,ati(:,3),'b','linewidth',2);
%             hold on;
%             plot(t,pred(:,3)/5,'r','linewidth',2);
%             ylabel('f_z (N)','fontweight','b');
        case 4
            subplot(5,1,3);
            plot(t,ati(:,3),'b','linewidth',2);
            hold on;
            plot(t,pred(:,3),'r','linewidth',2);
            ylabel('m_x(N.m)','fontweight','b');
            
        case 5
            subplot(5,1,4);
            plot(t,ati(:,4),'b','linewidth',2);
            hold on;
            plot(t,pred(:,4),'r','linewidth',2);
            ylabel('m_y(N.mm)','fontweight','b');
            
        case 6
            subplot(5,1,5);
            plot(t,ati(:,5),'b','linewidth',2);
            hold on;
            plot(t,pred(:,5),'r','linewidth',2);
            ylabel('m_z(N.mm)','fontweight','b');
    end
    
%     if axisNum == 3
%         delete(handle)
%     end
    
    grid on;
    ax = gca;
    ax.FontSize = 14;
    set(gca,'linewidth',2)
    
    if axisNum == 6
        set(gca,'XTickLabel',xticklabel);
        xlabel('Time(s)','fontweight','b');
    else
        xticklabel = get(gca,'XTickLabel');
        set(gca,'XTickLabel',[]);
    end
    
end
set(gcf,'color','w')
set(gcf, 'InvertHardCopy', 'off');

% ha=get(gcf,'children');
% set(ha(1),'position',[0.6105 0.15 0.2945 0.17])
% set(ha(2),'position',[0.6105 0.35 0.2945 0.17])
% set(ha(3),'position',[0.6105 0.55 0.2945 0.17])
% set(ha(4),'position',[0.1300 0.15 0.3347 0.17])
% set(ha(5),'position',[0.1300 0.35 0.3347 0.17])
% set(ha(6),'position',[0.1300 0.73 0.6105+0.2945-0.13 0.0359])
% set(ha(7),'position',[0.1300 0.55 0.3347 0.17])
if print_on
    print('overcoat','-djpeg','-r600');
end