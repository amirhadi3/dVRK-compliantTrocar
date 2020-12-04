clc;clear all;close all
%% add the location of the code
print_on = false;
proGraspMomentArm = 0.017; %m
%% read data
pinchLis = [0,1,2,3];
for i=1:4
%     [t{i},possig{i},effsig{i},pred{i}] = estimate(sprintf('data_shaft/pinch%d',i-1),'network_noGrip_5','NoGrip');
    [t{i},possig{i},effsig{i},pred{i}] = estimate(sprintf('data_shaft/pinch%d',i-1),'network_wGrip_5','wGrip');
end
%%
figure(1);
newcolors = ['#0072BD'; '#D95319'; '#7E2F8E';'#77AC30'];
colororder(newcolors)
set(gcf,'units','normalized','position',[0.1 0.1 0.4 0.8]);
for i = 1:4
    h{1} = subplot(4,2,1);
    plot(t{i},possig{i}(:,end)*180/pi,'linewidth',2);hold on;
    ylabel('q_7(deg.)','fontweight','b');
    l = legend('no-foam','1-layer','2-layers','3-layers','NumColumns',4,'fontsize',12);
    set(l,'position',[0.17    0.95    0.7    0.03]);
    setplot(gca)
    
    h{2}=subplot(4,2,2);
    plot(t{i},effsig{i}(:,end)/proGraspMomentArm,'linewidth',2);hold on;
    ylabel('f_{jaw}(N)','fontweight','b');
    setplot(gca)
    
    h{3} = subplot(4,2,3);
    plot(t{i},pred{i}(:,1),'linewidth',2);hold on;
    ylabel('f_x(N)','fontweight','b');
    setplot(gca)
    
    h{4} = subplot(4,2,4);
    plot(t{i},pred{i}(:,3),'linewidth',2);hold on;
    ylabel('m_x(N.mm)','fontweight','b');
    setplot(gca)
    
    h{5} = subplot(4,2,5);
    plot(t{i},pred{i}(:,2),'linewidth',2);hold on;
    ylabel('f_y(N)','fontweight','b');
    xlabel('time (s)','fontweight','b')
    setplot(gca)
    
    h{6} = subplot(4,2,6);
    plot(t{i},pred{i}(:,4),'linewidth',2);hold on;
    ylabel('m_y(N.mm)','fontweight','b');
    xlabel('time (s)','fontweight','b');
    setplot(gca)
    
    h{7} = subplot(4,2,[7,8]);
    plot(t{i},pred{i}(:,5),'linewidth',2);hold on;
    ylabel('m_z(N.mm)','fontweight','b');
    xlabel('time (s)','fontweight','b')
    setplot(gca)
    %     print('pinchtest','-djpeg','-r600');
end

function setplot(h)
hold on;
grid on;
set(h,'linewidth',2);
set(h,'FontSize',12);
end
%%
function [t,pos,eff,pred] = estimate(file,nname,option)
[pos,eff,x] = readdata(file, option);
load(sprintf('%s.mat',nname),'net','tr','mean_x','stddev_x','max_ati');
x = (x-mean_x)./stddev_x;
%% Find calibration results
net = network(net);
y = net(x');
pred = y'.*max_ati;
t = (1:length(pred))/1000;
% if strcmp(option,'wGrip')
%     pred = pred-mean(pred(1:10,:));
%     pred(:,1) = pred(:,1)/10;
% end
end

function [pos,eff,x] = readdata(file, option)
warning('off','all');
rmpath(genpath(sprintf('%s%s',pwd)));
warning('on','all');
addpath(sprintf('%s\\%s',pwd,file));
diffsig = table2array(readtable('diffData'));
sumsig = table2array(readtable('sumData'));
pos = table2array(readtable('posData'));
eff = table2array(readtable('effData'));
%% bias correction
diffsig = diffsig-mean(diffsig(1:10,:));
possig = pos(:,[3,4]);
effsig = eff(:,end);
Nsig = diffsig./sumsig;
Nsig = [Nsig,Nsig.^2];
if strcmp(option,'wGrip')
    x = [possig,effsig,Nsig];
elseif strcmp(option,'NoGrip')
    x = [possig,Nsig];
else
    error('the option is invalid')
end
end