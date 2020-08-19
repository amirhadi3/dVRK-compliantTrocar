function plot_Data(tplot,posplot_data,atiplot_data,Dplot_data,print_on,Rsq,plot_suffix)
eplot_data = atiplot_data-Dplot_data; 

% plot poistion data
figure;
[ha, pos] = tight_subplot(1, 1, 0.065, 0.1, 0.05);
set(gcf,'units','normalized','OuterPosition',[0 0.6 1 0.4]);
sfh = subplot(ha);
plot(tplot,posplot_data*1000,'b','linewidth',2);
grid on;
ax = gca;
ax.FontSize = 20;
set(gca,'linewidth',2)
grid on;
xlabel('Time(s)','fontweight','b');
ylabel('q_3(m)','fontweight','b');
set(gcf,'color','w')
set(gcf, 'InvertHardCopy', 'off');
if print_on
    print(sprintf('dispPlot%s',plot_suffix),'-djpeg','-r600');
end

figure;
forceLables = {'F_{x}(N)','F_{y}(N)','F_{z}(N)','M_{x}(N.mm)','M_{y}(N.mm)','M_{z}(N.mm)'};
errorLables = {'E_{Fx}','E_{Fy}','E_{Fz}','E_{Mx}','E_{My}','E_{Mz}'};
[ha, pos] = tight_subplot(6, 2, 0.065, 0.1, 0.05);
set(gcf,'units','normalized','OuterPosition',[0 0 1 1]);
index = reshape(1:12,2,6).';
for axisNum=1:12
    switch axisNum
        case {1,2}
            Dplot = Dplot_data(:,1);
            atiplot = atiplot_data(:,1);
            eplot = eplot_data(:,1);
        case {3,4}
            Dplot = Dplot_data(:,2);
            atiplot = atiplot_data(:,2);
            eplot = eplot_data(:,2);
        case {5,6}
            Dplot = Dplot_data(:,3);
            atiplot = atiplot_data(:,3);
            eplot = eplot_data(:,3);
        case {7,8}
            Dplot = Dplot_data(:,4);
            atiplot = atiplot_data(:,4);
            eplot = eplot_data(:,4);
        case {9,10}
            Dplot = Dplot_data(:,5);
            atiplot = atiplot_data(:,5);
            eplot = eplot_data(:,5);
        case {11,12}
            Dplot = Dplot_data(:,6);
            atiplot = atiplot_data(:,6);
            eplot = eplot_data(:,6);
    end
    
    sfh = subplot(ha(index(axisNum)));
    if (mod(index(axisNum),4) == 3 || mod(index(axisNum),4) == 0)
        sfh.Position = sfh.Position - [0 0.01 0 0.00];
        plot(tplot,atiplot-Dplot,'k','linewidth',2);
        ylabel(errorLables{ceil(axisNum/2)},'fontweight','b');
        hold on;
        set(gca,'XTickLabel',[]);
    else
        sfh.Position = sfh.Position + [0 -0.06 0 0.09];
        plot(tplot,atiplot-Dplot,'k','linewidth',2);
        hold on;
        plot(tplot,atiplot,'b','linewidth',2)
        plot(tplot,Dplot,'--r','linewidth',2)
        ylabel(forceLables{ceil(axisNum/2)},'fontweight','b');
        xticklabel = get(gca,'XTickLabel');
        set(gca,'XTickLabel',[]);
        ylim([-max(abs(atiplot)) max(abs(atiplot))]);
    end
    if axisNum == 1
        legend({'Error','ATI','OFS'},'Fontsize',16,'NumColumns',3,'location','southwest');
    end
    
    ax = gca;
    ax.FontSize = 20;
    set(gca,'linewidth',2)
    grid on;
    %     legend ATI OFS Error
    if axisNum == 6 || axisNum == 12
        set(gca,'XTickLabel',xticklabel);
        xlabel('Time(s)','fontweight','b');
    end
end
set(gcf,'color','w')
set(gcf, 'InvertHardCopy', 'off');
if print_on
    print(sprintf('calibrationPlot%s',plot_suffix),'-djpeg','-r600');
end

%% lineairty plot
titleList = {'F_x(N)','F_y(N)','F_z(N)','M_x(N.mm)','M_y(N.mm)','M_z(N.mm)'};
figure;
[ha, pos] = tight_subplot(1, 6, 0.028, 0.1, 0.05);
% set(gcf,'units','normalized','Position',[0.0031    0.0056    0.8078    0.9354]);
set(gcf,'units','normalized','OuterPosition',[0 0.55 1 0.45]);
for axisNum = 1:6
    sfh = subplot(ha(axisNum));
    Dplot = Dplot_data(:,axisNum);
    atiplot = atiplot_data(:,axisNum);
    xmax = max(abs(atiplot));
    ymax= max(abs(Dplot));
    auxVar = [ones(size(Dplot)) Dplot];
    linfit_Par = auxVar'*auxVar\(auxVar'*atiplot);
    limVal = max(abs(atiplot));
    xlinfit = (-xmax:xmax)';
    ylinfit = [ones(size(xlinfit)),xlinfit]*linfit_Par;
    plot(atiplot,Dplot,'b','linewidth',2);
    hold on
    plot(xlinfit,ylinfit,'--r','linewidth',2);
    
    a = get(gca,'XTick');
    set(gca,'YTick',a,'FontSize',14)
    set(gca,'YTick',a,'FontSize',14)
    
    xlim([-xmax xmax]);
    ylim([-ymax ymax]);
    
    grid on
    ax=gca;
    ax.FontSize = 20;
    set(gca,'linewidth',2);
    xlabel('ATI','Fontsize',20,'fontweight','b');
    if (axisNum == 1)
        ylabel('OFS','Fontsize',20,'fontweight','b');
    end
    txt = sprintf('R^2 = %4.3f',Rsq(axisNum));
    text(-xmax/2,ymax/2,txt,'Fontsize',18,'fontweight','b');
    title(titleList{axisNum},'Fontsize',24,'fontweight','b')
    axis square
end
set(gcf,'color','w')
set(gcf, 'InvertHardCopy', 'off');
if print_on
    print(sprintf('linearityPlot%s',plot_suffix),'-djpeg','-r600');
end
end