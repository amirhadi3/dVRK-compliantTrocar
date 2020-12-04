function plot_Data5(tplot,posplot_data,atiplot_data,Dplot_data,print_on,Rsq,plot_suffix)
eplot_data = atiplot_data-Dplot_data;

figure;
forceLables = {'f_{x}(N)','f_{y}(N)','f_{z}(N)','m_{x}(N.mm)','m_{y}(N.mm)','m_{z}(N.mm)'};
errorLables = {'e_{fx}','e_{fy}','e_{fz}','e_{mx}','e_{my}','e_{mz}'};
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
            %             Dplot = Dplot_data(:,3);
            %             atiplot = atiplot_data(:,3);
            %             eplot = eplot_data(:,3);
        case {7,8}
            Dplot = Dplot_data(:,3);
            atiplot = atiplot_data(:,3);
            eplot = eplot_data(:,3);
        case {9,10}
            Dplot = Dplot_data(:,4);
            atiplot = atiplot_data(:,4);
            eplot = eplot_data(:,4);
        case {11,12}
            Dplot = Dplot_data(:,5);
            atiplot = atiplot_data(:,5);
            eplot = eplot_data(:,5);
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
    set(gca,'linewidth',2)
    grid on;
    %     legend ATI OFS Error
    if axisNum == 4 || axisNum == 12
        set(gca,'XTickLabel',xticklabel);
        xlabel('Time(s)','fontweight','b');
    end
    ax.FontSize = 14;
    if axisNum == 5 || axisNum == 6
        delete(sfh)
    end
end
set(gcf,'color','w')
set(gcf, 'InvertHardCopy', 'off');
if print_on
    print(sprintf('calibrationPlot%s',plot_suffix),'-djpeg','-r600');
end
end