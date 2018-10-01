function [] = Sept_mar_line_correlation_plots_DYN_regional_fig7()
   % Sept_mar_line_correlation_plots_LE()
   % 
   % plots correlations generated in September_lead_correlations_DYN_regional.m
   %
   % See figure 7 in "Processes controlling Arctic and Antarctic sea ice predictability 
   % in the Community Earth System Model" by Ordonez, Bitz, and Blanchard-Wrigglesworth 
   % (in press).
   %
   % Ana Ordonez 10/2018

   xlabel_1 = {'Sep','','Nov','','Jan','','Mar',...
                '','May','','Jul','','Sep'};
   xlabel_2 = {'Mar','','May','','Jul','','Sep',...
            '','Nov','','Jan','','Mar'};

   figure (1)
   res = 1;
   mask_name = 'pacific';

   subplot(2,2,1)
   nhmonth = 3;
   nhx = xlabel_2;
   fname = ['DYN_sept_nh_b.e11.B1850C5CN._',mask_name,'_',num2str(nhmonth),'.mat'];
   load(fname)
   lines(R,res);
   set(gca,'xtick',[0:12],'xticklabel',nhx)
   ylabel('correlation coefficient')
   %title('Pacific DYN')
   axis([0 12 0 1]);

   legend_list = {'area','volume','area > 0.64 m','area > 1.39 m',...
               'area > 2.47 m','multiyear area'};
   l = legend(legend_list,'location','northwest');
   set(l,'fontsize',8)

   subplot(2,2,2)
   nhmonth = 9;
   nhx = xlabel_1;
   fname = ['DYN_sept_nh_b.e11.B1850C5CN._',mask_name,'_',num2str(nhmonth),'.mat'];
   load(fname)
   lines(R,res);
   set(gca,'xtick',[0:12],'xticklabel',nhx)
   ylabel('correlation coefficient')
   %title('Pacific')
   axis([0 12 0 1]);

   mask_name = 'atlantic';

   subplot(2,2,3)
   nhmonth = 3;
   nhx = xlabel_2;
   fname = ['DYN_sept_nh_b.e11.B1850C5CN._',mask_name,'_',num2str(nhmonth),'.mat'];
   load(fname)
   lines(R,res);
   set(gca,'xtick',[0:12],'xticklabel',nhx)
   ylabel('correlation coefficient')
   %title('Atlantic DYN')
   xlabel('month')
   axis([0 12 0 1]);

   subplot(2,2,4)
   nhmonth = 9;
   nhx = xlabel_1;
   fname = ['DYN_sept_nh_b.e11.B1850C5CN._',mask_name,'_',num2str(nhmonth),'.mat'];
   load(fname)
   lines(R,res);
   set(gca,'xtick',[0:12],'xticklabel',nhx)
   ylabel('correlation coefficient')
   %title('Atlantic SOM')
   axis([0 12 0 1]);
   xlabel('month');


   print(gcf,['Line_correlation_DYN_sio.png'],'-dpng');
   close
end


function [] = lines(R,res)
    % Use the default 2014b colors except 
    % for the last value
    time = [0:12];
    plot(time, squeeze(R(1,res,:)), 'linewidth',2);
    hold
    plot(time, squeeze(R(2,res,:)), 'linewidth',2);
    plot(time, squeeze(R(3,res,:)), 'linewidth',2, 'linestyle', '--');
    plot(time, squeeze(R(4,res,:)), 'linewidth',2, 'linestyle', '--');
    plot(time, squeeze(R(5,res,:)), 'linewidth',2, 'linestyle', '--');
    %plot(time, squeeze(R(6,res,:)), 'linewidth',2, 'linestyle', '--');
    plot(time, squeeze(R(7,res,:)), 'linewidth',2, ...
                       'linestyle','-.','Color',[ 0.600    0.600    0.600]);
end
