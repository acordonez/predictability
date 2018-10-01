% ice_area_hemi_seasonal_plot.m
%
% Makes monthly mean plots for ice extent and total volume
% for both hemispheres
%
% See figure 1 in "Processes controlling Arctic and Antarctic sea ice predictability 
% in the Community Earth System Model" by Ordonez, Bitz, and Blanchard-Wrigglesworth (in press).
%
% Ana Ordonez 10/2018

for sh = 0:1
   num = {'01','02','03','04','05','06','07','08','09','10','11','12'};

   % Change variables here for different plots
   vars = {'aice','hi'};
   plot_hi = 1; 

   disp(sh)

   if plot_hi
      maxind = 2;
   else
      maxind = 1;   
   end

   for theind = 1:maxind
      thevar = vars{theind};
      disp(thevar)

      % CESM seasonal cycle
      if sh
         hemi_str = 'sh';
         thestart = 1;
         theend = 76;
      else
         hemi_str = 'nh';
         thestart = 281;
         theend = 384;
      end
      ncf = ['/glade/p_old/cesmLE/CESM-CAM5-BGC-LE/ice/proc/tseries/monthly/',thevar];
      thesuff1 = {'.050001-059912','.060001-069912','.070001-079912','.080001-089912','.090001-100112'};
      thesuff2 = {'.170001-179912','.180001-189912','.190001-199912','.200001-209912','.210001-220012'};
      tmp = [];
      for n=1:length(thesuff1)
         % DYN
         f1 = ['/b.e11.B1850C5CN.f09_g16.005.cice.h.',thevar,'_',hemi_str,thesuff2{n},'.nc']; 
         tmp1 = ncread([ncf,f1],thevar,[1 1 1],[inf inf 100*12]); 
         tmp = cat(3,tmp,tmp1);
      end
      area_g = ncread([ncf,f1],'tarea');
      lim = size(tmp,2);

      tmp=reshape(tmp,320,lim,12,size(tmp,3)./12);
      tmp=nanmean(tmp,4);
      if (strcmp(thevar,'aice'))
         tmp(tmp < 15.) = 0;
         tmp(tmp >= 15.) = 1;
      end
      tmp=tmp .* repmat(area_g,1,1,12);
      tmp(tmp==0) = NaN;
      LE_dyn = squeeze(nansum(reshape(tmp,lim*320,12)));

      clear tmp; tmp = [];
      for n = 1:length(thesuff1)
         % SOM
         f4 = ['/e.e11.E1850C5CN.f09_g16.001.cice.h.',thevar,'_',hemi_str,thesuff1{n},'.nc']; 
         tmp1 = ncread([ncf,f4],thevar,[1 1 1],[inf inf 100*12]);
         tmp = cat(3,tmp,tmp1);
      end
      lim = size(tmp,2);
      tmp=reshape(tmp,320,lim,12,size(tmp,3)./12);
      tmp=nanmean(tmp,4);
      if (strcmp(thevar,'aice'))
         tmp(tmp < 15.) = 0;
         tmp(tmp >= 15.) = 1; 
      end
      tmp=tmp .* repmat(area_g,1,1,12);
      tmp(tmp==0) = NaN;
      LE_som = squeeze(nansum(reshape(tmp,lim*320,12)));

      % observations
      if sh 
         ice_obs = get_regional_monthly_sie_from_obs('s');
      else 
         ice_obs = get_regional_monthly_sie_from_obs('n');
      end
      ice_obs = nansum(ice_obs,1);

      % plot
      figure (sh+1)
      if strcmp(thevar,'aice')
         exp = 6; % for converting aice to km^2
         pltnum = 1;
      else
         exp = 0; % keep volume in m
         pltnum = 2;
      end

      (2,1,pltnum)
      plot(1:12,LE_dyn./10^exp,'color','black','linestyle','-','linewidth',2); 
      hold;
      plot(1:12,LE_som./10^exp,'color','black','linestyle','--','linewidth',2);
      if (strcmp(thevar,'aice'))
         plot(1:12,ice_obs,'color','black','linestyle',':','linewidth',2);
      end

      if sh & pltnum == 1
         h=legend('DYN','SOM','observations','location','northwest');
      elseif pltnum == 1
         h=legend('DYN','SOM','observations','location','southwest');
      end

      fsize = 14;

      if strcmp(thevar,'aice')
         ylabel('10^6 km^2','fontsize',fsize);
         set(gca,'xtick',[2:2:12]);
         if sh 
            axis([1 12 0.2E7 2.4E7]);
            set(gca,'TickDir','out','YMinorTick','on','ytick',[0.2E7:0.2E7:2.4E7]);
            set(gca,'ticklength',[0.02 0.02])
            labels = {'2','','6','','10','','14','','18','','22',''};
            set(gca,'yticklabel',labels,'fontsize',fsize);
            %title('Antarctic monthly mean ice area','fontsize',fsize);
         else
            axis([1 12 4.0E6 1.8E7]);
            set(gca,'TickDir','out','YMinorTick','on','ytick',[6E6:2E6:20E6]);
            set(gca,'TickLength',[0.02 0.02])
            labels = {'4','','8','','12','','16',''};
            set(gca,'yticklabel',[6:2:20],'fontsize',fsize)
            %title('Arctic monthly mean ice area','fontsize',fsize);
         end
         set(gca,'XMinorTick','on','xticklabel',[2:2:12],'fontsize',fsize);
         set(h,'FontSize',14);
         set(gca,'fontsize',fsize);
         title('Extent')
      else
         %title('Monthly mean total sea ice volume','fontsize',fsize);
         ylabel('10^1^3 m^3','fontsize',fsize)   
         if sh
            axis([1 12 0 2.5E13])
            set(gca,'TickDir','out','YMinorTick','on','ytick',[0E13:0.5E13:2.5E13]);
            set(gca,'TickLength',[0.02 0.02])
            set(gca,'yticklabel',[0:0.5:2.5],'fontsize',fsize);
         else
            axis([1 12 2E13 4.5E13]);
            set(gca,'TickDir','out','YMinorTick','on','ytick',[2E13:0.5E13:4.5E13]);
            set(gca,'TickLength',[0.02 0.02])
            set(gca,'yticklabel',[2:0.5:4.5],'fontsize',fsize);
         end
         set(gca,'XMinorTick','on','xtick',[2:2:12]);
         set(gca,'xticklabel',[2:2:12],'fontsize',fsize);
         xlabel('month');
         set(h,'FontSize',14);
         title('Volume')
      end
      set(gcf,'PaperPositionMode','manual')
      set(gcf, 'PaperUnits', 'inches');
      set(gcf, 'Position', [0 0 2.5 5]);
      print(gcf,['ice_seasonal',hemi_str,'.png'],'-dpng');
   end
end
