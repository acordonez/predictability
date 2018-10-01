% lag_correlations_aice_v4a
%
% Calculates monthly lag correlations for sea ice area and sea ice volume in the
% CESM large ensemble control runs. This analysis is done for different
% regions in the Arctic and Antarctic.
%
% This script can take an hour or more to run depending on the 
% number of regions and variables being plotted.
%
% See figure 4 in "Processes controlling Arctic and Antarctic sea ice predictability 
% in the Community Earth System Model" by Ordonez, Bitz, and Blanchard-Wrigglesworth (in press).
%
% Ana Ordonez 10/2018


for region = 1:2

   var_list = {'area','volume','volume-area'};
   if region == 1
      % Antarctica
       thehemi = 'sh';
       mod_ind_upper = 3;
   else
       % Arctic
       thehemi = 'nh';
       mod_ind_upper = 2;
   end
   disp(thehemi);

   % loop over DYN and SOM
   for mod_ind = 1:2 

      area_map = {'tx0.1v2','gx1v6'};
      thestart=[1 1];
      % loop over variables to correlate: 
      % 1=area->area, 2=volume->volume, 3=area->volume
      for var_ind = 1:3
         close all
         var = var_list{var_ind};
         map_stride = [1 1];
         % initialize correlation arrays
         R(12,13)=0;
         P=zeros(size(R));

         clear area hi_detrend hi_detrend2 tmp tot_area trnd 

         model_list = {'b.e11.B1850C5CN.','e.e11.E1850C5CN.'};
         model = model_list{res};
         if strcmp(model,'e.e11.E1850C5CN.')
            thesuff = '.090001-100112';
            first_year{res} = '0900';
            model_num='001';
         else 
             thesuff = '.210001-220012';
             first_year{res} = '2100';
             model_num ='005';
         end
         map_stride = [map_stride, 1];
         nyrs = 100;
         eq = [inf inf];
         dir = '/glade/p_old/cesmLE/CESM-CAM5-BGC-LE/ice/proc/tseries/monthly/';
         nc1 = ['f09_g16.',model_num,'.cice.h.'];
         nc2 = ['_',thehemi,thesuff,'.nc'];
         fh = [dir,'hi','/',model,nc1,'hi',nc2];
         fa = [dir,'aice','/',model,nc1,'aice',nc2];
         area = ncread(fh,'tarea',[1 1],[inf eq(res)],...
                       [map_stride(res) map_stride(res)]); 
         area = repmat(area,1,1,nyrs);
         lon = ncread(fa,'TLON',[1 1],[inf eq(res)],...
                      [map_stride(res) map_stride(res)]);
         lon = repmat(lon,1,1,nyrs);
         lat = ncread(fa,'TLAT',[1 1],[inf eq(res)],...
                      [map_stride(res) map_stride(res)]);
         lat = repmat(lat,1,1,nyrs);
         add_yr = 0*12;

         % Set region masks
         lon = wrapTo180(lon);
         seamask = zeros(size(area));
         if region == 1
            seamask(lon >= -60 & lon < 20) = 1;
            seamask(lon >= 20 & lon < 90) = 2;
            seamask(lon >= 90 & lon < 160) = 3;
            seamask(lon >= 160) = 4;
            seamask(lon < -130) = 4;
            seamask(lon >= -130 & lon < -60) = 5;
            ocean_list = {'Weddell','Indian','WPacific','Ross','ABSeas'};
   
         else
            seamask(lon >= -45 & lon < 20 & lat < 80) = 1;% Greenland
            seamask(lon >= 20 & lon < 100 & lat < 80) = 2;%Kara & barents
            seamask(lon >= 90 & lat < 65 | lon < -135 & lat < 65) = 3;%Bering sea and Seas of Okhotsk and japan
            seamask(lon < -65 & lon >= -100 & lat < 70) = 4;% Hudson Bay
            seamask(seamask ~= 4 & lat < 80 & lat > 40 & lon < -45 & lon > -90) = 5;% Baffin, Labrador, St. Lawrence
            seamask(seamask ==0 & lat > 65 & lon >90) = 6;% Arctic
            seamask(seamask == 0 & lat > 65 & lon < -90) = 6;% also arctic
            seamask(seamask == 0 & lat > 70 & lon <= 90 & lon >= -90) = 6;% Arctic 
            ocean_list = {'Greenland','BarentsKara','BerOkhJap',...
                             'Hudson','BafLabLaw','Arctic'};
         end

         % Loop over all the regions and load area
         for ocean = 1:length(ocean_list)
            R=zeros(size(R));
	    P=zeros(size(P));
            ocean_name = ocean_list{ocean};

            for mo=1:12     
            % get monthly time series
               if (strcmp(var,'area') | strcmp(var,'volume-area'))
                  tmp= ncread(fa,'aice',[thestart mo+add_yr],[inf eq(res) nyrs],...
                              [map_stride(res) map_stride(res) 12]);
                  tmp(tmp==0) = NaN;
                  tmp(seamask ~= ocean) = NaN;
                  tmp = tmp .* area ./ 100; % ice area
                  tmp = squeeze(nansum(nansum(tmp))); %total ice area
                  if strcmp(var,'volume-area')
                     tmp2 = tmp; clear tmp;
                  end
               end
               if (strcmp(var,'volume') | strcmp(var,'volume-area'))
                   tmp= ncread(fh,'hi',[thestart mo+add_yr],[inf eq(res) nyrs],...
                               [map_stride(res) map_stride(res) 12]);
                   tmp(tmp == 0) = NaN;
                   tmp(seamask ~= ocean) = NaN;
                   tmp = tmp .* area; % volume/grid area .* grid area = volume
                   tmp = squeeze(nansum(nansum(tmp,2),1)); % add after weighing
               end
               
               detrend before correlating
               trnd = gettrend(squeeze(tmp));
               hi_detrend(res,mo,:) = tmp'-[0:nyrs-1].*trnd./length(tmp);
               
               if strcmp(var,'volume-area')
                  trnd = gettrend(squeeze(tmp2));
                  hi_detrend2(res,mo,:) = tmp2'-[0:nyrs-1].*trnd./length(tmp2);
               end
              
            end
           
            set variables to correlate
            var1 = hi_detrend;
            var2 = hi_detrend;
            if strcmp(var,'volume-area')
               var2 = hi_detrend2;
            end
    
            compute correlations
            for mo = 1:12
              for step = 0:12
                  if mo+step > 12
                     [r,p] = corrcoef(var1(res,mo,1:nyrs-1),var2(res,mo+step-12,2:nyrs),...
                                      'rows','pairwise');
                      R(mo,step+1) = r(1,2);
                  elseif step == 12
                     [r,p] = corrcoef(var1(res,mo,1:nyrs-1),var2(res,mo,2:nyrs),...
                                      'rows','pairwise');
                      R(mo,step+1) = r(1,2);                    
                  else
                     [r,p] = corrcoef(var1(res,mo,1:nyrs-1),var2(res,mo+step,1:nyrs-1),...
                                         'rows','pairwise');
                      R(mo,step+1) = r(1,2);
                  end
                  Z = 0.5 .* log((1+r(1,2)) / (1-r(1,2)));
                  % 95% confidence
                  stand = 1.96 .*(1 ./ sqrt((nyrs-1) - 3));
                  % 99% confidence
                  %stand = 2.58 .*(1 ./ sqrt((nyrs-1) - 3));
                  if Z > 0
                     conf = Z - stand;
                  else
                     conf = Z + stand;
                  end
                     conf = tanh(conf);
                  if abs(r(1,2)-conf) < abs(r(1,2)) 
                     P(mo,step+1) = 1;
                  end
               end %step
            end %mo
    
            % Make correlation plot for this region
            figure (1)
            coordx=[0:12];
            coordy=[1:12];
            x,y]=meshgrid(coordx,coordy);
            Rplot=squeeze(R);
            Pplot=squeeze(P);
            months_sp = {'Jan','Feb','Mar','Apr','May','Jun','Jul','Aug','Sep','Oct','Nov','Dec'};
            if strcmp(thehemi,'sh')
               Rplot=Rplot([7:12,1:6],:);
               Pplot=Pplot([7:12,1:6],:);
               months_sp={'Jul','Aug','Sep','Oct','Nov','Dec','Jan','Feb','Mar','Apr','May','Jun'};
            end
            imagesc([0:1:12],[1:12],Rplot);
            hold
            k=find(Pplot==1); disp(size(k));
            plot(x(Pplot==1),y(Pplot==1),'k.','MarkerSize',20) 
            h = colorbar;
            caxis([-1 1]);
            colormap(flipud(lbmap(20,'RedBlue')));
            thelag = {'0','1','2','3','4','5','6','7','8','9','10','11','12'};
            set(gca,'xtick',coordx);
            set(gca,'ytick',coordy);
            set(gca,'yticklabel',months_sp,'fontsize',23);
            set(gca,'xticklabel',thelag,'fontsize',23);  
            title([model,' ',var,' ',ocean_name],'fontsize',23);
            print(gcf,['correlationplot_',ocean_name,'_',var,'_',model,'_',thehemi,'.png'],'-dpng','-r200');
            close

            end  %ocean
         end  %res
      end %var_ind
   end %mod_ind
end %region

