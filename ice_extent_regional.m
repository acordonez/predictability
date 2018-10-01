% ice_extent_regional.m
% 
% Plots the ice extent in CESM and observations for regions
% in the Arctic and Antarctic. 
% This version corrects for differences in land masks between grids.
%
% See figure 2 in "Processes controlling Arctic and Antarctic sea ice predictability 
% in the Community Earth System Model" by Ordonez, Bitz, and Blanchard-Wrigglesworth (in press).
%
% Ana Ordonez 10/2018

% line styles for plot
hemilinestyle = {'-k','--k'};
hemidotstyle = {'k.','kx'};

% loop through hemispheres
for h = 1:2
   if h == 1
      thehemi='nh';
      zone = [304,448];
   elseif h==2
      thehemi='sh';
      zone = [316,332];
   end
   % loads area data from observation grid
   datdir = '/glade/work/aordonez/';
   area_id = fopen([datdir,'sat_ice/ps',thehemi(1),'25area_v3.dat'],'r','l');
   area_obs = fread(area_id,zone,'int') ./ 1000; %sq km
   fclose(area_id);
   % look through DYN and SOM
   for res=1:2;
      disp('setting model data')

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
      nmnths = 1200; 
      thestart = [1 1];
      eq = [inf inf];
      dir = '/glade/p_old/cesmLE/CESM-CAM5-BGC-LE/ice/proc/tseries/monthly/';
      nc1 = ['f09_g16.',model_num,'.cice.h.'];
      nc2 = ['_',thehemi,thesuff,'.nc'];

      % read LE data
      if res == 1
          thesuff = {'.170001-179912','.180001-189912','.190001-199912',...
                     '.200001-209912','.210001-220012'};
      else
           thesuff = {'.050001-059912','.060001-069912','.070001-079912',...
                      '.080001-089912','.090001-100112'};
      end
      tmp = [];
      for n=1:length(thesuff)
          nc2 = ['_',thehemi,thesuff{n},'.nc'];
          f1 = [dir,'aice','/',model,nc1,'aice',nc2];
          tmp1 = ncread(f1,'aice',[thestart 1],[inf eq(res) nmnths]); 
          tmp = cat(3,tmp,tmp1); 
      end

      if strcmp(thehemi,'nh')
            ocean_list = {'N Pacific','Central Arctic','Barents & Kara',...
                 'Hudson Bay','NW Atlantic','GIN Seas'};
      else
            ocean_list = {'Weddell','Indian','W Pacific','Ross','AB Seas'};
      end

      disp('loading obs data')
      if h == 2 
         [ice,lat,lon] = get_monthly_sie_from_obs('s');
         [seamask,~] = get_seamask(['sh'],lat,lon);
         seamask2 = repmat(seamask,1,1,39);
         ice_obs = ice ;
      else
         [ice,lat,lon] = get_monthly_sie_from_obs('n');
         [seamask,~] = get_seamask(['nh'],lat,lon);
         seamask2 = repmat(seamask,1,1,39);
         ice_obs = ice;
      end

      disp('computing masked extent')

      % put LE data onto obs higher res grid
      le_regrid = zeros(size(ice,1),size(ice,2),size(tmp,3));
      if h == 2
          thefile = 'map_gx1v6SH_TO_SHstereo25km_blin.161213.nc';
      else
          thefile  = 'map_gx1v6NH_TO_stereo25km_blin.161123.nc';
      end
      for time = 1:size(tmp,3)
          le_regrid(:,:,time) = grid1togrid2(tmp(:,:,time)',thefile);
      end

      % Choose LE grid cells that have extent
      le_regrid(le_regrid < 15) = 0;
      le_regrid(le_regrid >=15) = 1;   
 
      % only count areas where both LE and obs have ice not land
      ice_obs(isnan(repmat(le_regrid(:,:,1),1,1,size(ice_obs,3)))) = NaN;

      % create masks that match size of each dataset
      seamask1 = repmat(seamask,1,1,size(le_regrid,3));
      seamask2 = repmat(seamask,1,1,size(ice_obs,3));

      disp('plotting figure')
      % set up figure
      if res == 1
         figure (1)
         fig = gcf;
         fig.PaperUnits = 'inches';
         fig.PaperPosition = [0 0 7 4] 
         fig.PaperPositionMode = 'manual';
      end

      % loop through each ocean basin
      % for totaling extent and plotting
      for ocean=1:length(ocean_list);
         %LE
         tmp = le_regrid;
         tmp(seamask1 ~= ocean) = 0;
         tmp = tmp .* repmat(area_obs,1,1,size(tmp,3)); %area x extent flag
         tmp=squeeze(nansum(nansum(tmp,2),1)); % total timeseries
         tmp=reshape(tmp,12,nmnths*5./12); % reshape to monthly
         tmp=nanmean(tmp,2); % monthy mean
         %Observations
         tmp2 = ice_obs;
         tmp2(seamask2 ~= ocean) = 0;
         tmp2 = reshape(tmp2,size(tmp2,1),size(tmp2,2),12,size(tmp2,3)/12);
         tmp2 = nansum(nansum(tmp2,1),2);
         tmp2 = squeeze(nanmean(tmp2,4));

         %tmp=tmp./1E12;
         subplot(2,3,ocean)
         plot(tmp./1e6,hemilinestyle{res},'linewidth',1.5);
         if res == 1
             hold
         end
         plot(tmp./1e6,hemidotstyle{res},'markersize',12);  
         if res == 2
             plot(tmp2./1e6,'k:','linewidth',2);
         end
         xlim([1 12]);
         title(ocean_list{ocean});
      end
   clear le_regrid ice_obs ice tmp tmp2
   end
  %legend(ocean_list)
  print(gcf,['DYN_and_SOM_ext_LE_',model_list{res},'_',thehemi,'500yrs_regions.png'],'-dpng','-r300');
  close
  clear lon lat area seamask 
 
end
