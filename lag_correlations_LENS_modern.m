% lag_correlations_LENS_modern
%
% This script computes and plots lag correlations for the cesm large ensemble
% years 1979-2009.
%
% Ana Ordonez 10/2018

LE = 1;

for region = 1:2
    disp(region)
    var_list = {'area','volume','volume-area'};
    if region == 1
       % Antarctica
       thehemi = 'sh';
       mod_ind_upper = 3;
       oceanlen = 6;
    else
       % Arctic
       thehemi = 'nh';
       mod_ind_upper = 2;
       oceanlen = 7;
    end
    disp(thehemi);
    thestart=[1 1];

    for ocean = 1:oceanlen
        disp('in ocean')
        for var_ind = 1:3
            disp(['var = ',var_list{var_ind}])
            R =zeros(32,12,13);
            P=zeros(size(R));
            for mod_num = 2:33
                disp(['model_num = ',num2str(mod_num)])
                var = var_list{var_ind};
                map_stride = [1 1];
                resmax=1;
                res = 1;
                clear area hi_detrend hi_detrend2 tmp tot_area trnd 
                model1 = 'b.e11.B20TRC5CNBDRD.';
                model2 = 'b.e11.BRCP85C5CNBDRD.';
                thesuff1 = '.192001-200512';
                thesuff2 = '.200601-208012';
                if mod_num < 10
                   model_num =['00',num2str(mod_num)];
                else
                   model_num = ['0',num2str(mod_num)];
                end
                map_stride = [map_stride, 1];
                nyrs1 = 27;
                nyrs2 = 3;
                eq = [inf inf];
                dir = '/glade/p/cesm0005/CESM-CAM5-BGC-LE/ice/proc/tseries/monthly/';
                nc1 = ['f09_g16.',model_num,'.cice.h.'];
                nc2 = ['_',thehemi,thesuff1,'.nc'];
                nc22 = ['_',thehemi,thesuff2,'.nc'];
                fh1 = [dir,'hi','/',model1,nc1,'hi',nc2];
                fh2 = [dir,'hi','/',model2,nc1,'hi',nc22];
                fa1 = [dir,'aice','/',model1,nc1,'aice',nc2];
                fa2 = [dir,'aice','/',model2,nc1,'aice',nc22];
                area = ncread(fh1,'tarea',[1 1],[inf eq(res)],...
                              [map_stride(res) map_stride(res)]); 
                area = repmat(area,1,1,30);
                lon = ncread(fa1,'TLON',[1 1],[inf eq(res)],...
                             [map_stride(res) map_stride(res)]);
                lon = repmat(lon,1,1,30);
                lat = ncread(fa1,'TLAT',[1 1],[inf eq(res)],...
                              [map_stride(res) map_stride(res)]);
                lat = repmat(lat,1,1,30);
                add_yr = 58*12; %start in year 1979
            
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
                if ocean > length(ocean_list)
                   seamask(seamask < ocean) = ocean;
                   if region == 1
                       ocean_name = 'SH';
                   else
                       ocean_name = 'NH';
                   end
                else
                   ocean_name = ocean_list{ocean};
                end
                for mo=1:12
                    % get time series for each month separately
                    if (strcmp(var,'area') | strcmp(var,'volume-area'))
                        tmp= ncread(fa1,'aice',[thestart mo+add_yr],[inf eq(res) nyrs1],...
                                    [map_stride(res) map_stride(res) 12]);
                        tmpcat= ncread(fa2,'aice',[thestart mo],[inf eq(res) nyrs2],...
                                    [map_stride(res) map_stride(res) 12]);
                        tmp = cat(3,tmp,tmpcat);
                        tmp(tmp==0) = NaN;
                        tmp(seamask ~= ocean) = NaN;
                        tmp = tmp .* area ./ 100; % ice area
                        tmp = squeeze(nansum(nansum(tmp))); %total ice area
                        if strcmp(var,'volume-area')
                           tmp2 = tmp; clear tmp;
                        end
                     end
                     if (strcmp(var,'volume') | strcmp(var,'volume-area'))
                        tmp= ncread(fh1,'hi',[thestart mo+add_yr],[inf eq(res) nyrs1],...
                                    [map_stride(res) map_stride(res) 12]);
                        tmpcat= ncread(fh2,'hi',[thestart mo],[inf eq(res) nyrs2],...
                                    [map_stride(res) map_stride(res) 12]);
                        tmp = cat(3,tmp,tmpcat);
                        tmp(tmp == 0) = NaN;
                        tmp(seamask ~= ocean) = NaN;
                        tmp = tmp .* area; % volume/grid area .* grid area = volume
                        tmp = squeeze(nansum(nansum(tmp,2),1)); % add after weighing
                     end
               
                     % detrend before correlating
                     nyrs = 30;
                     trnd = gettrend(squeeze(tmp));
                     hi_detrend(res,mo,:) = tmp'-[0:30-1].*trnd./length(tmp);
               
                     if strcmp(var,'volume-area')
                         trnd = gettrend(squeeze(tmp2));
                         hi_detrend2(res,mo,:) = tmp2'-[0:nyrs-1].*trnd./length(tmp2);
                     end
              
                 end % monthly time series
           
                 % set variables to correlate
                 var1 = hi_detrend;
                 var2 = hi_detrend;
                 if strcmp(var,'volume-area')
                    var2 = hi_detrend2;
                 end
    
                 % compute correlations
                 for mo = 1:12
                    for step = 0:12
                       if mo+step > 12
                          [r,p] = corrcoef(var1(res,mo,1:nyrs-1),var2(res,mo+step-12,2:nyrs),...
                                           'rows','pairwise');
                           R(mod_num-1,mo,step+1) = r(1,2);
                       elseif step == 12
                          [r,p] = corrcoef(var1(res,mo,1:nyrs-1),var2(res,mo,2:nyrs),...
                                           'rows','pairwise');
                           R(mod_num-1,mo,step+1) = r(1,2);                    
                       else
                          [r,p] = corrcoef(var1(res,mo,1:nyrs-1),var2(res,mo+step,1:nyrs-1),...
                                           'rows','pairwise');
                           R(mod_num-1,mo,step+1) = r(1,2);
                       end
                       Z = 0.5 .* log((1+r(1,2)) / (1-r(1,2)));
                       standd = 1.96 .*(1 ./ sqrt((nyrs-1) - 3));
                       if Z > 0
                          conf = Z - stand;
                       else
                          conf = Z + standd;
                       end
                       conf = tanh(conf);
                       if abs(r(1,2)-conf) < abs(r(1,2)) 
                         P(mod_num-1,mo,step+1) = 1;
                       end
                    end
                 end 
             end % model_num          

             figure (1)
             coordx=[0:12];
             coordy=[1:12];
             [x,y]=meshgrid(coordx,coordy);
             Rtmp = R;
             R = nanmean(R,1);
             P = nanmean(P,1);
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
             print(gcf,['modern_',ocean_name,'_',var,'_',thehemi,'.png'],'-dpng','-r200');
             close
         end % ocean
    end % var
end % region

