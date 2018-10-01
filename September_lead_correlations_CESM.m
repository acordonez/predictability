% September_lead_correlations_CESM
%
% Calculates Pearson's R correlations for Arctic 
% and Antarctic sea ice in DYN or SOM. Loops over
% hemispheres and models.
%
% The correlations can be plotted using the script
% 'Sept_mar_line_correlation_plots_CESM.m' to 
% produce figure 6.
%
% Ana Ordonez 10/2018

% loop over hemispheres
for region = 1:2

   var_list = {'total_area','volume','area2','area3','area4','area5','my','fy'};
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

   model_case = 'LE';
   res_max= 2;

   map_stride = [1 1];
   R(8,res_max,13)=0;
   P=zeros(size(R));
   thestart=[1 1];
   for res = 2:2 % 1 = DYN, 2 = SOM
       clear area hi_detrend hi_detrend2 tmp tot_area trnd 

       model_list = {'b.e11.B1850C5CN.','e.e11.E1850C5CN.'};
       model = model_list{res};
       if strcmp(model,'e.e11.E1850C5CN.')
          thesuff = {'.050001-059912','.060001-069912','.070001-079912','.080001-089912','.090001-100112'};
          first_year{res} = '0500';
          model_num='001';
       else
           thesuff = {'.170001-179912','.180001-189912','.190001-199912','.200001-209912','.210001-220012'};
           first_year{res} = '1700';
           model_num ='005';
       end
       %disp(first_year);
       map_stride = [map_stride, 1];
       nyrs =100;
       yr_factor = 5;
       yrs = nyrs*yr_factor;
       eq = [inf inf];
       dir = '/glade/p/cesm0005/CESM-CAM5-BGC-LE/ice/proc/tseries/monthly/';
       nc1 = ['f09_g16.',model_num,'.cice.h.'];
       for fnum = 1:yr_factor
           nc2 = ['_',thehemi,thesuff{fnum},'.nc'];
           fh{fnum} = [dir,'hi','/',model,nc1,'hi',nc2];
           fa{fnum} = [dir,'aice','/',model,nc1,'aice',nc2];
           fd{fnum} = [dir,'daidtd','/',model,nc1,'daidtd',nc2];
           fa5{fnum} = [dir,'aicen005','/',model,nc1,'aicen005',nc2];
           fa4{fnum} = [dir,'aicen004','/',model,nc1,'aicen004',nc2];
           fa3{fnum} = [dir,'aicen003','/',model,nc1,'aicen003',nc2];
           fa2{fnum} = [dir,'aicen002','/',model,nc1,'aicen002',nc2];
           ff{fnum} = [dir,'FYarea','/',model,nc1,'FYarea',nc2];
       end
       area = ncread(fh{1},'tarea',[1 1],[inf eq(res)],...
                    [map_stride(res) map_stride(res)]); 
       area = repmat(area,1,1,yrs);
       mask = repmat(mask,1,1,yrs);
       pacific = mask;
       pacific(mask ~= 3 & mask ~= 11 & mask ~= 12 & mask ~= 13) = NaN;
       pacific(~isnan(pacific)) = 1;
       atlantic = mask;
       atlantic(mask ~= 7 & mask ~= 8) = NaN;
       atlantic(~isnan(atlantic)) = 1;
       add_yr = 0*12;
              
       for mo=1:12
          % get monthly time series
              
          % total ice area and ice area categories 
          total_area = [];
          for fnum = 1:yr_factor
              tmp1= ncread(fa{fnum},'aice',[thestart mo+add_yr],[inf eq(res) nyrs],...
                           [map_stride(res) map_stride(res) 12]);
              total_area = cat(3,total_area,tmp1);
          end
          total_area(total_area==0) = NaN;
          area_frac_map = total_area;
          total_area = total_area .* area ./ 100; % ice area
          total_area = squeeze(nansum(nansum(total_area))); %total ice area

          area2 = []; area3 = []; area4 = []; area5 = [];
          for fnum = 1:yr_factor
              ftmp5= ncread(fa5{fnum},'aicen005',[thestart mo+add_yr],[inf eq(res) nyrs],...
                            [map_stride(res) map_stride(res) 12]);
              ftmp4= ncread(fa4{fnum},'aicen004',[thestart mo+add_yr],[inf eq(res) nyrs],...
                            [map_stride(res) map_stride(res) 12]);                       
              ftmp3= ncread(fa3{fnum},'aicen003',[thestart mo+add_yr],[inf eq(res) nyrs],...
                            [map_stride(res) map_stride(res) 12]);
              ftmp2= ncread(fa2{fnum},'aicen002',[thestart mo+add_yr],[inf eq(res) nyrs],...
                            [map_stride(res) map_stride(res) 12]);
              area2tmp = ftmp2 + ftmp3 + ftmp4 + ftmp5;
              area3tmp = ftmp3 + ftmp4 + ftmp5;
              area4tmp = ftmp4 + ftmp5;
              area5tmp = ftmp5;
              area2 = cat(3,area2,area2tmp);
              area3 = cat(3,area3,area3tmp);
              area4 = cat(3,area4,area4tmp);
              area5 = cat(3,area5,area5tmp);
          end
          area2(area2==0) = NaN;
          area2 = area2 .* area ./ 100; % ice area
          area2 = squeeze(nansum(nansum(area2)));                    
          area3(area3==0) = NaN;
          area3 = area3 .* area ./ 100; % ice area
          area3 = squeeze(nansum(nansum(area3))); 
          area4(area4==0) = NaN;
          area4 = area4 .* area ./ 100; % ice area
          area4 = squeeze(nansum(nansum(area4))); 
          area5(area5==0) = NaN;
          area5 = area5 .* area ./ 100; % ice area
          area5 = squeeze(nansum(nansum(area5))); 
         
          % first year and multiyear ice area
          fy = [];
          for fnum = 1:yr_factor
             tmp1= ncread(ff{fnum},'FYarea',[thestart mo+add_yr],[inf eq(res) nyrs],...
                            [map_stride(res) map_stride(res) 12]);
                 fy = cat(3,fy,tmp1);
          end
          fy(fy==0) = NaN;
          my = area_frac_map - (area_frac_map .* fy);
          fy = fy .* area ./ 100; % ice area
          my = my .* area ./ 100;
          my = squeeze(nansum(nansum(my))); %total multiyear ice area
          fy = squeeze(nansum(nansum(fy))); %total ice area

          % ice volume
          volume = [];
          for fnum = 1:yr_factor
                 tmp1= ncread(fh{fnum},'hi',[thestart mo+add_yr],[inf eq(res) nyrs],...
                              [map_stride(res) map_stride(res) 12]);
                 volume = cat(3,volume,tmp1);
          end
          volume(volume == 0) = NaN;
          volume = volume .* area; % volume/grid area .* grid area = volume
          volume = squeeze(nansum(nansum(volume,2),1)); % total volume
              
          for var_ind = 1:8
             vartmp = eval(var_list{var_ind});

             % detrend before correlating
             trnd = gettrend(squeeze(vartmp));
             hi_detrend(var_ind,res,mo,:) = vartmp'-[0:yrs-1].*trnd./length(vartmp);
          end       
       end % mo
    end % res
  
    for mo = 3:6:9
       for var_ind = 1:8
          % set variables to correlate
           var1 = squeeze(hi_detrend(1,res,:,:));
           var2 = squeeze(hi_detrend(var_ind,res,:,:));

           % compute correlations
           for step = -12:1:0
               if (mo+step) < 1
                  [r] = corrcoef(var1(mo,2:yrs),...
                                 var2(mo+step+12,1:yrs-1),...
                                 'rows','pairwise');
                  R(var_ind,res,step+13) = r(1,2);
               else
                  [r] = corrcoef(var1(mo,1:yrs-1),...
                                 var2(mo+step,1:yrs-1),...
                                 'rows','pairwise');
                  R(var_ind,res,step+13) = r(1,2);                   
               end
               Z = 0.5 .* log((1+r(1,2)) / (1-r(1,2)));
               stand = 1.96 .*(1 ./ sqrt((yrs-1) - 3));
               if Z > 0
                  conf = Z - stand;
               else
                  conf = Z + stand;
               end
               conf = tanh(conf);
               if abs(r(1,2)-conf) < abs(r(1,2)) 
                  P(var_ind,res,step+13) = 1;
               end
           end
       end
       save(['CESM_sept_',thehemi,'_',model,'_',num2str(mo),'.mat'],'R','P');
    end % mo
end %region
