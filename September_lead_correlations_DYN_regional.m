% September_lead_correlations_CESM.m
%
% Calculates Pearson's R correlations for Arctic 

% load sea ice outlook mask and regrid
ncfile = 'map_stereo25km_TO_gx1v6NH_blin.161123.nc';
mask1 = ncread('sio_2016_mask.nc','mask');
mask1=double(mask1);
mask1=mask1';
mask = grid1togrid2(mask1,ncfile);
mask = permute(mask,[2 1 3]);
mask = round(mask);
3,11,12,13 for Pacific
7,8 for Atlantic

mask_name = 'atlantic';
var_list = {'total_area','volume','area2','area3','area4','area5','my','fy','sst'};
thehemi = 'nh';

map_stride = [1 1];
R(9,1,13)=0;
P=zeros(size(R));
thestart=[1 1];

model_list = {'b.e11.B1850C5CN.','e.e11.E1850C5CN.'};
model = model_list{res};
if strcmp(model,'e.e11.E1850C5CN.')
          %thesuff = {'.090001-100112'};
          thesuff = {'.050001-059912','.060001-069912','.070001-079912','.080001-089912','.090001-100112'};
          first_year{res} = '0500';
          model_num='001';
else
           %thesuff = '.140001-149912';
           %first_year{res} = '1400'; 
           %thesuff = {'.210001-220012'};
           thesuff = {'.170001-179912','.180001-189912','.190001-199912','.200001-209912','.210001-220012'};
           first_year{res} = '1700';
           model_num ='005';
end
disp(first_year);
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
           fH{fnum} = ['/glade/p/cesm0005/CESM-CAM5-BGC-LE/ocn/proc/tseries/monthly/SST/',...
                    model,'f09_g16.',model_num,'.pop.h.SST',thesuff{fnum},'.nc'];
end
area = ncread(fh{1},'tarea',[1 1],[inf eq(res)],...
                     [map_stride(res) map_stride(res)]); 
area = repmat(area,1,1,yrs);
add_yr = 0*12;

if size(size(mask),2) < 3
          mask = repmat(mask,1,1,yrs);
          pacific = mask;
          pacific(mask ~= 3 & mask ~= 11 & mask ~= 12 & mask ~= 13) = NaN;
          pacific(~isnan(pacific)) = 1;
          atlantic = mask;
          atlantic(mask ~= 7 & mask ~= 8) = NaN;
          atlantic(~isnan(atlantic)) = 1;
end

if strcmp(mask_name,'atlantic')
           region_mask = atlantic;
           name1 = ['LE_sept_',thehemi,'_',model,'_atlantic_9_tmp.mat'];
           name2 = ['LE_sept_',thehemi,'_',model,'_atlantic_3_tmp.mat'];
else 
           region_mask = pacific;
           name1 = ['LE_sept_',thehemi,'_',model,'_pacific_9_tmp.mat'];
           name2 = ['LE_sept_',thehemi,'_',model,'_pacific_3_tmp.mat'];
end
              
disp('computing correlations')
for mo=1:12
    %get monthly time series
               
   total_area = [];
   for fnum = 1:yr_factor
      tmp1= ncread(fa{fnum},'aice',[thestart mo+add_yr],[inf eq(res) nyrs],...
                   [map_stride(res) map_stride(res) 12]);
      total_area = cat(3,total_area,tmp1);
   end
   total_area(total_area==0) = NaN;
   area_frac_map = total_area;
   total_area = region_mask .* total_area .* area ./ 100; % ice area
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
   area2 = region_mask .* area2 .* area ./ 100; % ice area
   area2 = squeeze(nansum(nansum(area2)));                    
   area3(area3==0) = NaN;
   area3 = region_mask .* area3 .* area ./ 100; % ice area
   area3 = squeeze(nansum(nansum(area3))); 
   area4(area4==0) = NaN;
   area4 = region_mask .* area4 .* area ./ 100; % ice area
   area4 = squeeze(nansum(nansum(area4))); 
   area5(area5==0) = NaN;
   area5 = region_mask .* area5 .* area ./ 100; % ice area
   area5 = squeeze(nansum(nansum(area5))); 
         
   fy = [];
   for fnum = 1:yr_factor
                 tmp1= ncread(ff{fnum},'FYarea',[thestart mo+add_yr],[inf eq(res) nyrs],...
                            [map_stride(res) map_stride(res) 12]);
                 fy = cat(3,fy,tmp1);
   end
   fy(fy==0) = NaN;
   my = area_frac_map - (area_frac_map .* fy);
   fy = region_mask .* fy .* area ./ 100; % ice area
   my = region_mask .* my .* area ./ 100;
   my = squeeze(nansum(nansum(my))); %total multiyear ice area
   fy = squeeze(nansum(nansum(fy))); %total ice area

   volume = [];
   for fnum = 1:yr_factor
                 tmp1= ncread(fh{fnum},'hi',[thestart mo+add_yr],[inf eq(res) nyrs],...
                              [map_stride(res) map_stride(res) 12]);
                 volume = cat(3,volume,tmp1);
   end
   volume(volume == 0) = NaN;
   volume = region_mask .* volume .* area; % volume/grid area .* grid area = volume
   volume = squeeze(nansum(nansum(volume,2),1)); % total volume

   if strcmp(model,'b.e11.B1850C5CN.')
                 sst = [];
                 for fnum = 1:yr_factor
                   tmp1= ncread(fH{fnum},'SST',[1 281 1 mo+add_yr],[inf eq(res) 1 nyrs],...
                                [map_stride(res) map_stride(res) 1 12]);
                   tmp1 = squeeze(tmp1);
                   sst = cat(3,sst,tmp1);
                 end
                 area_wgt = area .* region_mask ./ squeeze(nansum(nansum(area(:,:,1).*region_mask(:,:,1))));
                 sst = region_mask .* sst .* area_wgt;
                 sst = squeeze(nansum(nansum(sst))); %area weighted mean mixed layer depth
   else
                 sst = zeros(size(volume));
   end
              
   for var_ind = 1:9
                 vartmp = eval(var_list{var_ind});

                 % detrend before correlating
                 trnd = gettrend(squeeze(vartmp));
                 hi_detrend(var_ind,res,mo,:) = vartmp'-[0:yrs-1].*trnd./length(vartmp);
   end
              
end
        
  
for mo = 3:6:9
   for var_ind = 1:9
   %set variables to correlate
   var1 = squeeze(hi_detrend(1,res,:,:));
   var2 = squeeze(hi_detrend(var_ind,res,:,:));
   %compute correlations
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
   save(['Dyn_sept_',thehemi,'_',model,'_',mask_name,'_',num2str(mo),'.mat'],'R','P');
end



