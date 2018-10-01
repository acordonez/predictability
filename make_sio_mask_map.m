% make_sio_mask_map.m
% 
% Plots the Atlantic and Pacific regions for figure 7
%
% Ana Ordonez 10/2018

mask1 = ncread('sio_2016_mask.nc','mask');
masklat = ncread('sio_2016_mask.nc','lat');
masklon = ncread('sio_2016_mask.nc','lon');
mask1(mask1 == 1 | mask1 == 2 | mask1 == 4 | mask1 == 5 | mask1 == 6 | mask1 == 9 | mask1 == 10 | mask1 > 13) = 1;
mask1(mask1 == 7 | mask1 == 8)  = 2; %Atlantic
mask1(mask1 == 3 | mask1 == 11 | mask1 == 12 | mask1 == 13) = 3; %Pacific

load coast
figure (1)   
axesm('stereo','maplatlimit',[50 90], ...
   	'origin',[90 0],'Frame','on',...
        'Grid','on','GLineStyle',':','GLineWidth',0.5,...
        'GAltitude',0,...
        'MLineLocation',60,...
        'PLineLocation',[60],...
   	'MeridianLabel','off','ParallelLabel','off',...
   	'MLabelLocation',[180 120 60 0 -60 -120],...
   	'PLabelLocation',[60],... 
        'LabelRotation','off',...         
   	'fontsize',8,'fontname','arial'); axis off   
lat_bgrnd = 0:1:90;
lon_bgrnd = -180:1:180;
[LAT_bgrnd,LON_bgrnd] = meshgrid(lat_bgrnd,lon_bgrnd);
bgrnd = ones(size(LAT_bgrnd));
geoshow(LAT_bgrnd,LON_bgrnd,bgrnd,'DisplayType','texturemap');
geoshow(masklat,masklon,mask1,'Displaytype','texturemap');
geoshow(lat,long,'displaytype','polygon','FaceColor',[1 1 0.8],'EdgeColor',[0 0 0],'LineWidth',0.5);
map = [251, 154, 153; 51, 160, 44;  31, 120, 180; 227, 26, 28; 253, 191, 111; 166, 206, 227;178, 223, 138;] ./ 255; %from colorbrewer2.org
colormap(map)
print(gcf,'sio_map.jpeg','-djpeg','-r300'); 

