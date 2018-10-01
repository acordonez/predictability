function [ice,lat,lon] = get_monthly_sie_from_obs(thehemi)
% [ice, lat, lon] = get_monthly_sie_from_obs(thehemi)
%
% Reads ice concentration from satellite data
% and transforms into ice extent. Returns spatial
% area info on stereo polar grid. 
%
% Parameters:
% ice(longitude,latitude,month)
% thehemi: 'n' for Northern Hemisphere, 's' for Southern
%
% The satellite data is stored in a folder called 'sat_ice'. 
% More information on the data used can be found in our paper.
%
% Ana Ordonez 10/2018

if strcmp(thehemi,'n') 
    zone = [304,448];
else
    zone = [316,332];
end

% get coordinate info
datdir = '..';
lat_id = fopen([datdir,'sat_ice/ps',thehemi,'25lats_v3.dat'],'r','l');

lat = fread(lat_id,zone,'int')./ 100000;
fclose(lat_id);
lon_id = fopen([datdir,'sat_ice/ps',thehemi,'25lons_v3.dat']);
lon = fread(lon_id,zone,'int') ./ 100000;
area_id = fopen([datdir,'sat_ice/ps',thehemi,'25area_v3.dat'],'r','l');
area = fread(area_id,zone,'int') ./ 1000; %sq km
fclose(area_id);

% read all the files here
ice = [];

for num = 1:468
   %read data for 1979-2017
   fname = dir([datdir,'sat_ice/',thehemi,'/v3/bt*']);
   fname = {fname.name};
   if ((num == 108) | (num == 109 ))
       % Dec 1987 and Jan 1988 are missing in version 3
       tmp = NaN(zone);
   else 
       fnum = num;
       if num > 109
           % correcting indices for missing data
           fnum = num-2;
       end
       f1 = [datdir,'sat_ice/',thehemi,'/v3/',fname{fnum}];
       id = fopen(f1,'r','l');  
       tmp = fread(id,zone,'int16')./10;
       fclose(id);
   end
   ice = cat(3,ice,tmp);
end

% calculate extent
ice(ice < 15) = 0;
ice(ice == 120) = 0;
ice(ice >= 15) = 1;
nmnths = size(ice,3);
ice = ice .* repmat(area,1,1,nmnths);

end
