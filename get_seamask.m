function [seamask,ocean_list] = get_seamask(thehemi,lat,lon)
% [seamask,ocean_list] = get_seamask(thehemi,lat,lon)
% returns a grid with numbered regions
%
% Parameters:
% thehemi: 'nh' or 'sh'
% lat: 2d array of latitude
% lon: 2d array of longitude
% seamask: the region array
% ocean_list: list of names of each region

% Ana Ordonez 10/2018

   lon = wrapTo180(lon);
   seamask = zeros(size(lon));
   if strcmp(thehemi,'sh')
      seamask(lon >= -60 & lon < 20) = 1;
      seamask(lon >= 20 & lon < 90) = 2;
      seamask(lon >= 90 & lon < 160) = 3;
      seamask(lon >= 160) = 4;
      seamask(lon < -130) = 4;
      seamask(lon >= -130 & lon < -60) = 5;
      ocean_list = {'Weddell','Indian','WPacific','Ross','ABSeas'};
               
   elseif strcmp(thehemi,'nh')
      disp('seasmask called')
      seamask(lon >= -45 & lon < 20 & lat < 80) = 5;% Greenland
      seamask(lon >= 20 & lon < 100 & lat < 80) = 6;%Kara & barents
      seamask(lon >= 90 & lat < 65 | lon < -135 & lat < 65) = 1;%Bering sea and Seas of Okhotsk and japan
      seamask(lon < -65 & lon >= -100 & lat < 70) = 3;% Hudson Bay
      seamask(seamask ~= 3 & lat < 80 & lat > 40 & lon < -45 & lon > -90) = 4;% Baffin, Labrador, St. Lawrence
      seamask(seamask ==0 & lat > 65 & lon >90) = 2;% Arctic
      seamask(seamask == 0 & lat > 65 & lon < -90) = 2;% also arctic
      seamask(seamask == 0 & lat > 70 & lon <= 90 & lon >= -90) = 2;% Arctic 
      ocean_list = {'N Pacific','Central Arctic','Barents & Kara',...
                 'Hudson Bay','NW Atlantic','GIN Seas'};
    end
end
