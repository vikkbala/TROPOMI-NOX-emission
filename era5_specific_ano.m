
function [era5_ws_ano, era5_wd_ano, era5_wd_ano_org] =  era5_specific_ano(myFolder_era5, year_file, epoch_anom, pre_lev, lat_sou, lon_sou)

era5_file_name = ['ind_',year_file,'_7_9_900_1000.nc'];

era5_file = fullfile(myFolder_era5, era5_file_name);
era5_u = ncread(era5_file,'u'); era5_v = ncread(era5_file,'v'); 
era5_lat = ncread(era5_file,'latitude'); era5_lon = ncread(era5_file,'longitude'); 
era5_time = ncread(era5_file,'time'); era5_level = ncread(era5_file,'level');

time_mani = datetime(1970, 1, 1, 0, 0, 0, 'TimeZone','UTC') - datetime(1900, 1, 1, 0, 0, 0, 'TimeZone','UTC');
era5_time_cha = (era5_time - hours(time_mani))*60*60; % epoch time (seconds since 1970)
    
[~, era5_lat_near_ano] = (min(abs(lat_sou-era5_lat))); 
[~, era5_lon_near_ano] = (min(abs(lon_sou-era5_lon)));

[~, era5_time_near_ano] = (min(abs(epoch_anom-era5_time_cha)));

if abs(min(abs(epoch_anom-era5_time_cha))) <= 10800
    era5_u_ano = (era5_u(era5_lon_near_ano,era5_lat_near_ano,find(ismember(era5_level, pre_lev)),era5_time_near_ano-2:era5_time_near_ano+2)); 
    era5_v_ano = (era5_v(era5_lon_near_ano,era5_lat_near_ano,find(ismember(era5_level, pre_lev)),era5_time_near_ano-2:era5_time_near_ano+2));
    
    era5_u_ano_inter = interp1(double(era5_time_cha(era5_time_near_ano-2:era5_time_near_ano+2)), squeeze(era5_u_ano), epoch_anom, 'linear');
    era5_v_ano_inter = interp1(double(era5_time_cha(era5_time_near_ano-2:era5_time_near_ano+2)), squeeze(era5_v_ano), epoch_anom, 'linear');

    era5_u_ano_aprox = (era5_u(era5_lon_near_ano,era5_lat_near_ano,find(ismember(era5_level, pre_lev)),era5_time_near_ano)); 
    era5_v_ano_aprox = (era5_v(era5_lon_near_ano,era5_lat_near_ano,find(ismember(era5_level, pre_lev)),era5_time_near_ano));

    era5_ws_ano = ((era5_u_ano_inter.^2) + (era5_v_ano_inter.^2)).^(1/2); 
    era5_wd_ano_org = 180 + ((180/ pi) * atan2 (era5_u_ano_aprox, era5_v_ano_aprox));

    if era5_wd_ano_org >= 0 && era5_wd_ano_org <= 90
        era5_wd_ano = 270 - era5_wd_ano_org ;
    elseif era5_wd_ano_org >= 90 && era5_wd_ano_org <= 180
        era5_wd_ano = 180 - era5_wd_ano_org + 90;
    elseif era5_wd_ano_org >= 180 && era5_wd_ano_org <= 270
    era5_wd_ano = 270 - era5_wd_ano_org;
    elseif era5_wd_ano_org >= 270 && era5_wd_ano_org <= 360
        era5_wd_ano = 360 - era5_wd_ano_org + 270;
    end

else 
    era5_ws_ano = [];
    era5_wd_ano = [];
    era5_wd_ano_org = [];
end

end

% year_file = '2022';
% epoch_anom = nanmean(meas_time_qa);