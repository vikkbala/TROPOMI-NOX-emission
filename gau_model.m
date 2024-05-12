function [gau_model_bin, along_wind, no2_vc_summed] =  gau_model(prior_emi, lifetime_gau, lat_sou, lon_sou, era5_ws_ano, era5_wd_ano)

gau_model_bin = 1;

ws_sp = [2,5];
sp = [213, 104];

if era5_ws_ano <= 2
    sta_parm = 213;
elseif era5_ws_ano >= 5
    sta_parm = 104;
else
    sta_parm = interp1(ws_sp, sp, era5_ws_ano);
end
     
% along_wind = -111.195*2:11.1195/10:111.195*2; % in km
% across_wind = (-111.195*2:11.1195/10:111.195*2)*1000; % in meters

along_wind = -100*2:gau_model_bin:100*2; % in km
across_wind = (-100*2:gau_model_bin:100*2)*1000; % in meters

no2_vc = nan(length(across_wind),length(along_wind),length(prior_emi));

for rep_emi_xx = 1:length(prior_emi)
    for along_wind_xx = round(length(along_wind)/2):length(along_wind)
        for across_wind_xx = 1:length(across_wind)
            gaus_mod = (prior_emi(rep_emi_xx) / ((2*3.14)^(1/2) * sta_parm * along_wind(along_wind_xx)^(0.894) * era5_ws_ano)) * exp ((-1/2)*(across_wind(across_wind_xx) / (sta_parm * along_wind(along_wind_xx)^(0.894)))^2);
            no2_vc(across_wind_xx,along_wind_xx,rep_emi_xx) = gaus_mod * exp(-(along_wind(along_wind_xx)*1000)./ (lifetime_gau*3600 * era5_ws_ano));
        end
    end
end

lon_sou_main = (lon_sou(1) * cosd(era5_wd_ano))+(lat_sou(1) * sind(era5_wd_ano));
lat_sou_main = (-lon_sou(1) * sind(era5_wd_ano))+(lat_sou(1) * cosd(era5_wd_ano));

no2_vc_summed = no2_vc(:,:,1);

if length(prior_emi) > 1
    
    for rep_emi_xxx = 2:length(prior_emi)
        lon_sou_loop = (lon_sou(rep_emi_xxx) * cosd(era5_wd_ano))+(lat_sou(rep_emi_xxx) * sind(era5_wd_ano));
        lat_sou_loop = (-lon_sou(rep_emi_xxx) * sind(era5_wd_ano))+(lat_sou(rep_emi_xxx) * cosd(era5_wd_ano));
        
        sou_loop_dist_y = deg2km ((lat_sou_loop - lat_sou_main),'earth');
        sou_loop_dist_x = deg2km ((lon_sou_loop - lon_sou_main),'earth');
        
        x_add = round(sou_loop_dist_x/(gau_model_bin));
        y_add = round(sou_loop_dist_y/(gau_model_bin));
        
        if x_add <= 0
            no2_vc_summed = horzcat(zeros(length(across_wind),abs(x_add)),no2_vc_summed);
            no2_vc_loop = horzcat(no2_vc(:,:,rep_emi_xxx),zeros(length(across_wind),abs(x_add)));
        elseif x_add > 0
            no2_vc_summed = horzcat(no2_vc_summed,zeros(length(across_wind),abs(x_add))); 
            no2_vc_loop = horzcat(zeros(length(across_wind),abs(x_add)),no2_vc(:,:,rep_emi_xxx));
        end
        
        if y_add < 0 
            no2_vc_summed = vertcat(zeros(abs(y_add),length(no2_vc_summed)),no2_vc_summed);
            no2_vc_loop = vertcat(no2_vc_loop,zeros(abs(y_add),length(no2_vc_loop)));
        elseif y_add > 0
            no2_vc_summed = vertcat(no2_vc_summed,zeros(abs(y_add),length(no2_vc_summed)));
            no2_vc_loop = vertcat(zeros(abs(y_add),length(no2_vc_loop)),no2_vc_loop);
        end
        
        no2_vc_summed = nansum(cat(3, no2_vc_summed, no2_vc_loop),3);
        [R, C] = size(no2_vc_summed);
        
        if x_add < 0
            no2_vc_summed (:,1:abs(x_add)) = [];
        elseif x_add > 0
            no2_vc_summed (:,C-abs(x_add)+1:C) = [];
        end
        
        if y_add < 0
            no2_vc_summed (1:abs(y_add),:) = [];
        elseif y_add > 0
            no2_vc_summed (R-abs(y_add)+1:R,:) = [];
        end
    end
end


end