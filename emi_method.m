function emi_method (myFolder_tropomi, myFolder_era5, inves_date, lat_sou, lon_sou, inter_bin, tot_ran_y_abo, tot_ran_y_bel, tot_ran_x_upwind, tot_ran_x_dowwind, pre_lev, sour_sou, allow_win_rot, prior_emi, EMG_method, GAU_method, CSF_method)

[x_len, y_len, bin_lon_min, bin_lon_max, bin_lat_min_mat, bin_lat_max_mat, bin_lon_min_mat, bin_lon_max_mat, bin_lon_mat, bin_lat_mat, inves_tropomi_data]  = data_extract (myFolder_tropomi, inves_date, lat_sou(1), lon_sou(1), inter_bin, tot_ran_y_abo, tot_ran_y_bel, tot_ran_x_upwind, tot_ran_x_dowwind, sour_sou); 

if ~isempty(inves_tropomi_data)

    meas_time_inves_mean = nanmean(inves_tropomi_data.meas_time);
    lat_inves = inves_tropomi_data.lat; lon_inves = inves_tropomi_data.lon; 
    no2_inves = inves_tropomi_data.no2;

    lat_bou1_inves = inves_tropomi_data.lat_bou1; lat_bou2_inves = inves_tropomi_data.lat_bou2;
    lat_bou3_inves = inves_tropomi_data.lat_bou3; lat_bou4_inves = inves_tropomi_data.lat_bou4;

    lon_bou1_inves = inves_tropomi_data.lon_bou1; lon_bou2_inves = inves_tropomi_data.lon_bou2;
    lon_bou3_inves = inves_tropomi_data.lon_bou3; lon_bou4_inves = inves_tropomi_data.lon_bou4;

    inves_date_str = num2str(inves_date); inves_year = inves_date_str(1:4);

    [era5_ws_ano, era5_wd_ano, ~] = era5_specific_ano(myFolder_era5, inves_year, meas_time_inves_mean, pre_lev, lat_sou(1), lon_sou(1));
    a_r = 180+era5_wd_ano-(allow_win_rot); 

    if ~isempty (era5_ws_ano)
        no2_bin = data_grid(lat_sou(1), lon_sou(1), inter_bin, tot_ran_y_abo, tot_ran_x_upwind, no2_inves, lat_inves, lon_inves, lat_bou1_inves, lat_bou2_inves, lat_bou3_inves, lat_bou4_inves, lon_bou1_inves, lon_bou2_inves, lon_bou3_inves, lon_bou4_inves, x_len, y_len, a_r, bin_lat_min_mat, bin_lat_max_mat, bin_lon_min_mat, bin_lon_max_mat, bin_lon_mat, bin_lat_mat);

        if contains(EMG_method ,'yes')
            [lifetime_est_emg, emission_est_emg, al, mu, x0, s, b] = emg_fitting(no2_bin, inter_bin, bin_lon_min, bin_lon_max, era5_ws_ano);
            disp('EMG method is performed: Results are below')
            disp(['Estimated lifetime (hr): ', num2str(lifetime_est_emg)])
            disp(['Estimated NO_X emission (g/s): ', num2str(emission_est_emg)])
            disp(['Estimated NO_X emission (Kt/year): ', num2str(emission_est_emg*10^(-9)*365*24*60*60)])
            disp(['Wind speed (m/s): ', num2str(era5_ws_ano)])
            disp(['Five fitted parameters: ', [' al: ', num2str(al)], [' mu: ', num2str(mu)], [' x0: ', num2str(x0)], [' s: ', num2str(s)], [' b: ', num2str(b)]])
        end

        if contains(GAU_method ,'yes')
            [lifetime_est_gau, emission_est_gau] = gau_fitting(no2_bin, lat_sou, lon_sou, inter_bin, bin_lon_min, bin_lon_max, era5_ws_ano, prior_emi, a_r, tot_ran_x_dowwind, tot_ran_x_upwind);
            
            disp('Gau method is performed: Results are below')
            disp(['Estimated lifetime (hr): ', num2str(lifetime_est_gau)])
            disp(['Estimated NO_X emission (g/s): ', num2str(emission_est_gau')])
        end

        if contains(CSF_method ,'yes')

            bkg_csf = mean(no2_bin(:,tot_ran_x_upwind/inter_bin-3:tot_ran_x_upwind/inter_bin-1), 'all');
            no2_bin_enh = (no2_bin - bkg_csf) * 46.01 * inter_bin*1000;

            csf_emi_downwind = nan(1,5);

            if ~exist('lifetime_est', 'var')
                disp('Lifetime is not available from other methods. So user should provide the Lifetime for CSF method.');
                lifetime_est_emg = input('NOX lifetime (hr): ');
            end

            for csf_dis_xx = 1:5
                csf_emi = nansum(no2_bin_enh(:,(tot_ran_x_upwind/inter_bin)+csf_dis_xx), 'all') * era5_ws_ano;
                csf_emi_downwind(csf_dis_xx) = csf_emi/ exp(-(csf_dis_xx*inter_bin*1000)/(era5_ws_ano*lifetime_est_emg*3600));
            end

            disp('CSF method is performed: Results are below')
            disp(['Estimated NO_X emission (g/s): ', num2str(csf_emi_downwind)])
            disp(['Estimated mean NO_X emission (g/s): ', num2str(nanmean(csf_emi_downwind, 'all'))])

        end

    elseif isempty (era5_ws_ano)
        disp('There is no ERA5 wind data available for this power plant/investigation date. Therefore, the program stopped')
    end


elseif isempty(inves_tropomi_data)
    disp('There is no TROPOMI data available for this power plant/investigation date. Therefore, the program stopped')
end

end