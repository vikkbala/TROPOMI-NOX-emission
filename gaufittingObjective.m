function error = gaufittingObjective(params, Y, lat_sou, lon_sou, era5_ws_ano, a_r, tot_ran_x_dowwind, tot_ran_x_upwind, inter_bin, prior_emi)

    emission_gau = [];
    for parm_len = 1:length(prior_emi)
        emission_gau = [emission_gau; params(parm_len)];
    end
    
    % emission_gau = params(1);
    lifetime_gau = params(length(prior_emi)+1);

    [gau_model_bin, along_wind, no2_vc_summed] =  gau_model(emission_gau, lifetime_gau, lat_sou, lon_sou, era5_ws_ano, a_r);

    no2_ld = nansum(no2_vc_summed,1);

    % along_wind_req = along_wind(202:find(along_wind == tot_ran_x_dowwind));

    no2_ld_req = no2_ld(((length(along_wind)-1)/2)+2:find(along_wind == tot_ran_x_dowwind));
    no2_ld_req_mean =  nanmean(reshape(no2_ld_req, [inter_bin*gau_model_bin, tot_ran_x_dowwind/inter_bin]),1);

    model = zeros(1,(tot_ran_x_upwind+tot_ran_x_dowwind)/inter_bin);

    model(1:tot_ran_x_upwind/inter_bin) = Y(1:tot_ran_x_upwind/inter_bin);
    model(tot_ran_x_upwind/inter_bin+1:(tot_ran_x_upwind/inter_bin)+tot_ran_x_dowwind/inter_bin) = no2_ld_req_mean;

    error = sum((model - Y).^2);

end