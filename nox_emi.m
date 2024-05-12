clear; clc
%% Power plant locations 

% lat_sou = 21.4125018; lon_sou = 79.9673098;
% lat_sou = 14.697282; lon_sou = 78.459364; % Rayalasemma 
% lat_sou = 16.499546; lon_sou = 75.834632; % Kudgi Super Thermal Power Project
%% Required fields

myFolder_tropomi = 'D:\TROPOMI_NO2_india\downloads'; % TROPOMI NO2 folder
myFolder_era5 = 'D:\ERA5';     

inves_date = 20221113; % investigation date of power plant (in YYYYMMDD)

% lat_sou = [18.7550]; lon_sou = [79.4561]; 18.827988546820343, 79.5751079532966
lat_sou = [18.7550, 18.8270561]; lon_sou = [79.4561, 79.5702839]; % source location i.e., location of power plant
% lat_sou = [18.75910671989932, 18.827988546820343]; lon_sou = [79.45518625981245, 79.5751079532966];
% lat_sou = 16.499546-0.1; lon_sou = 75.834632-0.1;

prior_emi = [1000, 1000]; % in g/s (if it is unknown, better give 100) 

inter_bin = 1;            % the grid size (1 = 1km by 1km)
tot_ran_y_abo = 40;       % distance (km) above the emission source considered in emission calculation, i.e., tot_ran_y = 20 --> then 20km above the emission source)
tot_ran_y_bel = 40;       % distance (km) below the emission source considered in emission calculation.
tot_ran_x_upwind = 50;    % distance (km) in upwind direction of the emission source considered in emission calculation.
tot_ran_x_dowwind = 100;  % distance (km) in downwind direction of the emission source considered in emission calculation.

pre_lev = 1000;  % pressure level of ERA5 wind information used in this analysis and emission estimation (it could be multiple level, e.g., [1000, 975, 950], then mean value will be considered)

sour_sou = 1;    % rectangular domain around the source considered in this analysis (1 = 1 degree): Don't required to be change

allow_win_rot = 0; % angle to rotate plume in addition to wind info from ERA-5 (negative value makes plume rotate clockwise and vice versa for positive value)

EMG_method = 'no'; GAU_method = 'yes'; CSF_method = 'no'; 
%% Emission estimation: Exponentially modified gaussian method and/or Cross-sectional emission flux method

emi_method (myFolder_tropomi, myFolder_era5, inves_date, lat_sou, lon_sou, inter_bin, tot_ran_y_abo, tot_ran_y_bel, tot_ran_x_upwind, tot_ran_x_dowwind, pre_lev, sour_sou, allow_win_rot, prior_emi, EMG_method, GAU_method, CSF_method);


