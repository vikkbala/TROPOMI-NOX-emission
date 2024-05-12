function [x_len, y_len, bin_lon_min, bin_lon_max, bin_lat_min_mat, bin_lat_max_mat, bin_lon_min_mat, bin_lon_max_mat, bin_lon_mat, bin_lat_mat, inves_tropomi_data] = data_extract (myFolder_tropomi, inves_date, lat_sou, lon_sou, inter_bin, tot_ran_y_abo, tot_ran_y_bel, tot_ran_x_upwind, tot_ran_x_dowwind, sour_sou) 

% Data manipulation (i.e, TROPOMI NO2 data extraction)

filePattern = fullfile(myFolder_tropomi, '*.nc4');
ncFiles   = dir(filePattern);

lon_min = lon_sou-sour_sou; lon_max = lon_sou+sour_sou;    
lat_min = lat_sou-sour_sou; lat_max = lat_sou+sour_sou;

bin_lat_min = flip(-tot_ran_y_bel:inter_bin:tot_ran_y_abo-inter_bin)';
bin_lat_max = flip ((-tot_ran_y_bel+inter_bin):inter_bin:tot_ran_y_abo)';

bin_lon_min = -tot_ran_x_upwind:inter_bin:(tot_ran_x_dowwind-inter_bin);
bin_lon_max = (-tot_ran_x_upwind+inter_bin):inter_bin:tot_ran_x_dowwind;

x_len = length(bin_lon_min); y_len = length(bin_lat_min);

bin_lat_min_mat = repmat(bin_lat_min,1,x_len);
bin_lat_max_mat = repmat(bin_lat_max,1,x_len);
bin_lon_min_mat = repmat(bin_lon_min,y_len,1);
bin_lon_max_mat = repmat(bin_lon_max,y_len,1);

bin_lon_mat = (bin_lon_min_mat + bin_lon_max_mat)/2;
bin_lat_mat = (bin_lat_min_mat + bin_lat_max_mat)/2;

file_idx_inves = [];
for ncFiles_xx = 1:length(ncFiles)
    baseFileName = ncFiles(ncFiles_xx).name; 
    date_FileName = str2double(baseFileName(21:28));

    if date_FileName == inves_date
        file_idx_inves = [file_idx_inves; ncFiles_xx];
    end
end

tropomi_struct = struct();

if ~isempty (file_idx_inves)
    
    for file_idx_inves_xx = 1:length(file_idx_inves)
        baseFileName_inves = ncFiles(file_idx_inves(file_idx_inves_xx)).name;
        file = fullfile(myFolder_tropomi, baseFileName_inves);
        lons = ncread(file,'/PRODUCT/longitude');  lats = ncread(file,'/PRODUCT/latitude');
        lons_bou = ncread(file,'/PRODUCT/SUPPORT_DATA/GEOLOCATIONS/longitude_bounds');  
        lats_bou = ncread(file,'/PRODUCT/SUPPORT_DATA/GEOLOCATIONS/latitude_bounds');
        
        no2_data = ncread(file,'/PRODUCT/nitrogendioxide_tropospheric_column');
        qa_data  = ncread(file,'/PRODUCT/qa_value');
        
        del_time = ncread(file,'/PRODUCT/delta_time');
        start_time = ncread(file,'/PRODUCT/time');
        time_mani = datetime(2010, 1, 1, 0, 0, 0, 'TimeZone','UTC') - datetime(1970, 1, 1, 0, 0, 0, 'TimeZone','UTC');
        start_time_cha = (start_time + seconds(time_mani)); % epoch time (seconds since 1970)
        meas_time = del_time/1000 + start_time_cha; meas_time = repmat(meas_time',size(lons,1),1);
        
        lat_bou1 = squeeze(lats_bou(1,:,:)); lat_bou2 = squeeze(lats_bou(2,:,:));
        lat_bou3 = squeeze(lats_bou(3,:,:)); lat_bou4 = squeeze(lats_bou(4,:,:));
        
        lon_bou1 = squeeze(lons_bou(1,:,:)); lon_bou2 = squeeze(lons_bou(2,:,:));
        lon_bou3 = squeeze(lons_bou(3,:,:)); lon_bou4 = squeeze(lons_bou(4,:,:));
        
        qa_do = find(lats >= lat_min & lats <= lat_max & lons >= lon_min & lons <= lon_max & qa_data >= 0.75 & no2_data > 0);
        
        
        if length(qa_do) >= 500
            lat_qa = lats(qa_do); lon_qa = lons(qa_do); 
            no2_qa = no2_data(qa_do); meas_time_qa = meas_time(qa_do);
            lat_bou1_qa = lat_bou1(qa_do); lat_bou2_qa = lat_bou2(qa_do); 
            lat_bou3_qa = lat_bou3(qa_do); lat_bou4_qa = lat_bou4(qa_do); 
            lon_bou1_qa = lon_bou1(qa_do); lon_bou2_qa = lon_bou2(qa_do); 
            lon_bou3_qa = lon_bou3(qa_do); lon_bou4_qa = lon_bou4(qa_do); 
            
            fie_name = genvarname(baseFileName_inves(21:35));
            
            tropomi_struct.(fie_name) = struct();
            tropomi_struct.(fie_name).lat = lat_qa;
            tropomi_struct.(fie_name).lon = lon_qa;
            tropomi_struct.(fie_name).no2 = no2_qa;
            tropomi_struct.(fie_name).meas_time = meas_time_qa;
            
            tropomi_struct.(fie_name).lat_bou1 = lat_bou1_qa;
            tropomi_struct.(fie_name).lat_bou2 = lat_bou2_qa;
            tropomi_struct.(fie_name).lat_bou3 = lat_bou3_qa;
            tropomi_struct.(fie_name).lat_bou4 = lat_bou4_qa;
            
            tropomi_struct.(fie_name).lon_bou1 = lon_bou1_qa;
            tropomi_struct.(fie_name).lon_bou2 = lon_bou2_qa;
            tropomi_struct.(fie_name).lon_bou3 = lon_bou3_qa;
            tropomi_struct.(fie_name).lon_bou4 = lon_bou4_qa;
        end

    end
end

if numel(fieldnames(tropomi_struct)) >=1
    fie_names = fieldnames(tropomi_struct);
    if length(fie_names) >= 2
        disp('There is more than one TROPOMI data available for this investigation date. Below are starting of filenames')
        disp(fie_names)
        user_option = input('Choose which file to investiate (1 or 2 or 3 or etc..): ');
    else
        user_option = 1;
    end
    inves_tropomi_data = tropomi_struct.(char(fie_names(user_option)));
else
    inves_tropomi_data = [];

end

end

