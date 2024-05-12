function [xdis, ydis] = diskm (a_r, lat_sou, lon_sou, lat_req, lon_req)

lon_sou_rot = (lon_sou*cosd(a_r))+(lat_sou*sind(a_r));
lat_sou_rot = (-lon_sou*sind(a_r))+(lat_sou*cosd(a_r));

lon_req_rot = (lon_req*cosd(a_r))+(lat_req*sind(a_r));
lat_req_rot = (-lon_req*sind(a_r))+(lat_req*cosd(a_r));

xdis = deg2km(lon_sou_rot - lon_req_rot);
ydis = deg2km(lat_sou_rot - lat_req_rot);

end