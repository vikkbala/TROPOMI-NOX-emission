function [lifetime_est, emission_est, al, mu, x0, s, b] = emg_fitting(no2_bin, inter_bin, bin_lon_min, bin_lon_max, era5_ws_ano)


% fitting the unknown parameters, and estimating NOX emission and lifetime
no2_bin_conv = no2_bin.* inter_bin*1000;
lin_den = movmean(nansum(no2_bin_conv,1), 2);

[max_ld, max_ld_ind] = max(lin_den);

xData = (bin_lon_min + bin_lon_max)./2;
yData = lin_den;

% lower and upper bound of 5 unknown fitted parameters
lb = [0, xData(1), inter_bin/3, min(xData), 0]; 
ub = [inf, xData(end), inf, xData(max_ld_ind), max_ld]; 

% Initial parameter guess
initialParams = [sum(lin_den,'all'), xData(max_ld_ind), era5_ws_ano*(4*60*60)/1000, xData(max_ld_ind)/4, min(lin_den)];

% Options for fmincon
options = optimoptions('fmincon', 'Algorithm', 'interior-point', 'Display', 'iter');

% Perform the optimization
params = fmincon(@(params) emgfittingObjective(params, xData, yData), initialParams, [], [], [], [], lb, ub, [], options);

% fitted parameters
al = params(1);
mu = params(2);  
x0 = params(3);
s = params(4);
b = params(5); 
x = xData;

fitted_curve = (al/2) * exp((s^2 ./ (2 * x0^2)) - ((x - mu)/x0)) .* erfc (((s^2) - (x0*(x-mu))) / (2^(1/2) * s *x0)) + b;

% Results
lifetime_est = x0*1000/era5_ws_ano/3600;
emission_est = 1.32 * (al * 46.01* era5_ws_ano);

figure     
plot (xData, yData, 'color','b','linewidth',3,'LineStyle','-','DisplayName', 'Observed')
hold on
plot (xData, fitted_curve, 'color','r','linewidth',3,'LineStyle','-','DisplayName', 'Fitted')
grid minor
legend()
ylim([0 max(lin_den)+1])
% text(0,1, strcat('NO_X Lifetime (hr):',num2str(lifetime_est)),'fontweight','bold','FontSize', 14,'FontName', 'Times New Roman')
% text(0,2, strcat('NO_X Emission (g/s):',num2str(emission_est)),'fontweight','bold','FontSize', 14,'FontName', 'Times New Roman')
% text(0,3, strcat('Wind speed (m/s):',num2str(era5_ws_ano)),'fontweight','bold','FontSize', 14,'FontName', 'Times New Roman')
xlabel('distance (km)')
ylabel('Line density (mole/m)')
set(gca, 'YDir', 'normal','FontSize', 14,'fontweight','bold','FontName', 'Times New Roman')

end