%-------------------------------------------------------------------------%
%                                                                         %                             
%                    SPACECRAFT GUIDANCE AND NAVIGATION                   %
%                              A.Y. 2021/2022                             %
%                                                                         %
%                        GIAN MARCO PALDINO - 968731                      %
%                                                                         %
%                               ASSIGNMENT 2                              %
%                                                                         %
%-------------------------------------------------------------------------%
% DISCLAIMERS:                                                            %
%          - This code was developed and tested using a Mac/Intel OSX     %
%            64bit machine.                                               %
%          - The file 'map.jpg' must be added to the path in order to     %
%            correctly visualize Earth texture in the plots.              %
%                                                                         %
%-------------------------------------------------------------------------%
%                                                                         %
%              RUN THIS SECTION BEFORE RUNNING THE EXERCISES              %
%                                                                         %
% Add all the files and folders in the current working directory to the   %
% MATLAB path (kernels, sgp4, tdm, tle)                                   %
addpath(genpath(pwd));                                                    %
%                                                                         %
% Load SPICE kernels:                                                     %
cspice_furnsh('assignment02.tm');                                         %
% Load ExoMars kernel  (CHANGE THE FILE PATH TO RUN THIS ON WINDOWS)      %          
cspice_furnsh('kernels/exomars.bsp'); % change / to \ for Windows         % 
%                                                                         %                                                                                
%-------------------------------------------------------------------------%

%% EXERCISE 1

% Clear memory workspace
clearvars; close all; clc

disp(' <strong> Exercise 1 </strong>')
disp(' ')

% Count total kernels number
fprintf('Total kernels number: %d\n', cspice_ktotal('ALL'));

%-------------------------------------------------------------------------%

% REQUEST 1

% Set initial and final epochs
t0 = '2021-Nov-19 14:25:39.652 UTC';
tf = '2021-Nov-23 14:00:00.000 UTC';
t0_et = cspice_str2et(t0);
tf_et = cspice_str2et(tf);

% Define time grid
time_grid = (t0_et:60:tf_et);

% Define initial mean state: position [km] and velocity [km/s]
rr0_mean = [23721.610, 17903.673, -49.918]';
vv0_mean = [-1.150987, 1.529718, 3.122389]';
x0_mean = [rr0_mean;vv0_mean];

% Get Earth's standard gravitational parameter from SPICE [km^3/s^2]
mu = cspice_bodvrd('EARTH','GM',1);

% Compute orbital elements
[a,e,inc,OM,om,th] = car2kep(rr0_mean,vv0_mean,mu); 

% Compute orbital period
T = 2*pi*sqrt(a^3/mu);

% Propagate the orbit for 1 orbital period 
tspan1 = (0:60:T);
[timevec1,X_tf1] = twobody_propagator(rr0_mean,vv0_mean,mu,tspan1,'Unperturbed');

% Plot orbit
figure()
plot3(rr0_mean(1),rr0_mean(2),rr0_mean(3),'o','MarkerSize',5,'linewidth',4)
hold on
plot3(X_tf1(:,1),X_tf1(:,2),X_tf1(:,3),'linewidth',1);
grid on 
box on
hold on
PlotEarth
xlabel('X [km]')
ylabel('Y [km]')
zlabel('Z [km]')
legend('Initial mean state')

% Propagate the orbit from initial epoch to final epoch
[timevec2,rv_eci] = twobody_propagator(rr0_mean,vv0_mean,mu,time_grid,'Unperturbed');

% Compute azimuth [rad], elevation [rad], range [km] and range rate [km/s]
% of the satellite for each station
[sat_azimuth.MI, sat_elevation.MI, sat_range.MI, sat_range_rate.MI] = ...
    antenna_pointing('MILANO', time_grid, rv_eci');

[sat_azimuth.WE, sat_elevation.WE, sat_range.WE, sat_range_rate.WE] = ...
    antenna_pointing('WELLINGTON', time_grid, rv_eci');

[sat_azimuth.LS, sat_elevation.LS, sat_range.LS, sat_range_rate.LS] = ...
    antenna_pointing('LA-SILLA', time_grid, rv_eci');

% Define minimum elevation treshold for each station [deg]
min_elevation = struct('MI', 30, 'WE',20, 'LS',10); 

% Compute visibility windows 
[t0_win_et, tf_win_et,vis_wind,vis_wind_et,vis_wind_date,ind_t0_win,ind_tf_win] = ...
 visibility_windows(sat_elevation,min_elevation,time_grid);

% Plot visibility windows
% Milan
figure()
subplot(1,3,1)
plot(time_grid/cspice_spd(), sat_elevation.MI*cspice_dpr(),'linewidth',1)
yl = yline(min_elevation.MI,'-','Min El','linewidth',1);
yl.LabelHorizontalAlignment = 'left';
yl.Color = [.85 .33 .10];
xlabel('Epoch [MJD2000]')
ylabel('Elevation [deg]')
xlim([time_grid(1)/cspice_spd(),time_grid(end)/cspice_spd()])
ylim([-90,90])
yticks(-90:30:90)
title('MILAN')
grid on

% Wellington
subplot(1,3,2)
plot(time_grid/cspice_spd(), sat_elevation.WE*cspice_dpr(),'linewidth',1)
yl = yline(min_elevation.WE,'-','Min El','linewidth',1);
yl.LabelHorizontalAlignment = 'left';
yl.Color = [.85 .33 .10];
xlabel('Epoch [MJD2000]')
ylabel('Elevation [deg]')
xlim([time_grid(1)/cspice_spd(),time_grid(end)/cspice_spd()])
ylim([-90,90])
yticks(-90:30:90)
title('WELLINGTON')
grid on

% La Silla
subplot(1,3,3)
plot(time_grid/cspice_spd(), sat_elevation.LS*cspice_dpr(),'linewidth',1)
yl = yline(min_elevation.LS,'-','Min El','linewidth',1);
yl.LabelHorizontalAlignment = 'left';
yl.Color = [.85 .33 .10];
xlabel('Epoch [MJD2000]')
ylabel('Elevation [deg]')
xlim([time_grid(1)/cspice_spd(),time_grid(end)/cspice_spd()])
ylim([-90,90])
yticks(-90:30:90)
title('LA SILLA')
grid on


%-------------------------------------------------------------------------%

% REQUEST 2

% Initial covariance matrix [km^2, km^2/s, km^2/s^2]
P0 = [ 2.6e-2  1.4e-2 -1.8e-3    0       0       0    ;
       1.4e-2  1.8e-2  2.3e-3    0       0       0    ;
      -1.8e-3  2.3e-3  1.0e-2    0       0       0    ;
         0       0       0     1.6e-7    0       0    ;
         0       0       0       0     1.6e-7    0    ;
         0       0       0       0       0     1.6e-7 ];

% Propagate the initial mean and covariance to the last epoch of each 
% visibility window using a linearized approach (LinCov)
[x_hat_lincov,P_lincov,sqrt_trace_lincov] = lincov(x0_mean,P0,t0_et,tf_win_et,mu);

% Propagate the initial mean and covariance to the last epoch of each 
% visibility window using the Unscented Transform (UT)
[Y_hat_UT,P_UT,sqrt_trace_UT] = UT(x0_mean,P0,t0_et,tf_win_et,mu);

% Plot the square roots of the traces of position and velocity covariance
% submatrices for each method
figure()
subplot(2,2,1)
plot(tf_win_et.MI/cspice_spd(),sqrt_trace_lincov.pos.MI,'ob',tf_win_et.WE/cspice_spd(),sqrt_trace_lincov.pos.WE,'or',tf_win_et.LS/cspice_spd(),sqrt_trace_lincov.pos.LS,'og','linewidth',4)
grid on
legend('MILANO','WELLINGTON','LA-SILLA','Location','SouthEast')
xlabel('Epoch [MJD2000]')
ylabel('$$\sqrt{Tr(P)}$$ [km]','Interpreter','latex')
title('Position covariance submatrix (lincov)')

subplot(2,2,3)
plot(tf_win_et.MI/cspice_spd(),sqrt_trace_lincov.vel.MI,'ob',tf_win_et.WE/cspice_spd(),sqrt_trace_lincov.vel.WE,'or',tf_win_et.LS/cspice_spd(),sqrt_trace_lincov.vel.LS,'og','linewidth',4)
grid on
legend('MILANO','WELLINGTON','LA-SILLA','Location','SouthEast')
xlabel('Epoch [MJD2000]')
ylabel('$$\sqrt{Tr(P)}$$ [km/s]','Interpreter','latex')
title('Velocity covariance submatrix (lincov)')

subplot(2,2,2)
plot(tf_win_et.MI/cspice_spd(),sqrt_trace_UT.pos.MI,'ob',tf_win_et.WE/cspice_spd(),sqrt_trace_UT.pos.WE,'or',tf_win_et.LS/cspice_spd(),sqrt_trace_UT.pos.LS,'og','linewidth',4)
grid on
legend('MILANO','WELLINGTON','LA-SILLA','Location','SouthEast')
xlabel('Epoch [MJD2000]')
ylabel('$$\sqrt{Tr(P)}$$ [km]','Interpreter','latex')
title('Position covariance submatrix (UT)')

subplot(2,2,4)
plot(tf_win_et.MI/cspice_spd(),sqrt_trace_UT.vel.MI,'ob',tf_win_et.WE/cspice_spd(),sqrt_trace_UT.vel.WE,'or',tf_win_et.LS/cspice_spd(),sqrt_trace_UT.vel.LS,'og','linewidth',4)
grid on
legend('MILANO','WELLINGTON','LA-SILLA','Location','SouthEast')
xlabel('Epoch [MJD2000]')
ylabel('$$\sqrt{Tr(P)}$$ [km/s]','Interpreter','latex')
title('Velocity covariance submatrix (UT)')

% Display results of requests 1 and 2 in the command window
disp(' ')
disp('<strong>REQUESTS 1 and 2</strong>')
disp(' ')
disp('<strong> Station: MILANO </strong> ')

for i = 1:length(t0_win_et.MI)
fprintf('<strong>  Visibility window: </strong> %0.f \n',i)
fprintf('Start visibility epoch: %s UTC\n',vis_wind_date.MI(i,:))
fprintf('End visibility epoch: %s UTC\n',vis_wind_date.MI(i+length(t0_win_et.MI),:))

fprintf('Estimated state (lincov): [')
fprintf('%f, ',x_hat_lincov.MI(1:5,i))
fprintf('%f]'' km, km/s\n',x_hat_lincov.MI(end,i))

F = [repmat(' %d',1,6),'\n'];
fprintf('Estimated covariance matrix (lincov) [km^2,km^2/s,km^2/s^2]:\n')
fprintf(F,P_lincov.MI(:,:,i).')

fprintf('Estimated state (UT): [')
fprintf('%f, ',Y_hat_UT.MI(1:5,i))
fprintf('%f]'' km, km/s\n',Y_hat_UT.MI(end,i))

F = [repmat(' %d',1,6),'\n'];
fprintf('Estimated covariance matrix (UT) [km^2,km^2/s,km^2/s^2]:\n')
fprintf(F,P_UT.MI(:,:,i).')
end

disp(' ')
disp('<strong> Station: WELLINGTON </strong> ')

for i = 1:length(t0_win_et.WE)
fprintf('<strong>  Visibility window: </strong> %0.f \n',i)
fprintf('Start visibility epoch: %s UTC\n',vis_wind_date.WE(i,:))
fprintf('End visibility epoch: %s UTC\n',vis_wind_date.WE(i+length(t0_win_et.WE),:))

fprintf('Estimated state (lincov): [')
fprintf('%f, ',x_hat_lincov.WE(1:5,i))
fprintf('%f]'' km, km/s\n',x_hat_lincov.WE(end,i))

fprintf('Estimated covariance matrix (lincov) [km^2,km^2/s,km^2/s^2]:\n')
fprintf(F,P_lincov.WE(:,:,i).');

fprintf('Estimated state (UT): [')
fprintf('%f, ',Y_hat_UT.WE(1:5,i))
fprintf('%f]'' km, km/s\n',Y_hat_UT.WE(end,i))

F = [repmat(' %d',1,6),'\n'];
fprintf('Estimated covariance matrix (UT) [km^2,km^2/s,km^2/s^2]:\n')
fprintf(F,P_UT.WE(:,:,i).')
end

disp(' ')
disp('<strong> Station: LA SILLA </strong> ')

for i = 1:length(t0_win_et.LS)
fprintf('<strong>  Visibility window: </strong> %0.f \n',i)
fprintf('Start visibility epoch: %s UTC\n',vis_wind_date.LS(i,:))
fprintf('End visibility epoch: %s UTC\n',vis_wind_date.LS(i+length(t0_win_et.LS),:))

fprintf('Estimated state (lincov): [')
fprintf('%f, ',x_hat_lincov.LS(1:5,i))
fprintf('%f]'' km, km/s\n',x_hat_lincov.LS(end,i))

fprintf('Estimated covariance matrix (lincov) [km^2,km^2/s,km^2/s^2]:\n')
fprintf(F,P_lincov.LS(:,:,i).');

fprintf('Estimated state (UT): [')
fprintf('%f, ',Y_hat_UT.LS(1:5,i))
fprintf('%f]'' km, km/s\n',Y_hat_UT.LS(end,i))

fprintf('Estimated covariance matrix (UT) [km^2,km^2/s,km^2/s^2]:\n')
fprintf(F,P_UT.LS(:,:,i).');
end

%-------------------------------------------------------------------------%

% REQUEST 3

% Monte Carlo simulation

% Last visibility window is associated to La Silla station
et_last_vis_win = tf_win_et.LS(end); % ephemeris time

% Generate MC samples
n_samples = 100;
R = mvnrnd(x0_mean,P0,n_samples);

% Propagate samples
y = zeros(6,n_samples); % initialize variable
for i=1:n_samples
y(:,i) = flow_2BP(R(i,1:3)',R(i,4:6)',t0_et,et_last_vis_win,mu,1);
end

% Compute sample mean (MC)
x_MC = mean(y,2);

% Compute sample covariance (MC)
P_MC = cov(y');

% Compute square root of the trace of position and velocity covariance
% submatrices (MC)
sqrt_trace_MC.pos = sqrt(trace(P_MC(1:3,1:3)));
sqrt_trace_MC.vel = sqrt(trace(P_MC(4:6,4:6)));

% Display results of Monte Carlo analysis
figure()
subplot(1,2,1)
plot(et_last_vis_win/cspice_spd,sqrt_trace_MC.pos,'ob',et_last_vis_win/cspice_spd,sqrt_trace_lincov.pos.LS(5),'or',et_last_vis_win/cspice_spd,sqrt_trace_UT.pos.LS(5),'og','linewidth',4)
xlabel('Epoch [MJD2000]')
ylabel('$$\sqrt{Tr(P)}$$ [km]','Interpreter','latex')
legend('MC','LINCOV','UT')
title('Position covariance submatrix')
grid on

subplot(1,2,2)
plot(et_last_vis_win/cspice_spd,sqrt_trace_MC.vel,'ob',et_last_vis_win/cspice_spd,sqrt_trace_lincov.vel.LS(5),'or',et_last_vis_win/cspice_spd,sqrt_trace_UT.vel.LS(5),'og','linewidth',4)
xlabel('Epoch [MJD2000]')
ylabel('$$\sqrt{Tr(P)}$$ [km/s]','Interpreter','latex')
legend('MC','LINCOV','UT')
title('Velocity covariance submatrix')
grid on

% Find angle between the topocentric direction of the satellite and the 
% topocentric direction of the center of the field-of-view of the sensor
% for each sample
% Define sensor FoV (circular) 
FoV_LS = 1; % [deg]

% Get from SPK kernel the station state in ECI reference frame
rv_station_eci = cspice_spkezr('LA-SILLA', et_last_vis_win, 'J2000', 'NONE', 'EARTH');

% Get state rotation matrix from ECI (J2000) to TOPOCENTRIC
ROT_ECI2TOPO = cspice_sxform('J2000','LA-SILLA_TOPO', et_last_vis_win);

% Compute the topocentric direction of the center of the field-of-view 
% of the sensor
% Compute reference state at the final visibility window
X_ref = flow_2BP(x0_mean(1:3),x0_mean(4:6),t0_et,et_last_vis_win,mu,1); 
rv_station_sat_ref_eci = X_ref - rv_station_eci;
rv_station_sat_ref_topo = ROT_ECI2TOPO*rv_station_sat_ref_eci; % rotate to topo
dir_sensor_topo = rv_station_sat_ref_topo(1:3)/norm(rv_station_sat_ref_topo(1:3)); % compute direction

% Final state vectors of the satellite in ECI frame for the MC points
rv_sat_eci = y;

% Compute final state vectors of satellite as seen from the station in J2000
rv_station_sat_eci = rv_sat_eci - rv_station_eci;

% Convert ECI state into topocentric
rv_station_sat_topo = zeros(size(rv_station_sat_eci));
angle = zeros(1,n_samples);
for i = 1:size(rv_station_sat_eci,2)
    rv_station_sat_topo(:,i) = ROT_ECI2TOPO*rv_station_sat_eci(:,i);
    % Find angle between the topocentric direction of the satellite and
    % topocentric direction of the center of the field-of-view of the
    % sensor
    angle(:,i) = acos(dot((rv_station_sat_topo(1:3,i)/norm(rv_station_sat_topo(1:3,i))),dir_sensor_topo));
end

% Compute the percentage of samples inside the sensor's FoV
percentage_FoV = sum(rad2deg(angle) < FoV_LS);

% Plot MC points distribution
% Initial distribution
figure()
plot3(R(:,1),R(:,2),R(:,3),'ob')
hold on
plot3(rr0_mean(1),rr0_mean(2),rr0_mean(3),'or','MarkerFaceColor','r')
legend('MC samples','Initial mean state')
xlabel('X [km]')
ylabel('Y [km]')
zlabel('Z [km]')
grid on 
box on

% Final distribution
figure()
plot3(y(1,:),y(2,:),y(3,:),'ob')
hold on
plot3(x_MC(1),x_MC(2),x_MC(3),'or','MarkerFaceColor','r')
legend('Propagated MC points','Final mean state (MC)')
xlabel('X [km]')
ylabel('Y [km]')
zlabel('Z [km]')
grid on 
box on

% Bar plot for the samples that lie inside the FoV
% Plot data
figure()
X_tf1 = categorical({'Samples inside FoV'});
Y_tf1 = percentage_FoV;
b_tf1 = bar(X_tf1,Y_tf1,0.3);
hold on 
X_tf2 = categorical({'Samples outside FoV'});
Y_tf2 = 100 - percentage_FoV;
b_tf2 = bar(X_tf2,Y_tf2,0.3);
% Add percentage value on the tip of the bars
xtips_in = b_tf1.XEndPoints;
ytips_in = b_tf1.YEndPoints;
ylabels_in = string(b_tf1.YData);
text(xtips_in,ytips_in,ylabels_in,'HorizontalAlignment','center','VerticalAlignment','bottom')
xtips_out = b_tf2.XEndPoints;
ytips_out = b_tf2.YEndPoints;
ylabels_out = string(b_tf2.YData);
text(xtips_out,ytips_out,ylabels_out,'HorizontalAlignment','center','VerticalAlignment','bottom')
ylabel('Samples %')
ylim([0,100])
grid on

% Display request 3 results in the command window
disp(' ')
disp('<strong> REQUEST 3: Monte Carlo analysis </strong>')
fprintf('Last visibility epoch: %s UTC \n',vis_wind_date.LS(end,:))
fprintf('Estimated state (MC): [')
fprintf('%f, ',x_MC(1:5))
fprintf('%f]'' km, km/s\n',x_MC(end))
fprintf('Estimated covariance matrix (MC) [km^2,km^2/s,km^2/s^2]:\n')
fprintf(F,P_MC.');

fprintf('Estimated state (lincov): [')
fprintf('%f, ',x_hat_lincov.LS(1:5,end))
fprintf('%f]'' km, km/s\n',x_hat_lincov.LS(end,end))
fprintf('Estimated covariance matrix (lincov) [km^2,km^2/s,km^2/s^2]:\n')
fprintf(F,P_lincov.LS(:,:,end).');

fprintf('Estimated state (UT): [')
fprintf('%f, ',Y_hat_UT.LS(1:5,end))
fprintf('%f]'' km, km/s\n',Y_hat_UT.LS(end,end))
fprintf('Estimated covariance matrix (UT) [km^2,km^2/s,km^2/s^2]:\n')
fprintf(F,P_UT.LS(:,:,end).');

disp(' ')
fprintf('Sqrt of the trace of position covariance submatrix (lincov): %f km \n',sqrt_trace_lincov.pos.LS(end))
fprintf('Sqrt of the trace of position covariance submatrix (UT): %f km \n', sqrt_trace_UT.pos.LS(end))
fprintf('Sqrt of the trace of position covariance submatrix (MC): %f km \n', sqrt_trace_MC.pos)
fprintf('Sqrt of the trace of velocity covariance submatrix (lincov): %f km \n',sqrt_trace_lincov.vel.LS(end))
fprintf('Sqrt of the trace of velocity covariance submatrix (UT): %f km \n', sqrt_trace_UT.vel.LS(end))
fprintf('Sqrt of the trace of velocity covariance submatrix (MC): %f km \n', sqrt_trace_MC.vel)
disp(' ')
fprintf('Percentage of the samples inside sensor FoV: %0.f %% \n',percentage_FoV)

%-------------------------------------------------------------------------%

%-------------------------------------------------------------------------%

%% EXERCISE 2
%-------------------------------------------------------------------------%
%                                                                         %
%                       Run only after running Ex 1                       %
%                                                                         %
%-------------------------------------------------------------------------%

disp(' <strong> Exercise 2 </strong>')
disp(' ')

%-------------------------------------------------------------------------%
% Set options for SGP4
typerun    = 'u';  % user-provided inputs to SGP4 Matlab function
opsmode    = 'a';  % afspc approach ('air force space command')
whichconst =  72;  % WGS72 constants (radius, gravitational parameter)

% Define TLE
longstr1 = '1 28922U 05051A   21323.60115338 -.00000088  00000-0  00000+0 0  9998';
longstr2 = '2 28922  58.3917  37.3090 0004589  75.2965 284.7978  1.69421077 98469';

% Get satrec from TLE
satrec = twoline2rv(longstr1, longstr2, typerun,'e',opsmode, whichconst);

% Get TLE epoch
[year,mon,day,hr,min,sec] = invjday(satrec.jdsatepoch, satrec.jdsatepochf);
sat_epoch_str = sprintf('%d-%02d-%02dT%02d:%02d:%02.6f', [year,mon,day,hr,min,sec]);
sat_epoch_et = cspice_str2et(sat_epoch_str);

fprintf('Satellite num ID: %d\n', satrec.satnum);
fprintf('TLE reference epoch: UTC %s\n', sat_epoch_str);

% Compute minutes since initial epoch
t_grid_min = (time_grid-t0_et)/60; % [min]

% arcsec to radians conversion factor
arcsec2rad = pi/(180*3600);

% Get the delta-Psi and delta-Epsilon corrections from EOP-Last5years.txt
% and convert them from arcsed to radians
ddpsi = -0.112684*arcsec2rad; %  [rad]
ddeps = -0.006350*arcsec2rad; %  [rad]

% Initialize variables to improve computational speed
rteme = zeros(3,length(t_grid_min));
vteme = zeros(3,length(t_grid_min));
ttt = zeros(1,length(t_grid_min));
reci = zeros(3,length(t_grid_min));
veci = zeros(3,length(t_grid_min));

% Propagate satellite state using SGP4
for i = 1:length(t_grid_min) 
% Propagate satellite state in TEME reference frame 
[satrec,rteme(:,i),vteme(:,i)] = sgp4(satrec, t_grid_min(i));
% Compute centuries from TDT 2000 January 1 00:00:00.000
ttt(i) = cspice_unitim(time_grid(i), 'ET', 'TDT')/cspice_jyear()/100;
% Put a dummy value for the acceleration in TEME
ateme = zeros(3,length(t_grid_min));
% Convert from TEME to ECI
[reci(:,i), veci(:,i), ~] = teme2eci(rteme(:,i), vteme(:,i), ateme(:,i), ttt(i), ddpsi, ddeps);
end

% Define measurements noise
sigma_az = 0.1;      % [deg]
sigma_el = 0.1;      % [deg]
sigma_range = 0.01;  % [km]
sigma_ra.WE = 5e-4;  % [deg]
sigma_dec.WE = 5e-4; % [deg]
sigma_ra.LS = 1e-3;  % [deg]
sigma_dec.LS = 1e-3; % [deg]

% Define measurements noise matrices
Sigma.MI = diag([sigma_az^2, sigma_el^2, sigma_range^2]);
Sigma.WE = diag([sigma_ra.WE^2, sigma_dec.WE^2]);
Sigma.LS = diag([sigma_ra.LS^2, sigma_dec.LS^2]);

% Compute the ideal and noisy measurements over the given visibility 
% windows for each station
[meas,meas_noisy] = measurements(reci,veci,time_grid,vis_wind,Sigma);

% Define final epochs until which simulated measurements will be used
tf1 = '2021-Nov-21 14:00:00 UTC';
tf2 = '2021-Nov-22 14:00:00 UTC';
tf3 = '2021-Nov-23 14:00:00 UTC';

% Convert date stings to ephemeris times
tf1_et = cspice_str2et(tf1);
tf2_et = cspice_str2et(tf2);
tf3_et = cspice_str2et(tf3); 

% Plot measurements
% Milano: Az, El, Range
% Az
figure()
subplot(3,1,1)
plot(time_grid/cspice_spd,meas_noisy.az.MI,'.','linewidth',1.5)
xlabel('Epoch [MJD2000]')
ylabel('Azimuth [deg]')
xl1 = xline(tf1_et/cspice_spd,'-','tf_1','linewidth',1);
xl1.LabelHorizontalAlignment = 'left';
xl1.LabelOrientation = 'horizontal';
xl2 = xline(tf2_et/cspice_spd,'-','tf_2','linewidth',1);
xl2.LabelHorizontalAlignment = 'left';
xl2.LabelOrientation = 'horizontal';
xlim([t0_et/cspice_spd,tf_et/cspice_spd])
title('MILANO')
grid on
% El
subplot(3,1,2)
plot(time_grid/cspice_spd,meas_noisy.el.MI,'.','Color',[0.85,0.33,0.10],'linewidth',1.5)
xlabel('Epoch [MJD2000]')
ylabel('Elevation [deg]')
xl1 = xline(tf1_et/cspice_spd,'-','tf_1','linewidth',1);
xl1.LabelHorizontalAlignment = 'left';
xl1.LabelOrientation = 'horizontal';
xl2 = xline(tf2_et/cspice_spd,'-','tf_2','linewidth',1);
xl2.LabelHorizontalAlignment = 'left';
xl2.LabelOrientation = 'horizontal';
xlim([t0_et/cspice_spd,tf_et/cspice_spd])
title('MILANO')
grid on
% Range
subplot(3,1,3)
plot(time_grid/cspice_spd,meas_noisy.range.MI,'.','Color',[0.93,0.69,0.13],'linewidth',1.5)
xlabel('Epoch [MJD2000]')
ylabel('Range [km]')
xl1 = xline(tf1_et/cspice_spd,'-','tf_1','linewidth',1);
xl1.LabelHorizontalAlignment = 'left';
xl1.LabelOrientation = 'horizontal';
xl2 = xline(tf2_et/cspice_spd,'-','tf_2','linewidth',1);
xl2.LabelHorizontalAlignment = 'left';
xl2.LabelOrientation = 'horizontal';
xlim([t0_et/cspice_spd,tf_et/cspice_spd])
title('MILANO')
grid on

% Wellington: Ra, Dec
% Ra
figure()
subplot(2,1,1)
plot(time_grid/cspice_spd,meas_noisy.ra.WE,'.')
xlabel('Epoch [MJD2000]')
ylabel('Right Ascension [deg]')
xl1 = xline(tf1_et/cspice_spd,'-','tf_1','linewidth',1);
xl1.LabelHorizontalAlignment = 'left';
xl1.LabelOrientation = 'horizontal';
xl2 = xline(tf2_et/cspice_spd,'-','tf_2','linewidth',1);
xl2.LabelHorizontalAlignment = 'left';
xl2.LabelOrientation = 'horizontal';
xlim([t0_et/cspice_spd,tf_et/cspice_spd])
ylim([0,360])
title('WELLINGTON')
grid on
% Dec
subplot(2,1,2)
plot(time_grid/cspice_spd,meas_noisy.dec.WE,'.','Color',[0.85,0.33,0.10])
xlabel('Epoch [MJD2000]')
ylabel('Declination [deg]')
xl1 = xline(tf1_et/cspice_spd,'-','tf_1','linewidth',1);
xl1.LabelHorizontalAlignment = 'left';
xl1.LabelOrientation = 'horizontal';
xl2 = xline(tf2_et/cspice_spd,'-','tf_2','linewidth',1);
xl2.LabelHorizontalAlignment = 'left';
xl2.LabelOrientation = 'horizontal';
xlim([t0_et/cspice_spd,tf_et/cspice_spd])
title('WELLINGTON')
grid on

% La Silla: Ra, Dec
% Ra
figure()
subplot(2,1,1)
plot(time_grid/cspice_spd,meas_noisy.ra.LS,'.')
xlabel('Epoch [MJD2000]')
ylabel('Right Ascension [deg]')
title('LA SILLA')
xl1 = xline(tf1_et/cspice_spd,'-','tf_1','linewidth',1);
xl1.LabelHorizontalAlignment = 'left';
xl1.LabelOrientation = 'horizontal';
xl2 = xline(tf2_et/cspice_spd,'-','tf_2','linewidth',1);
xl2.LabelHorizontalAlignment = 'left';
xl2.LabelOrientation = 'horizontal';
xlim([t0_et/cspice_spd,tf_et/cspice_spd])
ylim([0,360])
grid on
% Dec
subplot(2,1,2)
plot(time_grid/cspice_spd,meas_noisy.dec.LS,'.','Color',[0.85,0.33,0.10])
xl1 = xline(tf1_et/cspice_spd,'-','tf_1','linewidth',1);
xl1.LabelHorizontalAlignment = 'left';
xl1.LabelOrientation = 'horizontal';
xl2 = xline(tf2_et/cspice_spd,'-','tf_2','linewidth',1);
xl2.LabelHorizontalAlignment = 'left';
xl2.LabelOrientation = 'horizontal';
xlim([t0_et/cspice_spd,tf_et/cspice_spd])
xlabel('Epoch [MJD2000]')
ylabel('Declination [deg]')
title('LA SILLA')
grid on

%-------------------------------------------------------------------------%

% REQUEST 2

% Define time spans over which the measurements will be used
tspan_tf1 = (t0_et:60:tf1_et);
tspan_tf2 = (t0_et:60:tf2_et);
tspan_tf3 = (t0_et:60:tf3_et);

% Define function handles of the cost function for all the different cases:
% Case 1: measurements used until 2021-Nov-21 14:00:00 UTC, unperturbed
% motion, no a priori info
fun_tf1_unpert_noapriori = @(x) cost_fun(x,vis_wind,meas_noisy,tspan_tf1,mu,x0_mean,P0,'Unperturbed',0);
% Case 2: measurements used until 2021-Nov-21 14:00:00 UTC, J2 perturbed
% motion, no a priori info
fun_tf1_J2_noapriori = @(x) cost_fun(x,vis_wind,meas_noisy,tspan_tf1,mu,x0_mean,P0,'J2',0);
% Case 3: measurements used until 2021-Nov-21 14:00:00 UTC, unperturbed
% motion, with a priori info
fun_tf1_unpert_apriori = @(x) cost_fun(x,vis_wind,meas_noisy,tspan_tf1,mu,x0_mean,P0,'Unperturbed',1);
% Case 4: measurements used until 2021-Nov-21 14:00:00 UTC, J2 perturbed
% motion, with a priori info
fun_tf1_J2_apriori = @(x) cost_fun(x,vis_wind,meas_noisy,tspan_tf1,mu,x0_mean,P0,'J2',1);
% Case 5: measurements used until 2021-Nov-22 14:00:00 UTC, unperturbed
% motion, no a priori info
fun_tf2_unpert_noapriori = @(x) cost_fun(x,vis_wind,meas_noisy,tspan_tf2,mu,x0_mean,P0,'Unperturbed',0);
% Case 6: measurements used until 2021-Nov-22 14:00:00 UTC, J2 perturbed
% motion, no a priori info
fun_tf2_J2_noapriori = @(x) cost_fun(x,vis_wind,meas_noisy,tspan_tf2,mu,x0_mean,P0,'J2',0);
% Case 7: measurements used until 2021-Nov-23 14:00:00 UTC, unperturbed
% motion, no a priori info
fun_tf3_unpert_noapriori = @(x) cost_fun(x,vis_wind,meas_noisy,tspan_tf3,mu,x0_mean,P0,'Unperturbed',0);
% Case 8: measurements used until 2021-Nov-23 14:00:00 UTC, J2 perturbed
% motion, no a priori info
fun_tf3_J2_noapriori = @(x) cost_fun(x,vis_wind,meas_noisy,tspan_tf3,mu,x0_mean,P0,'J2',0);

% Define cell array containing all the function handles
fun_handles = {fun_tf1_unpert_noapriori; fun_tf1_J2_noapriori; ...
                   fun_tf1_unpert_apriori; fun_tf1_J2_apriori;  ...
                   fun_tf2_unpert_noapriori; fun_tf2_J2_noapriori; ...
                   fun_tf3_unpert_noapriori; fun_tf3_J2_noapriori};

% Solve the navigation problem for the different cases               
disp('Find the least squares (minimum variance) solution to the navigation problem')

% Set lsqnonlin options
opt = optimoptions('lsqnonlin', 'Algorithm', 'levenberg-marquardt', 'Display', 'iter');

% Initialize variables for computational speed
x_ls = zeros(6,length(fun_handles));
resnorm = zeros(1,length(fun_handles));
P_ls = zeros(6,6,length(fun_handles));
sqrt_trace_ls_pos = zeros(1,length(fun_handles));
sqrt_trace_ls_vel = zeros(1,length(fun_handles));

for i=1:length(fun_handles)
    fprintf('Case %0.f:\n',i)
    % Solve nonlinear least-squares problem
    [x_ls(:,i), resnorm(i), residual, exitflag, ~, ~, jac] = lsqnonlin(fun_handles{i}, x0_mean, [], [], opt);
    % Compute resulting covariance
    Jac = full(jac);
    P_ls(:,:,i) = resnorm(i) / (length(residual)-length(x0_mean)) .* inv(Jac'*Jac);
    % Compute sqrt of the trace of position and velocity covariance submatrices
    sqrt_trace_ls_pos(i) = sqrt(trace(P_ls(1:3,1:3,i)));
    sqrt_trace_ls_vel(i) = sqrt(trace(P_ls(4:6,4:6,i)));
    disp(' ')
end

% Display results of Exercise 2 in the command window
for i=1:length(fun_handles)
    if i==1 
       disp('<strong> Request 2a </strong>')
       fprintf('Case %0.f:\n',i)
       disp('Measurements used until 2021-Nov-21 14:00:00 UTC, unperturbed motion, no a priori info.')
    elseif i==2
        fprintf('Case %0.f:\n',i)
        disp('Measurements used until 2021-Nov-21 14:00:00 UTC, J2 perturbed motion, no a priori info.')
    elseif i==3
        disp('<strong> Request 2b </strong>')
        fprintf('Case %0.f:\n',i)
        disp('Measurements used until 2021-Nov-21 14:00:00 UTC, unperturbed motion, with a priori info.')
    elseif i==4
        fprintf('Case %0.f:\n',i)
        disp('Measurements used until 2021-Nov-21 14:00:00 UTC, J2 perturbed motion, with a priori info.')
    elseif i==5
        disp('<strong> Request 2c </strong>')
        fprintf('Case %0.f:\n',i)
        disp('Measurements used until 2021-Nov-22 14:00:00 UTC, unperturbed motion, no a priori info.')
    elseif i==6
        fprintf('Case %0.f:\n',i)
        disp('Measurements used until 2021-Nov-22 14:00:00 UTC, J2 perturbed motion, no a priori info.')
    elseif i==7
        fprintf('Case %0.f:\n',i)
        disp('Measurements used until 2021-Nov-23 14:00:00 UTC, unperturbed motion, no a priori info.')
    elseif i==8
        fprintf('Case %0.f:\n',i)
        disp('Measurements used until 2021-Nov-23 14:00:00 UTC, J2 perturbed motion, no a priori info.')
    end
fprintf('Estimated state: [')
fprintf('%f, ',x_ls(1:5,i))
fprintf('%f]'' km, km/s\n',x_ls(end,i))
fprintf('Covariance matrix lsq [km^2,km^2/s,km^2/s^2]: \n');
fprintf(F,P_ls(:,:,i));
fprintf('Squared norm of the residual: %g \n',resnorm(i))
fprintf('Square root of the trace of position covariance submatrix: %f km\n',sqrt_trace_ls_pos(i))
fprintf('Square root of the trace of velocity covariance submatrix: %f km/s\n',sqrt_trace_ls_vel(i))
disp(' ')
end

% Plot square roots of the traces of position and velocity covariance 
% submatrices for the different cases considered
% No a priori VS a priori
figure()
subplot(4,1,1)
plot(t0_et/cspice_spd,sqrt_trace_ls_pos(1),'o',t0_et/cspice_spd,sqrt_trace_ls_pos(3),'o','linewidth',4)
legend('No a priori info','With a priori info')
xlabel('Epoch [MJD2000]')
ylabel('$$\sqrt{Tr(P)}$$ [km]','Interpreter','latex')
title('Position covariance submatrix (unperturbed motion)')
grid on

subplot(4,1,3)
plot(t0_et/cspice_spd,sqrt_trace_ls_pos(2),'o',t0_et/cspice_spd,sqrt_trace_ls_pos(4),'o','linewidth',4)
legend('No a priori info','With a priori info')
ylabel('$$\sqrt{Tr(P)}$$ [km]','Interpreter','latex')
title('Position covariance submatrix (J2 perturbed motion)')
grid on

subplot(4,1,2)
plot(t0_et/cspice_spd,sqrt_trace_ls_vel(1),'o',t0_et/cspice_spd,sqrt_trace_ls_vel(3),'o','linewidth',4)
legend('No a priori info','With a priori info')
xlabel('Epoch [MJD2000]')
ylabel('$$\sqrt{Tr(P)}$$ [km/s]','Interpreter','latex')
title('Velocity covariance submatrix (unperturbed motion)')
grid on

subplot(4,1,4)
plot(t0_et/cspice_spd,sqrt_trace_ls_vel(2),'o',t0_et/cspice_spd,sqrt_trace_ls_vel(4),'o','linewidth',4)
xlabel('Epoch [MJD2000]')
ylabel('$$\sqrt{Tr(P)}$$ [km/s]','Interpreter','latex')
legend('No a priori info','With a priori info')
title('Velocity covariance submatrix (J2 perturbed motion)')
grid on

% Plot tf1 VS tf2 VS tf3
figure()
subplot(4,1,1)
plot(t0_et/cspice_spd,sqrt_trace_ls_pos(1),'o',t0_et/cspice_spd,sqrt_trace_ls_pos(5),'o',t0_et/cspice_spd,sqrt_trace_ls_pos(7),'o','linewidth',4)
legend('tf1','tf2','tf3')
xlabel('Epoch [MJD2000]')
ylabel('$$\sqrt{Tr(P)}$$ [km]','Interpreter','latex')
title('Position covariance submatrix (unperturbed motion)')
grid on

subplot(4,1,2)
plot(t0_et/cspice_spd,sqrt_trace_ls_vel(1),'o',t0_et/cspice_spd,sqrt_trace_ls_vel(5),'o',t0_et/cspice_spd,sqrt_trace_ls_vel(7),'o','linewidth',4)
legend('tf1','tf2','tf3')
xlabel('Epoch [MJD2000]')
ylabel('$$\sqrt{Tr(P)}$$ [km/s]','Interpreter','latex')
title('Velocity covariance submatrix (unperturbed motion)')
grid on

subplot(4,1,3)
plot(t0_et/cspice_spd,sqrt_trace_ls_pos(2),'o',t0_et/cspice_spd,sqrt_trace_ls_pos(6),'o',t0_et/cspice_spd,sqrt_trace_ls_pos(8),'o','linewidth',4)
legend('tf1','tf2','tf3')
xlabel('Epoch [MJD2000]')
ylabel('$$\sqrt{Tr(P)}$$ [km]','Interpreter','latex')
title('Position covariance submatrix (J2 perturbed motion)')
grid on

subplot(4,1,4)
plot(t0_et/cspice_spd,sqrt_trace_ls_vel(2),'o',t0_et/cspice_spd,sqrt_trace_ls_vel(6),'o',t0_et/cspice_spd,sqrt_trace_ls_vel(8),'o','linewidth',4)
legend('tf1','tf2','tf3')
xlabel('Epoch [MJD2000]')
ylabel('$$\sqrt{Tr(P)}$$ [km/s]','Interpreter','latex')
title('Velocity covariance submatrix (J2 perturbed motion)')
grid on

% Compute error wrt SGP4 propagation for unperturbed and J2 

% Perform the Integration in the unperturbed 2BP and J2 perturbed 2BP 
[~,rv_ref] = twobody_propagator(rr0_mean,vv0_mean,mu,tspan_tf3,'Unperturbed');
[~,rv_ref_J2] = twobody_propagator(rr0_mean,vv0_mean,mu,tspan_tf3,'J2');

% Retrieve the state propagated with SGP4
rv_sgp4 = [reci; veci]';

% Compute the orbit propagation error 
err = zeros(1,length(rv_sgp4)); % Prellocate variables for speed
err_J2 = zeros(1,length(rv_sgp4));
for i=1:length(tspan_tf3)
err(i) = norm(rv_sgp4(i,1:3)-rv_ref(i,1:3));
err_J2(i) = norm(rv_sgp4(i,1:3)-rv_ref_J2(i,1:3));
end

% Plot the error
figure()
plot(tspan_tf3/cspice_spd,err,tspan_tf3/cspice_spd,err_J2,'linewidth',1.5)
xl1 = xline(tf1_et/cspice_spd,'-','tf_1','linewidth',1);
xl1.LabelHorizontalAlignment = 'left';
xl1.LabelOrientation = 'horizontal';
xl2 = xline(tf2_et/cspice_spd,'-','tf_2','linewidth',1);
xl2.LabelHorizontalAlignment = 'left';
xl2.LabelOrientation = 'horizontal';
xlim([t0_et/cspice_spd,tf_et/cspice_spd])
legend('Unperturbed','J2','Location','NorthWest')
xlabel('Epoch [MJD2000]')
ylabel('Error [km]')
grid on

%-------------------------------------------------------------------------%

%-------------------------------------------------------------------------%

%% EXERCISE 3

% Clear memory workspace and set path
clearvars; close all; clc
addpath('tdm')

disp(' <strong> Exercise 3 </strong>')
disp(' ')

%-------------------------------------------------------------------------%

% Input the tdm file names
tdm_filename = string;
tdm_filename(1) = 'NEW_NORCIA_EXM_20161010T041435_20161010T150333.tdm';
tdm_filename(2) = 'MALARGUE_EXM_20161010T162859_20161011T033236.tdm';
tdm_filename(3) = 'NEW_NORCIA_EXM_20161011T041336_20161011T150208.tdm';
tdm_filename(4) = 'MALARGUE_EXM_20161011T162803_20161012T033108.tdm';
tdm_filename(5) = 'NEW_NORCIA_EXM_20161012T041238_20161012T150043.tdm';
tdm_filename(6) = 'MALARGUE_EXM_20161012T162708_20161013T032940.tdm';
tdm_filename(7) = 'NEW_NORCIA_EXM_20161013T041139_20161013T145918.tdm';
tdm_filename(8) = 'MALARGUE_EXM_20161013T162614_20161014T032811.tdm';
tdm_filename(9) = 'NEW_NORCIA_EXM_20161014T041041_20161014T145753.tdm';
tdm_filename(10) = 'MALARGUE_EXM_20161014T162519_20161015T032642.tdm';
tdm_filename(11) = 'NEW_NORCIA_EXM_20161015T040944_20161015T145627.tdm';
tdm_filename(12) = 'MALARGUE_EXM_20161015T162426_20161016T032513.tdm';

% Retrieve measurements data from tdm files

% Preallocate variables for computational speed
Meas.tdm_data = [];
Meas.Az = [];
Meas.El = [];
Meas.Range = [];
Meas.ET = [];
Meas.Station_Name = [];
for i=1:numel(tdm_filename) % For each tdm file
% Read the tdm file and extract all useful data 
% (station name, start time, stop time, measurements times and values)
Meas.tdm_data = [Meas.tdm_data; read_tdm(tdm_filename(i))];
% Retrieve Azimuth [deg], Elevation [deg] and Range [km] measurements
Meas.Az = [Meas.Az; Meas.tdm_data(i).Az];
Meas.El = [Meas.El; Meas.tdm_data(i).El];
Meas.Range = [Meas.Range; Meas.tdm_data(i).Range];
% Retrieve measurement ephemeris time
Meas.ET = [Meas.ET; Meas.tdm_data(i).Measurement_ET];
% Retrieve station name
Meas.Station_Name = [Meas.Station_Name; Meas.tdm_data(i).Station_Name];
end

% Set counterclockwise convention for Azimuth measurements
Meas.Az = wrapTo360(-Meas.Az(:)); 

% Reorganize all the measurements into a single array
for i=1:numel(tdm_filename)
Meas.tdm_data(i).All = [wrapTo360(-Meas.tdm_data(i).Az)'; Meas.tdm_data(i).El'; Meas.tdm_data(i).Range'];
end
Meas.All = [0, 0, 0];
for i = 1:length(Meas.Az)
Meas.All = [Meas.All; Meas.Az(i), Meas.El(i), Meas.Range(i)];
end
% Rearrange by columns 
Meas.All = Meas.All'; 
Meas.Station_Name = [0; Meas.Station_Name];

% Define measurements noise
sigma_az = 1.5e-3;   % [deg]
sigma_el = 1.3e-3;   % [deg]
sigma_range = 0.075; % [km]
rho_az_el = 0.1;
rho_az_range = 0;
rho_el_range = 0;

% Define measurements noise matrix
R = [sigma_az^2, rho_az_el*sigma_az*sigma_el, rho_az_range*sigma_az*sigma_range;
     rho_az_el*sigma_az*sigma_el, sigma_el^2, rho_el_range*sigma_el*sigma_range;
     rho_az_range*sigma_az*sigma_range, rho_el_range*sigma_el*sigma_range, sigma_range^2];

% Reference epoch
t0 = '2016-10-10T00:00:00.000'; % UTC
t0_et = cspice_str2et(t0); % ephemeris time

% Define measurements ephemeris time vector
Meas.ET = [t0_et; Meas.ET];

% Final epoch
t_format = 'YYYY-MM-DDTHR:MN:SC.###::UTC';
tf = cspice_timout(Meas.ET(end),t_format); % convert from ET to date string

% Mean state: position [km] and velocity [km/s]
rr0_mean = [+1.68467660241E+08 -1.07050886902E+08 -5.47243873455E+07]';
vv0_mean = [+1.34362486580E+01 +1.68723391839E+01 +8.66147058235E+00]';
x0_mean = [rr0_mean; vv0_mean];

% Initial covariance [km^2, km^2/s, km^2/s^2]
P0 = [+2.01E+04 -7.90E+00 -4.05E+00 -5.39E-03 +6.37E-06 +3.26E-06;
      -7.90E+00 +2.01E+04 +2.64E+00 +6.37E-06 -5.38E-03 -2.07E-06;
      -4.05E+00 +2.64E+00 +2.01E+04 +3.25E-06 -2.03E-06 -5.38E-03;
      -5.39E-03 +6.37E-06 +3.25E-06 +1.92E-07 -2.28E-09 -1.16E-09;
      +6.37E-06 -5.38E-03 -2.03E-06 -2.28E-09 +1.91E-07 +7.31E-10;
      +3.26E-06 -2.07E-06 -5.38E-03 -1.16E-09 +7.31E-10 +1.91E-07];

% Sun standard gravitational parameter [km^3/s^2]
mu_sun = cspice_bodvrd('SUN','GM',1);

% REQUEST 1
% Use an unscented Kalman filter to update sequentially the spacecraft state
tk_1 = t0_et;
Pk(:,:,1) = P0;
xk(:,1) = x0_mean;
Station_Name = string;
for i = 2:length(Meas.ET)
    % Retrieve sequentially all measurements, ET and station name
    yk = Meas.All(:,i);
    tk = Meas.ET(i);
    Station_Name(i) = Meas.Station_Name(i);
    % Process the measurements with an Unscented Kalman Filter to update
    % sequentially the spacecraft state
    [xk(:,i),Pk(:,:,i)] = UTKalmanFilter(xk(:,i-1),Pk(:,:,i-1),tk_1,tk,yk,R,convertStringsToChars(Station_Name(i)),mu_sun);
    tk_1 = tk;    
end

%-------------------------------------------------------------------------%

% REQUEST 2

% Compute the error along the tracking campaign between the estimated mean 
% states and the true trajectory

% Set target ID (ExoMars)
target_ID = '-143';
% Retrieve true ExoMars trajectory 
rv_EXM = cspice_spkezr(target_ID,Meas.ET','J2000','NONE','SUN');

% Compute the error between the estimated xk and the true trajectory
% Preallocate variables for computational speed
err_EXM = zeros(1,length(rv_EXM));
err_EXM_pos = zeros(1,length(rv_EXM));
err_EXM_vel = zeros(1,length(rv_EXM));
for i=1:length(rv_EXM)
err_EXM(i) = norm(rv_EXM(:,i) - xk(:,i));
err_EXM_pos(i) = norm(rv_EXM(1:3,i) - xk(1:3,i));
err_EXM_vel(i) = norm(rv_EXM(4:6,i) - xk(4:6,i));
end

%-------------------------------------------------------------------------%

% REQUEST 3

% Compute the square roots of the traces of the position and velocity 
% covariance submatrices 
% Preallocate variables for computational speed
sqrt_trace = zeros(1,length(Meas.ET));
sqrt_trace_pos = zeros(1,length(Meas.ET));
sqrt_trace_vel = zeros(1,length(Meas.ET));
for i = 1:length(Meas.ET)
sqrt_trace(i) = sqrt(trace(Pk(:,:,i)));
sqrt_trace_pos(i) = sqrt(trace(Pk(1:3,1:3,i)));
sqrt_trace_vel(i) = sqrt(trace(Pk(4:6,4:6,i)));
end

% Plot the covariance matrices and the error wrt true trajectory
figure()
% subplot(2,1,1)
plot(Meas.ET/cspice_spd,err_EXM_pos,Meas.ET/cspice_spd,sqrt_trace_pos,Meas.ET/cspice_spd,3*sqrt_trace_pos,'linewidth',1.8)
xlabel('Epoch [MJD2000]')
xlim([t0_et/cspice_spd, Meas.ET(end)/cspice_spd])
yl1 = yline(75,'-','75 km','linewidth',1);
ylabel('Position error [km]')
legend('ExoMars error','Covariance trace root','3\sigma')
grid on

figure()
% subplot(2,1,2)
plot(Meas.ET/cspice_spd,err_EXM_vel,Meas.ET/cspice_spd,sqrt_trace_vel,Meas.ET/cspice_spd,3*sqrt_trace_vel,'linewidth',1.8)
xlabel('Epoch [MJD2000]')
xlim([t0_et/cspice_spd, Meas.ET(end)/cspice_spd])
yl2 = yline(2e-4,'-','2E-04 km/s','linewidth',1);
ylabel('Velocity error [km/s]')
legend('ExoMars error','Covariance trace root','3\sigma')
grid on

% Display results of Exercise 3 in the command window
fprintf('Initial epoch: UTC %s \n',t0)
fprintf('Tracking campaign final epoch: UTC %s \n',tf)
fprintf('Number of available measurements (Az, El, Range): %0.f \n', length(Meas.Az))
disp(' ')
fprintf('Square root of the trace of the position covariance submatrix at the last measurement epoch: \n %f km\n',sqrt_trace_pos(end))
fprintf('Square root of the trace of the velocity covariance submatrix at the last measurement epoch: \n %f km/s\n',sqrt_trace_vel(end))


% Clear kernel pool
cspice_kclear

% Check that kernel pool is clear
fprintf('\nTotal kernels number after kclear: %d\n', cspice_ktotal('ALL'));

%-------------------------------------------------------------------------%

%-------------------------------------------------------------------------%

%% FUNCTIONS

% function: car2kep.m

function [a,e,i,OM,om,th] = car2kep(rr,vv,mu) 
%-------------------------------------------------------------------------%
%
% car2kep.m transforms cartesian coordinates into keplerian parameters.
%
%-------------------------------------------------------------------------%
% PROTOTYPE:
%  [a,e,i,OM,om,th] = car2kep(rr,vv,mu) 
%
%-------------------------------------------------------------------------%
% INPUT ARGUMENTS:
%  rr            [3x1]   Position vector                   [km]
%  vv            [3x1]   Velocity vector                   [km/s]
%  mu            [1]     Standard gravitational parameter  [km^3/s^2]
%                        of the primary body
%
%-------------------------------------------------------------------------%
% OUTPUT ARGUMENTS:
%  a             [1]     Semi-major axis                   [km]
%  e             [1]     Eccentricity                      [-]
%  i             [1]     Inclination                       [rad]
%  OM            [1]     RAAN                              [rad]
%  om            [1]     Argument of periapsis             [rad]
%  th            [1]     True anomaly                      [rad]
%
% ------------------------------------------------------------------------%
% AUTHOR:
%  Gian Marco Paldino
%
%-------------------------------------------------------------------------%
% VERSIONS: 
%  23/10/2021: First version
%
%-------------------------------------------------------------------------%

% Semi-major axis is calculated
r = norm(rr);
v = norm(vv);

E = (1/2*v^2)-(mu/r); % E = Specific mechanic energy

a = -(mu/(2*E));

% Eccentricity is calculated

hh = cross(rr,vv); % hh = Angular momentum

h = norm(hh);

ee = ((cross(vv,hh))/mu)-(rr/r); % ee = Eccentricity vector

e = norm(ee);

% Inclination is calculated

i = acos(hh(3)/h);

% RAAN is calculated

kk = [0 0 1]';

NN = ((cross(kk,hh))/(norm(cross(kk,hh)))); % NN = Node line

if NN(2) >= 0
    OM = acos(NN(1));
else 
    OM = 2*pi-acos(NN(1));
end

% Argument of periapsis is calculated 

if ee(3) >= 0
    om = acos(dot(NN,ee)/e);
else 
    om = 2*pi-acos(dot(NN,ee)/e);
end

% True anomaly is calculated

vr = dot(vv,rr)/r;

if vr >= 0
    th = acos(dot(rr,ee)/(r*e));
else
    th = 2*pi-acos(dot(rr,ee)/(r*e));
end

end

%-------------------------------------------------------------------------%    

%-------------------------------------------------------------------------%

% function: odetwobody.m

function Dx = odetwobody(~,x,mu)
%-------------------------------------------------------------------------%
%
% odetwobody.m provides the odefun for the integration of the 
% two-body problem using Cartesian coordinates.
%
%-------------------------------------------------------------------------%
% PROTOTYPE:
%  Dx = odetwobody(t,x,mu)
%
%-------------------------------------------------------------------------%
% INPUT ARGUMENTS:
%  t                [1]    Time                              [s]
%  x                [6x1]  State of the system:              [km, km/s]
%                          Position and velocity            
%  mu               [1]    Standard gravitational parameter  [km^3/s^2]
%                          of the primary body
%                   
%-------------------------------------------------------------------------%
% OUTPUT ARGUMENTS: 
%  Dx               [6x1]   Derivative of the state          [km/s, km/s^2]
% 
%-------------------------------------------------------------------------%
% AUTHOR:
%  Gian Marco Paldino
%
%-------------------------------------------------------------------------%
% VERSIONS: 
%  23/10/2021: First version
%
%-------------------------------------------------------------------------%

% Get position vector [km] and velocity vector [km/s] from the state vector
rr = x(1:3);
vv = x(4:6);

% Norm of the position vector
r = norm(rr);

% Get x, y, z components of the position vector [km]
x = rr(1);
y = rr(2);
z = rr(3);      

% Compute the derivative of the state
Dx = [vv(1); vv(2); vv(3); -mu/r^3 * x; -mu/r^3 * y; -mu/r^3 * z];
    
end

%-------------------------------------------------------------------------%

%-------------------------------------------------------------------------%

% function: odetwobody.m

function Dx = odetwobody_J2(~,x,mu)
%-------------------------------------------------------------------------%
%
% odetwobody.m provides the odefun for the integration of the 
% two-body problem using Cartesian coordinates accounting for J2
% perturbation effect.
%
%-------------------------------------------------------------------------%
% PROTOTYPE:
%  Dx = odetwobody(t,x,mu)
%
%-------------------------------------------------------------------------%
% INPUT ARGUMENTS:
%  t                [1]    Time                              [s]
%  x                [6x1]  State of the system:              [km, km/s]
%                          Position and velocity            
%  mu               [1]    Standard gravitational parameter  [km^3/s^2]
%                          of the primary body
%                   
%-------------------------------------------------------------------------%
% OUTPUT ARGUMENTS: 
%  Dx               [6x1]   Derivative of the state          [km/s, km/s^2]
% 
%-------------------------------------------------------------------------%
% AUTHOR:
%  Gian Marco Paldino
%
%-------------------------------------------------------------------------%
% VERSIONS: 
%  23/10/2021: First version
%
%-------------------------------------------------------------------------%

% Get position vector [km] and velocity vector [km/s] from the state vector
rr = x(1:3);
vv = x(4:6);

% Norm of the position vector
r = norm(rr);

% Get x, y, z components of the position vector [km]
x = rr(1);
y = rr(2);
z = rr(3);    

% Earth radius [km] and second zonal harmonic J2
R_e = cspice_bodvrd( 'EARTH', 'RADII', 3);
R_e = R_e(1);
J2 = 0.00108263;

% Compute perturbing acceleration due to J2 in cartesian coordinates
kJ2 = (1.5*J2*mu*R_e^2)/r^4;
a_J2_x = kJ2 * (x/r)*(5*(z^2/r^2)-1);
a_J2_y = kJ2 * (y/r)*(5*(z^2/r^2)-1);
a_J2_z = kJ2 * (z/r)*(5*(z^2/r^2)-3);

% Compute the derivative of the state
Dx = [  vv(1)                   ;
        vv(2)                   ;
        vv(3)                   ;
        -mu/r^3 * x  +  a_J2_x  ;
        -mu/r^3 * y  +  a_J2_y  ;
        -mu/r^3 * z  +  a_J2_z  ];
    
end

%-------------------------------------------------------------------------%

%-------------------------------------------------------------------------%

% function: twobody_propagator.m

function [T,X] = twobody_propagator(rr0,vv0,mu,tspan,PerturbationType)
%-------------------------------------------------------------------------%
%
% orbit_propagation.m performs the numerical integration of the equations 
% of motion for the 2BP. 
%
%-------------------------------------------------------------------------%
% PROTOTYPE:
%  [T,X] = twobody_propagator(rr0,vv0,mu,tspan,'Unperturbed')
%  [T,X] = twobody_propagator(rr0,vv0,mu,tspan,'J2')
%
%-------------------------------------------------------------------------%
% INPUT ARGUMENTS:
%  rr0              [3x1]  Initial position                     [km]
%  vv0              [3x1]  Initial velocity                     [km/s]
%  mu               [1]    Standard gravitational parameter     [km^3/s^2]
%                          of the primary body
%  t_span           [1xn]  Integration time steps vector        [s]  
%  PerturbationType [char] Type of perturbation of the motion:                        
%                          'Unperturbed' for unperturbed motion
%                          'J2'          for J2 perturbation
%
%-------------------------------------------------------------------------%
% OUTPUT ARGUMENTS: 
%  T               [1xn]   Integrated orbit time steps          [s]
%  X               [6x1]   State of the system:                 [km, km/s]
%                          Position and velocity            
%
%-------------------------------------------------------------------------%
% AUTHOR:
%  Gian Marco Paldino
%
%-------------------------------------------------------------------------%
% VERSIONS: 
%  23/10/2021: First version
%
%-------------------------------------------------------------------------%

% Set integration initial conditions
X0 = [rr0 ; vv0];

% Set integrator options
ode_options = odeset( 'RelTol', 1e-12, 'AbsTol', 1e-14 );

switch PerturbationType
    
    case 'Unperturbed'
    % Perform the Integration in the unperturbed 2BP
    [T, X] = ode45( @(t,x) odetwobody(t,x,mu), tspan, X0, ode_options);
    case 'J2'
    % Perform the Integration in the J2 perturbed 2BP 
    [T, X] = ode45( @(t,x) odetwobody_J2(t,x,mu), tspan, X0, ode_options); 
end

end

%-------------------------------------------------------------------------%

%-------------------------------------------------------------------------%

% function: PlotEarth.m

function PlotEarth
%-------------------------------------------------------------------------%
%
% plotEarth.m plots the Earth inside a 3d graphic.
%
% ------------------------------------------------------------------------%
% AUTHOR:
%  Gian Marco Paldino
%
%-------------------------------------------------------------------------%
% VERSIONS: 
%  23/10/2021: first version
%
%-------------------------------------------------------------------------%

C = imread('map.jpg');
theta = 0;
[x,y,z] = ellipsoid(0,0,0, 6378.135, 6378.135,6356.750,1E2);
surf(x,y,z,circshift(flip(C),[0,ceil(size(C,2)/360*theta)]),'FaceColor','texturemap','EdgeColor','none');
axis equal;
end

%-------------------------------------------------------------------------%

%-------------------------------------------------------------------------%

% function: flow_2BP.m

function Xt = flow_2BP(rr0,vv0,ti,tf,mu,outflag)
%-------------------------------------------------------------------------%
%
% flow_2BP.m computes the flow of the 2BP.
%
%-------------------------------------------------------------------------%
% PROTOTYPE:
%  Xt = flow_2BP(rr0,vv0,ti,tf,mu,0)
%  Xt = flow_2BP(rr0,vv0,ti,tf,mu,1)
%
%-------------------------------------------------------------------------%
% INPUT ARGUMENTS:
%  rr0             [3x1]  Initial position                   [km]
%  vv0             [3x1]  Initial velocity                   [km/s]
%  ti              [1]    Initial time                       [s]
%  tf              [1]    Final time                         [s]
%  mu              [1]    Standard gravitational parameter   [km^3/s^2]
%                         of the primary body
%  outflag         [1]    Output flag: position only (0)     [-]
%                         or full state (1)     
%                   
%-------------------------------------------------------------------------%
% OUTPUT ARGUMENTS: 
%  Xt              [6x1]  Flow of the 2BP at time tf        [km, km/s]
% 
%-------------------------------------------------------------------------%
% AUTHOR:
%  Gian Marco Paldino
%
%-------------------------------------------------------------------------%
% VERSIONS: 
%  29/10/2021: First version
%
%-------------------------------------------------------------------------%

% Set integration timespan 
tspan = [ti tf];

% Set integration options
flow_options = odeset('RelTol', 1e-12, 'AbsTol', 1e-14);

% Set initial conditions for the integration of the ODE  
X0 = [rr0; vv0];

% Perform the integration
[~, X] = ode45(@(t,x) odetwobody(t,x,mu), tspan, X0, flow_options);

% Retrieve the state (flow) at the last time step 
% Decide whether to retrieve position only or full state (position and
% velocity)
if nargin>4
    switch outflag
        case 0  % Position only as output
             Xt = X(end,1:3)';  
        case 1  % Full state as output
             Xt = X(end,:)';
        otherwise    
           warning('Input 0 for position only or 1 for full state output')
    end
end
    
end

%-------------------------------------------------------------------------%

%-------------------------------------------------------------------------%

% visibility_windows.m 

function [t0_win_et, tf_win_et,vis_wind,vis_wind_et,vis_wind_date,ind_t0_win,ind_tf_win] = visibility_windows(sat_elevation,min_elevation,time_grid)
%-------------------------------------------------------------------------%
%
% visibility_windows.m computes the visibility windows for Milano, 
% Wellington and La Silla stations.
%
%-------------------------------------------------------------------------%
% PROTOTYPE:
% [t0_win_et, tf_win_et,vis_wind,vis_wind_et,vis_wind_date,ind_t0_win,ind_tf_win] = visibility_windows(sat_elevation,min_elevation,time_grid)
%
%-------------------------------------------------------------------------%
% INPUT ARGUMENTS:
%  sat_elevation  [1x1 struct]  Satellite elevation               [deg]           
%  min_elevation  [1x1 struct]  Minimum elevation                 [deg]
%  time_grid      [1xn]         Time span (Ephemeris Time)        [s]
%
%  (NOTE: n = number of points in the time span)
%
%-------------------------------------------------------------------------%
% OUTPUT ARGUMENTS: 
%  t0_win_et      [struct]   Initial epochs of the visibility     [s]
%                            window (Ephemeris Time)
%  tf_win_et      [struct]   Final epochs of the visibility       [s]
%                            window (Ephemeris Time)
%  vis_wind       [struct]   Array with value:                    [-]
%                             1   if corresponding epoch belongs 
%                                 to the visibility window
%                             NaN if corresponding epoch does not  
%                                 belong to the visibility window
%  vis_wind_et    [struct]   Initial and final epochs of the      [s]
%                            visibility window (Ephemeris Time)
%  vis_wind_date  [struct]   Initial and final epochs of the      [-]
%                            visibility window (date string)
%  ind_t0_win     [struct]   Index associated to initial epochs   [-]
%                            of each visibility window
%  ind_tf_win     [struct]   Index associated to final epochs     [-]
%                            of each visibility window
%
%-------------------------------------------------------------------------%
% Each structure has the following fields: 'MI' for Milano station
%                                          'WE' for Wellington station
%                                          'LS' for La Silla station
%
%-------------------------------------------------------------------------%
% AUTHOR:
%  Gian Marco Paldino
%
%-------------------------------------------------------------------------%
% VERSIONS: 
%  30/11/2021: First version
%
%-------------------------------------------------------------------------%

% Preallocate variables before for loop to improve computational speed
vis_wind.MI = NaN(1,length(time_grid));
vis_wind.WE = NaN(1,length(time_grid));
vis_wind.LS = NaN(1,length(time_grid));
% If the epoch is in the visibility window --> vis_wind = 1
% If the epoch is not in the visibility window --> vis_wind = NaN
for k = 1:length(time_grid)
    if rad2deg(sat_elevation.MI(k)) > (min_elevation.MI)
    vis_wind.MI(k) = 1;
    end
    if rad2deg(sat_elevation.WE(k)) > (min_elevation.WE)
    vis_wind.WE(k) = 1;
    end
    if rad2deg(sat_elevation.LS(k)) > (min_elevation.LS)
    vis_wind.LS(k) = 1;
    end
end

% Preallocate variables before for loop to improve computational speed
flag_t0.MI = zeros(1,length(sat_elevation.MI));
flag_t0.WE = zeros(1,length(sat_elevation.WE));
flag_t0.LS = zeros(1,length(sat_elevation.LS));
flag_tf.MI = zeros(1,length(sat_elevation.MI));
flag_tf.WE = zeros(1,length(sat_elevation.WE));
flag_tf.LS = zeros(1,length(sat_elevation.LS));

for i=1:length(time_grid)-1
% If flag_t0 == 1 --> its position in the array is associated to the index 
% of the initial epoch of the visibility window
    for j = i+1
        if rad2deg(sat_elevation.MI(i)) < min_elevation.MI && ...
           rad2deg(sat_elevation.MI(j)) > min_elevation.MI
            flag_t0.MI(j) = 1;
        else 
            flag_t0.MI(j) = 0;
        end
        if rad2deg(sat_elevation.WE(i)) < min_elevation.WE && ...
           rad2deg(sat_elevation.WE(j)) > min_elevation.WE
            flag_t0.WE(j) = 1;
        else 
            flag_t0.WE(j) = 0;
        end
        if rad2deg(sat_elevation.LS(i)) < min_elevation.LS && ...
           rad2deg(sat_elevation.LS(j)) > min_elevation.LS
            flag_t0.LS(j) = 1;
        else 
            flag_t0.LS(j) = 0;
        end
% If flag_tf == 1 --> its position in the array is associated to the index 
% of the final epoch of the visibility window
        if rad2deg(sat_elevation.MI(i)) > min_elevation.MI && ...
           rad2deg(sat_elevation.MI(j)) < min_elevation.MI
            flag_tf.MI(i) = 1;
        else 
            flag_tf.MI(i) = 0;
        end
        if rad2deg(sat_elevation.WE(i)) > min_elevation.WE && ...
           rad2deg(sat_elevation.WE(j)) < min_elevation.WE
            flag_tf.WE(i) = 1;
        else 
            flag_tf.WE(i) = 0;
        end
        if rad2deg(sat_elevation.LS(i)) > min_elevation.LS && ...
           rad2deg(sat_elevation.LS(j)) < min_elevation.LS
            flag_tf.LS(i) = 1;
        else 
            flag_tf.LS(i) = 0;
        end
    end
end

% Retrieve index associated to initial epoch of each visibility window
ind_t0_win = struct('MI',find(flag_t0.MI==1),'WE',find(flag_t0.WE==1),...
                    'LS',find(flag_t0.LS==1));

% Retrieve initial epoch of each visibility window (Ephemeris Time)              
t0_win_et = struct('MI',time_grid(ind_t0_win.MI),'WE',time_grid(ind_t0_win.WE),...
                   'LS',time_grid(ind_t0_win.LS));
               
% Retrieve index associated to final epoch of each visibility window
ind_tf_win = struct('MI',find(flag_tf.MI==1),'WE',find(flag_tf.WE==1),...
                    'LS',find(flag_tf.LS==1));

% Retrieve final epoch of each visibility window (Ephemeris Time)              
tf_win_et = struct('MI',time_grid(ind_tf_win.MI),'WE',time_grid(ind_tf_win.WE),...
                   'LS',time_grid(ind_tf_win.LS));
               
% Assemble array with ET of the initial epoch of each visibility window 
% in the first column and ET of the final epoch of each visibility window 
% in the second column
vis_wind_et = struct('MI',[t0_win_et.MI; tf_win_et.MI]',...
                  'WE',[t0_win_et.WE; tf_win_et.WE]',...
                  'LS',[t0_win_et.LS; tf_win_et.LS]');
              
% Set time format for SPICE time conversions
t_format = 'DD-Mon-YYYY HR:MN:SC.####::UTC';          

% Convert ET to date string for the initial epoch of each visibility window
t0_win_date = struct('MI',cspice_timout(t0_win_et.MI,t_format),...
                     'WE',cspice_timout(t0_win_et.WE,t_format),...
                     'LS',cspice_timout(t0_win_et.LS,t_format));
                 
% Convert ET to date string for the final epoch of each visibility window              
tf_win_date = struct('MI',cspice_timout(tf_win_et.MI,t_format),...
                     'WE',cspice_timout(tf_win_et.WE,t_format),...
                     'LS',cspice_timout(tf_win_et.LS,t_format));
                 
% Date strings of the initial and final epochs for each visibility window
vis_wind_date = struct('MI',[t0_win_date.MI; tf_win_date.MI],...
                  'WE',[t0_win_date.WE; tf_win_date.WE],...
                  'LS',[t0_win_date.LS; tf_win_date.LS]);
               
end

%-------------------------------------------------------------------------%

%-------------------------------------------------------------------------%

% stm.m

function STM = stm(x0,t0,tf,mu)
%-------------------------------------------------------------------------%
%
% stm.m computes the State Transition Matrix of the 2BP characterized by
% pure Keplerian motion.
%
%-------------------------------------------------------------------------%
% PROTOTYPE:
%  STM = stm(x0,t0,tf,mu)
%
%-------------------------------------------------------------------------%
% INPUT ARGUMENTS:
%  x               [6x1]  State of the system:               [km, km/s]
%                          Position and velocity            
%  t0              [1]    Initial time                       [s]
%  tf              [1]    Final time                         [s]
%  mu              [1]    Standard gravitational parameter   [km^3/s^2]
%                         of the primary body  
%                   
%-------------------------------------------------------------------------%
% OUTPUT ARGUMENTS: 
%  STM              [6x6]  State Transition Matrix           [km, km/s]
% 
%-------------------------------------------------------------------------%
% AUTHOR:
%  Gian Marco Paldino
%
%-------------------------------------------------------------------------%
% VERSIONS: 
%  30/11/2021: First version
%
%-------------------------------------------------------------------------%

% Retrieve position and velocity components from the initial state
x = x0(1);
y = x0(2);
z = x0(3);
vx = x0(4);
vy = x0(5);
vz = x0(6);

% Compute the flow of the 2BP starting from the initial unperturbed state
flow = flow_2BP(x0(1:3),x0(4:6),t0,tf,mu,1);

% Define the perturbations
eps_x0 = [sqrt(eps)*max(1,abs(x)); 0; 0; 0; 0; 0];
eps_y0 = [0; sqrt(eps)*max(1,abs(y)); 0; 0; 0; 0];
eps_z0 = [0; 0; sqrt(eps)*max(1,abs(z)); 0; 0; 0];
eps_vx0 = [0; 0; 0; sqrt(eps)*max(1,abs(vx)); 0; 0];
eps_vy0 = [0; 0; 0; 0; sqrt(eps)*max(1,abs(vy)); 0];
eps_vz0 = [0; 0; 0; 0; 0; sqrt(eps)*max(1,abs(vz))];

% Perturb the initial condition
x0_pert_x0 = x0 + eps_x0;
x0_pert_y0 = x0 + eps_y0;
x0_pert_z0 = x0 + eps_z0;
x0_pert_vx0 = x0 + eps_vx0;
x0_pert_vy0 = x0 + eps_vy0;
x0_pert_vz0 = x0 + eps_vz0;

% Compute the flow of the 2BP using the perturbed initial conditions
PHI_eps_x0 = flow_2BP(x0_pert_x0(1:3),x0_pert_x0(4:6),t0,tf,mu,1);
PHI_eps_y0 = flow_2BP(x0_pert_y0(1:3),x0_pert_y0(4:6),t0,tf,mu,1);
PHI_eps_z0 = flow_2BP(x0_pert_z0(1:3),x0_pert_z0(4:6),t0,tf,mu,1);
PHI_eps_vx0 = flow_2BP(x0_pert_vx0(1:3),x0_pert_vx0(4:6),t0,tf,mu,1);
PHI_eps_vy0 = flow_2BP(x0_pert_vy0(1:3),x0_pert_vy0(4:6),t0,tf,mu,1);
PHI_eps_vz0 = flow_2BP(x0_pert_vz0(1:3),x0_pert_vz0(4:6),t0,tf,mu,1);

% Compute the State Transition Matrix
STM = zeros(6,6);
for j = 1:6
    STM(j,1) = (PHI_eps_x0(j)-flow(j))/(sqrt(eps)*max(1,abs(x)));
    STM(j,2) = (PHI_eps_y0(j)-flow(j))/(sqrt(eps)*max(1,abs(y)));
    STM(j,3) = (PHI_eps_z0(j)-flow(j))/(sqrt(eps)*max(1,abs(z)));
    STM(j,4) = (PHI_eps_vx0(j)-flow(j))/(sqrt(eps)*max(1,abs(vx)));
    STM(j,5) = (PHI_eps_vy0(j)-flow(j))/(sqrt(eps)*max(1,abs(vy)));
    STM(j,6) = (PHI_eps_vz0(j)-flow(j))/(sqrt(eps)*max(1,abs(vz)));
end

end

%-------------------------------------------------------------------------%

%-------------------------------------------------------------------------%

% lincov.m

function [x_hat_lincov,P_lincov,sqrt_trace] = lincov(x0_mean,P0,t0_et,tf_win_et,mu)
%-------------------------------------------------------------------------%
%
% lincov.m propagates the mean state and the covariance matrix to the last 
% epoch of each visibility window using a linearized approach.
%
%-------------------------------------------------------------------------%
% PROTOTYPE:
%  [x_hat_lincov,P_lincov,sqrt_trace] = lincov(x0_mean,P0,t0_et,tf_win_et,mu)
%
%-------------------------------------------------------------------------%
% INPUT ARGUMENTS:
%  x0_mean       [6x1]     Initial mean state of the system:    [km, km/s]
%                          Position and velocity     
%  P0            [6x6]     Initial covariance matrix            [km^2, km^2/s, km^2/s^2]  
%  t0_et         [1]       Initial epoch (Ephemeris Time)       [s]
%  tf_win_et     [struct]  Final epoch of each visibility       [s]
%                          window (Ephemeris Time) 
%  mu            [1]       Standard gravitational parameter     [km^3/s^2]
%                          of the primary body  
%                   
%-------------------------------------------------------------------------%
% OUTPUT ARGUMENTS: 
%  x_hat_lincov  [struct]  Propagated mean state of the system: [km, km/s]
%                          Position and velocity     
%  P_lincov      [struct]  Propagated covariance matrix         [km^2, km^2/s, km^2/s^2]  
%  sqrt_trace    [struct]  Square root of the trace of position [km, km/s]
%                          and velocity covariance submatrices
%
%-------------------------------------------------------------------------%
% Each structure has the following fields: 'MI' for Milano station
%                                          'WE' for Wellington station
%                                          'LS' for La Silla station
%
%-------------------------------------------------------------------------%
% AUTHOR:
%  Gian Marco Paldino
%
%-------------------------------------------------------------------------%
% VERSIONS: 
%  01/12/2021: First version
%
%-------------------------------------------------------------------------%

% Compute number of visibility windows for each station
num_vis_win = struct('MI',length(tf_win_et.MI),...
                     'WE',length(tf_win_et.WE),...
                     'LS',length(tf_win_et.LS));
 
% Propagate mean and covariance to the last epoch of each visibility window 
% and compute sqrt of the trace of position and velocity covariance submatrices

% Milano
for n=1:num_vis_win.MI
 % Propagate initial mean 
 x_hat_lincov.MI(:,n) = flow_2BP(x0_mean(1:3),x0_mean(4:6),t0_et,tf_win_et.MI(n),mu,1);
 % Compute State Transition Matrix
 STM.MI(:,:,n) = stm(x0_mean,t0_et,tf_win_et.MI(n),mu);
 % Compute covariance
 P_lincov.MI(:,:,n) = STM.MI(:,:,n)*P0*STM.MI(:,:,n)';
 % Compute covariance trace root
 sqrt_trace.pos.MI(n) = sqrt(trace(P_lincov.MI(1:3,1:3,n)));
 sqrt_trace.vel.MI(n) = sqrt(trace(P_lincov.MI(4:6,4:6,n)));
end
% Wellington
for n=1:num_vis_win.WE
 % Propagate initial mean 
 x_hat_lincov.WE(:,n) = flow_2BP(x0_mean(1:3),x0_mean(4:6),t0_et,tf_win_et.WE(n),mu,1);
 % Compute State Transition Matrix
 STM.WE(:,:,n) = stm(x0_mean,t0_et,tf_win_et.WE(n),mu);
 % Compute covariance
 P_lincov.WE(:,:,n) = STM.WE(:,:,n)*P0*STM.WE(:,:,n)';
 % Compute covariance trace root
 sqrt_trace.pos.WE(n) = sqrt(trace(P_lincov.WE(1:3,1:3,n)));
 sqrt_trace.vel.WE(n) = sqrt(trace(P_lincov.WE(4:6,4:6,n)));
end
% La Silla
for n=1:num_vis_win.LS
 % Propagate initial mean 
 x_hat_lincov.LS(:,n) = flow_2BP(x0_mean(1:3),x0_mean(4:6),t0_et,tf_win_et.LS(n),mu,1);
 % Compute State Transition Matrix
 STM.LS(:,:,n) = stm(x0_mean,t0_et,tf_win_et.LS(n),mu);
 % Compute covariance
 P_lincov.LS(:,:,n) = STM.LS(:,:,n)*P0*STM.LS(:,:,n)';
 % Compute covariance trace root
 sqrt_trace.pos.LS(n) = sqrt(trace(P_lincov.LS(1:3,1:3,n)));
 sqrt_trace.vel.LS(n) = sqrt(trace(P_lincov.LS(4:6,4:6,n)));
end

end

%-------------------------------------------------------------------------%

%-------------------------------------------------------------------------%

% UT.m

function [Y_hat,PY,sqrt_trace] = UT(x0_mean,P0,t0_et,tf_win_et,mu)
%-------------------------------------------------------------------------%
%
% UT.m propagates the mean state and the covariance matrix to the last 
% epoch of each visibility window using the unscented transform method.
%
%-------------------------------------------------------------------------%
% PROTOTYPE:
%  [Y_hat,PY,sqrt_trace] = UT(x0_mean,P0,t0_et,tf_win_et,mu)
%
%-------------------------------------------------------------------------%
% INPUT ARGUMENTS:
%  x0_mean       [6x1]     Initial mean state of the system:  [km, km/s]
%                          Position and velocity     
%  P0            [6x6]     Initial covariance matrix            [km^2, km^2/s, km^2/s^2]  
%  t0_et         [1]       Initial epoch (Ephemeris Time)       [s]
%  tf_win_et     [struct]  Final epoch of each visibility       [s]
%                          window (Ephemeris Time) 
%  mu            [1]       Standard gravitational parameter     [km^3/s^2]
%                          of the primary body  
%                   
%-------------------------------------------------------------------------%
% OUTPUT ARGUMENTS: 
%  Y_hat        [struct]   Propagated mean state of the system: [km, km/s]
%                          Position and velocity     
%  P_UT         [struct]   Propagated covariance matrix         [km^2, km^2/s, km^2/s^2]  
%  sqrt_trace   [struct]   Square root of the trace of position [km, km/s]
%                          and velocity covariance submatrices
%
%-------------------------------------------------------------------------%
% Each structure has the following fields: 'MI' for Milano station
%                                          'WE' for Wellington station
%                                          'LS' for La Silla station
%
%-------------------------------------------------------------------------%
% AUTHOR:
%  Gian Marco Paldino
%
%-------------------------------------------------------------------------%
% VERSIONS: 
%  01/12/2021: First version
%
%-------------------------------------------------------------------------%

% Compute the scaling parameter lambda
alpha = 10^-2; % determines the spread of sigma points
beta = 2; % usually set to 2 for optimality
n = numel(x0_mean); % number of states
k = 0; % usually set to 0
lambda = alpha^2 * (n+k) - n;

% Compute the square root matrix
sqrt_matrix = sqrtm((n+lambda)*P0);

% Compose the sigma points
X0 = x0_mean';
X = zeros(2*n,n);
for i = 1:(2*n)
    if i > n
        X(i,:) = x0_mean' - sqrt_matrix(i-n,:);
    else
        X(i,:) = x0_mean' + sqrt_matrix(i,:);
    end
end

% Compute the weights
W0_m = lambda/(n+lambda);
W0_c = W0_m + (1 + alpha^2 + beta);
Wi = 1/(2*(n+lambda));

% Compute number of visibility windows for each station
num_vis_win = struct('MI',length(tf_win_et.MI),...
                     'WE',length(tf_win_et.WE),...
                     'LS',length(tf_win_et.LS));

% Initialize variables
PY = struct('MI',zeros(6,6,num_vis_win.MI),'WE',zeros(6,6,num_vis_win.WE),'LS',zeros(6,6,num_vis_win.LS));
PYi = struct('MI',zeros(6,6,num_vis_win.MI),'WE',zeros(6,6,num_vis_win.WE),'LS',zeros(6,6,num_vis_win.LS));
summation = struct('MI',zeros(6,6,num_vis_win.MI),'WE',zeros(6,6,num_vis_win.WE),'LS',zeros(6,6,num_vis_win.LS));

% MILANO
for l=1:num_vis_win.MI
    % Propagate sigma points
    Y0.MI(:,l) = flow_2BP(X0(1:3)',X0(4:6)',t0_et,tf_win_et.MI(l),mu,1);
    WY0.MI(:,l) = W0_m*Y0.MI(:,l); % add weights
    for i = 1:(2*n)
        Y.MI(:,i,l) = flow_2BP(X(i,1:3)',X(i,4:6)',t0_et,tf_win_et.MI(l),mu,1);
        WY.MI(:,i,l) = Wi*Y.MI(:,i,l); % add weights
        % Compute weighted sample mean
        Y_hat.MI(:,l) = (WY0.MI(:,l) + sum(WY.MI(:,:,l),2));
    end
    % Compute weighted sample covariance (for i=0)
    PY0.MI(:,:,l) = W0_c*(Y0.MI(:,l)-Y_hat.MI(:,l))*(Y0.MI(:,l)-Y_hat.MI(:,l))';
    % Compute weighted sample covariance (for i=1:12)
    for i = 1:(2*n)
        summation.MI(:,:,l) = Wi*(Y.MI(:,i,l)-Y_hat.MI(:,l))*(Y.MI(:,i,l)-Y_hat.MI(:,l))';
        PYi.MI(:,:,l) = PYi.MI(:,:,l) + summation.MI(:,:,l);
    end
    % Sum to compute overall weighted sample covariance
    PY.MI(:,:,l) = PY0.MI(:,:,l) + PYi.MI(:,:,l);
    
    % Compute square root of the trace of position and velocity covariance
    % submatrices
    sqrt_trace.pos.MI(l) = sqrt(trace(PY.MI(1:3,1:3,l)));
    sqrt_trace.vel.MI(l) = sqrt(trace(PY.MI(4:6,4:6,l)));
end                 
                 
% WELLINGTON
for l=1:num_vis_win.WE
    % Propagate sigma points
    Y0.WE(:,l) = flow_2BP(X0(1:3)',X0(4:6)',t0_et,tf_win_et.WE(l),mu,1);
    WY0.WE(:,l) = W0_m*Y0.WE(:,l); % add weights
    for i = 1:(2*n)
        Y.WE(:,i,l) = flow_2BP(X(i,1:3)',X(i,4:6)',t0_et,tf_win_et.WE(l),mu,1);
        WY.WE(:,i,l) = Wi*Y.WE(:,i,l); % add weights
        % Compute weighted sample mean
        Y_hat.WE(:,l) = (WY0.WE(:,l) + sum(WY.WE(:,:,l),2));
    end
    % Compute weighted sample covariance (for i=0)
    PY0.WE(:,:,l) = W0_c*(Y0.WE(:,l)-Y_hat.WE(:,l))*(Y0.WE(:,l)-Y_hat.WE(:,l))';
    % Compute weighted sample covariance (for i=1:12)
    for i = 1:(2*n)
        summation.WE(:,:,l) = Wi*(Y.WE(:,i,l)-Y_hat.WE(:,l))*(Y.WE(:,i,l)-Y_hat.WE(:,l))';
        PYi.WE(:,:,l) = PYi.WE(:,:,l) + summation.WE(:,:,l);
    end
    % Sum to compute overall weighted sample covariance
    PY.WE(:,:,l) = PY0.WE(:,:,l) + PYi.WE(:,:,l);
    
    % Compute square root of the trace of position and velocity covariance
    % submatrices
    sqrt_trace.pos.WE(l) = sqrt(trace(PY.WE(1:3,1:3,l)));
    sqrt_trace.vel.WE(l) = sqrt(trace(PY.WE(4:6,4:6,l)));
end  

% LA SILLA
for l=1:num_vis_win.LS
    % Propagate sigma points
    Y0.LS(:,l) = flow_2BP(X0(1:3)',X0(4:6)',t0_et,tf_win_et.LS(l),mu,1);
    WY0.LS(:,l) = W0_m*Y0.LS(:,l); % add weights
    for i = 1:(2*n)
        Y.LS(:,i,l) = flow_2BP(X(i,1:3)',X(i,4:6)',t0_et,tf_win_et.LS(l),mu,1);
        WY.LS(:,i,l) = Wi*Y.LS(:,i,l); % add weights
        % Compute weighted sample mean
        Y_hat.LS(:,l) = (WY0.LS(:,l) + sum(WY.LS(:,:,l),2));
    end
    % Compute weighted sample covariance (for i=0)
    PY0.LS(:,:,l) = W0_c*(Y0.LS(:,l)-Y_hat.LS(:,l))*(Y0.LS(:,l)-Y_hat.LS(:,l))';
    % Compute weighted sample covariance (for i=1:12)
    for i = 1:(2*n)
        summation.LS(:,:,l) = Wi*(Y.LS(:,i,l)-Y_hat.LS(:,l))*(Y.LS(:,i,l)-Y_hat.LS(:,l))';
        PYi.LS(:,:,l) = PYi.LS(:,:,l) + summation.LS(:,:,l);
    end
    % Sum to compute overall weighted sample covariance
    PY.LS(:,:,l) = PY0.LS(:,:,l) + PYi.LS(:,:,l);
    
    % Compute square root of the trace of position and velocity covariance
    % submatrices
    sqrt_trace.pos.LS(l) = sqrt(trace(PY.LS(1:3,1:3,l)));
    sqrt_trace.vel.LS(l) = sqrt(trace(PY.LS(4:6,4:6,l)));
    
end  

end

%-------------------------------------------------------------------------%

%-------------------------------------------------------------------------%

% measurements.m

function [meas,meas_noisy] = measurements(reci,veci,time_grid,vis_wind,Sigma)
%-------------------------------------------------------------------------%
%
% measurements.m provides the expected and noisy measurements associated
% to each visibility window.
%
%-------------------------------------------------------------------------%
% PROTOTYPE:
% [meas,meas_noisy] = measurements(reci,veci,time_grid,vis_wind,Sigma)
%
%-------------------------------------------------------------------------%
% INPUT ARGUMENTS:
%  reci          [3xn]     ECI position vectors                [km]    
%  veci          [3xn]     ECI velocity vector s               [km/s]  
%  time_grid     [1xn]     Reference time interval             [s]
%                          (Ephemeris Time)
%  vis_wind      [struct]  Array with value:                   [-]
%                             1   if corresponding epoch  
%                                 belongs to the vis window
%                             NaN if corresponding epoch does   
%                                 not belong to the vis window
%  Sigma         [struct]  Measurements noise matrix           [deg^2,km^2]
%                   
%-------------------------------------------------------------------------%
% OUTPUT ARGUMENTS: 
%  meas          [struct]  Expected measurements:
%                           Az, El, Range for Milano           [deg,deg,km]
%                           Ra, Dec for Wellington             [deg,deg]
%                           Ra, Dec for La Silla               [deg,deg]
%  meas_noisy    [struct]  Noisy measurements:
%                           Az, El, Range for Milano           [deg,deg,km]
%                           Ra, Dec for Wellington             [deg,deg]
%                           Ra, Dec for La Silla               [deg,deg]
%
%-------------------------------------------------------------------------%
% Each structure has the following fields: 'MI' for Milano station
%                                          'WE' for Wellington station
%                                          'LS' for La Silla station
%
%-------------------------------------------------------------------------%
% AUTHOR:
%  Gian Marco Paldino
%
%-------------------------------------------------------------------------%
% VERSIONS: 
%  05/12/2021: First version
%
%-------------------------------------------------------------------------%

% Compute Az, El, Range for MILANO station
[meas.az.MI,meas.el.MI, meas.range.MI] = antenna_pointing('MILANO',time_grid,[reci; veci]);

% Compute Ra, Dec for WELLINGTON and LA SILLA stations
[~,meas.ra.LS,meas.dec.LS] = cspice_recrad(reci);
[~,meas.ra.WE,meas.dec.WE] = cspice_recrad(reci);

% Set measurement to NaN if the corresponding epoch is not in the 
% visibility window of the station
for k = 1:length(vis_wind.MI) % MILANO
    if isnan(vis_wind.MI(k))  
        meas.az.MI(k) = NaN;
        meas.el.MI(k) = NaN;
        meas.range.MI(k) = NaN;
    end
end
for k = 1:length(vis_wind.WE) % WELLINGTON
    if isnan(vis_wind.WE(k))
    meas.ra.WE(k) = NaN; meas.dec.WE(k) = NaN;
    end
end
for k = 1:length(vis_wind.LS) % LA SILLA
    if isnan(vis_wind.LS(k))
    meas.ra.LS(k) = NaN; meas.dec.LS(k) = NaN;
    end
end

% Convert angular measurements from radians to degrees
meas.az.MI = rad2deg(meas.az.MI);
meas.el.MI = rad2deg(meas.el.MI);
meas.ra.WE = rad2deg(meas.ra.WE);
meas.dec.WE = rad2deg(meas.dec.WE);
meas.ra.LS = rad2deg(meas.ra.LS);
meas.dec.LS = rad2deg(meas.dec.LS);

% Define mean measurements for each station
mu_meas.MI = [meas.az.MI', meas.el.MI', meas.range.MI'];
mu_meas.WE = [meas.ra.WE', meas.dec.WE'];
mu_meas.LS = [meas.ra.LS', meas.dec.LS'];

% Generate random noisy measurements for each station
meas_noisy.MI = mvnrnd(mu_meas.MI,Sigma.MI);
meas_noisy.WE = mvnrnd(mu_meas.WE,Sigma.WE);
meas_noisy.LS = mvnrnd(mu_meas.LS,Sigma.LS);

% Extract noisy measurements
meas_noisy.az.MI = meas_noisy.MI(:,1);
meas_noisy.el.MI = meas_noisy.MI(:,2);
meas_noisy.range.MI = meas_noisy.MI(:,3);
meas_noisy.ra.WE = meas_noisy.WE(:,1);
meas_noisy.dec.WE = meas_noisy.WE(:,2);
meas_noisy.ra.LS = meas_noisy.LS(:,1);
meas_noisy.dec.LS = meas_noisy.LS(:,2);

end

%-------------------------------------------------------------------------%

%-------------------------------------------------------------------------%

% cost_fun.m

function residual = cost_fun(x,vis_wind,meas_noisy,tspan,mu,x0_mean,P0,PerturbationType,APrioriInfo)
%-------------------------------------------------------------------------%
%
% measurements.m provides the expected and noisy measurements associated
% to each visibility window.
%
%-------------------------------------------------------------------------%
% PROTOTYPE:
% residual = cost_fun(x,vis_wind,meas_noisy,tspan,mu,x0_mean,P0,'Unperturbed',0)
% residual = cost_fun(x,vis_wind,meas_noisy,tspan,mu,x0_mean,P0,'Unperturbed',1)
% residual = cost_fun(x,vis_wind,meas_noisy,tspan,mu,x0_mean,P0,'J2',0)
% residual = cost_fun(x,vis_wind,meas_noisy,tspan,mu,x0_mean,P0,'J2',1)
%
%-------------------------------------------------------------------------%
% INPUT ARGUMENTS:
%  x                [6x1]     State of the system:                [km,km/s]
%                             Position and velocity 
%  vis_wind         [struct]  Array with value:                   [-]
%                              1   if corresponding epoch  
%                                  belongs to the vis window
%                              NaN if corresponding epoch does   
%                                  not belong to the vis window
%  meas_noisy       [struct]  Noisy measurements:
%                              Az, El, Range for Milano           [deg,deg,km]
%                              Ra, Dec for Wellington             [deg,deg]
%                              Ra, Dec for La Silla               [deg,deg]
%  tspan            [1xn]     Reference time interval             [s]
%                             (Ephemeris Time)
%  mu               [1]       Standard gravitational parameter    [km^3/s^2]
%                             of the Earth
%  x0_mean          [6x1]     Initial mean state of the system:   [km, km/s]
%                             Position and velocity     
%  P0               [6x6]     Initial covariance matrix           [km^2, km^2/s, km^2/s^2] 
%  PerturbationType [char]    'Unperturbed' for unperturbed 2BP   [-]
%                             'J2' for J2 perturbed 2BP
%  APrioriInfo      [1]        0 --> LSQ with no a priori info    [-]
%                              1 --> LSQ with a priori info
%
%-------------------------------------------------------------------------%
% OUTPUT ARGUMENTS: 
%  residual      [nx1]      Weighted difference of the            [deg, km]
%                           measurements
%
%-------------------------------------------------------------------------%
% Each structure has the following fields: 'MI' for Milano station
%                                          'WE' for Wellington station
%                                          'LS' for La Silla station
%
%-------------------------------------------------------------------------%
% AUTHOR:
%  Gian Marco Paldino
%
%-------------------------------------------------------------------------%
% VERSIONS: 
%  05/12/2021: First version
%
%-------------------------------------------------------------------------%

% Retrieve position [km] and velocity [km/s] from the state
rr0 = x(1:3);
vv0 = x(4:6);

% Compute the reference trajectory
switch PerturbationType   
    case 'Unperturbed'
    % Perform the Integration in the unperturbed 2BP
    [~,rv_ref] = twobody_propagator(rr0,vv0,mu,tspan,'Unperturbed');
    case 'J2'
    % Perform the Integration in the J2 perturbed 2BP
     [~,rv_ref] = twobody_propagator(rr0,vv0,mu,tspan,'J2');    
end

% Define measurements noise
sigma_az = 0.1;      % [deg]
sigma_el = 0.1;      % [deg]
sigma_range = 0.01;  % [km]

sigma_ra.WE = 5e-4;  % [deg]
sigma_dec.WE = 5e-4; % [deg]

sigma_ra.LS = 1e-3;  % [deg]
sigma_dec.LS = 1e-3; % [deg]

% Compute weights
W_m.MI = diag([1/sigma_az^2, 1/sigma_el^2, 1/sigma_range^2]);
W_m.WE = diag([1/sigma_ra.WE^2, 1/sigma_dec.WE^2]);
W_m.LS = diag([1/sigma_ra.LS^2, 1/sigma_dec.LS^2]);

% Initialize the residual
residual = [];

% Check if a priori info is present
if APrioriInfo == 1
W_ap = sqrtm(inv(P0));
diff_apriori_weighted = W_ap * (x - x0_mean);
residual = [residual; diff_apriori_weighted(:)];  
elseif APrioriInfo == 0
    residual = [];
else 
    error('Input 1 for LSQ with a priori info or 0 for LSQ without a priori info')
end

% Define time span length
n = length(tspan);

% Compute the residuals
for i = 1:n
    % Milano
    if ~isnan(vis_wind.MI(i)) % If the epoch belongs to the visibility window
    % Compute predicted Az, El, Range 
    [sat_azimuth.MI, sat_elevation.MI, sat_range.MI, ~] = ...
    antenna_pointing('MILANO', tspan(i), rv_ref(i,:)');
    % Rearrange the predicted measurements into column vector
    meas_pred = [rad2deg(sat_azimuth.MI); rad2deg(sat_elevation.MI); sat_range.MI];
    % Retrieve noisy measurements
    meas_real = meas_noisy.MI(i,:)';
    % Compute weighted measurements difference
    diff_meas_weighted = W_m.MI * (meas_pred - meas_real);
    % Compute the residual
    residual = [residual; diff_meas_weighted(:)]; % cannot be preallocated 
                                                  % because its size
                                                  % changes for every case
                                                  % considered
    end
    % Wellington
    if ~isnan(vis_wind.WE(i)) % If the epoch belongs to the visibility window
    % Compute Ra, Dec  
    [~,sat_ra.WE,sat_dec.WE] = cspice_recrad(rv_ref(i,1:3)');    
    % Rearrange the predicted measurements into column vector
    meas_pred = [rad2deg(sat_ra.WE); rad2deg(sat_dec.WE)];
    % Retrieve noisy measurements
    meas_real = meas_noisy.WE(i,:)';
    % Compute weighted measurements difference
    diff_meas_weighted = W_m.WE * (meas_pred - meas_real);
    % Compute the residual
    residual = [residual; diff_meas_weighted(:)];
    end
    % La Silla
    if ~isnan(vis_wind.LS(i)) % If the epoch belongs to the visibility window
    % Compute Ra, Dec
    [~,sat_ra.LS,sat_dec.LS] = cspice_recrad(rv_ref(i,1:3)');    
    % Rearrange the predicted measurements into column vector 
    meas_pred = [rad2deg(sat_ra.LS); rad2deg(sat_dec.LS)];
    % Retrieve noisy measurements
    meas_real = meas_noisy.LS(i,:)';
    % Compute weighted measurements difference
    diff_meas_weighted = W_m.LS * (meas_pred - meas_real);
    % Compute the residual
    residual = [residual; diff_meas_weighted(:)];
    end

end
              
end

%-------------------------------------------------------------------------%

%-------------------------------------------------------------------------%

% read_tdm.m

function Measurements = read_tdm(tdm_filename)
%-------------------------------------------------------------------------%
%
% read_tdm.m retrieves the measurements data from a tdm file.
%
%-------------------------------------------------------------------------%
% PROTOTYPE:
%  Measurements = read_tdm(tdm_filename)
%
%-------------------------------------------------------------------------%
% INPUT ARGUMENTS:
%  tdm_filename   [char]      Input tdm file name
%                   
%-------------------------------------------------------------------------%
% OUTPUT ARGUMENTS: 
%  Measurements   [struct]    Data extracted from the tdm file
%              --> Fields:    'Station_name' [char]   Name of the station
%
%                             'Start_Time'   [struct] Tdm start epoch 
%                                         --> Fields: 'ET' Ephemeris Time
%                                                     'Date' Date string
%
%                             'Stop_Time'    [struct] Tdm stop epoch 
%                                         --> Fields: 'ET' Ephemeris Time
%                                                     'Date' Date string
%
%                              'Data'           [nx4]  ET [s], Az [deg], 
%                                                      El [deg], Range [km]
%
%                              'Measurement_ET' [nx1]  Measurement epoch [s]
%                                                      (Ephemeris Time)
%   
%                              'Az'             [nx1]  Azimuth   [deg]
%                           
%                              'El'             [nx1]  Elevation [deg]
% 
%                              'Range'          [nx1]  Range     [km]
%
% (NOTE: n is the number of measurements contained in the tdm file)
%
%-------------------------------------------------------------------------%
% AUTHOR:
%  Gian Marco Paldino
%
%-------------------------------------------------------------------------%
% VERSIONS: 
%  10/12/2021: First version
%
%-------------------------------------------------------------------------%

% Import the lines of the tdm file
lines = readlines(tdm_filename); 
% Total number of lines in the tdm file
tot_num_lines = length(lines); 

% Retrieve start time
start_time_linenum = 10; % number of the line of the start time
start_time_line = split(lines(start_time_linenum),' = ');
start_time_datestr = convertStringsToChars(start_time_line(2)); % start time date string
start_time_et = cspice_str2et(start_time_datestr); % start time ephemeris time

% Retrieve stop time
stop_time_linenum = 11; % number of the line of the stop time
stop_time_line = split(lines(stop_time_linenum),' = ');
stop_time_datestr = convertStringsToChars(stop_time_line(2)); % stop time date string
stop_time_et = cspice_str2et(stop_time_datestr); % stop time ephemeris time

% Retrieve station name
station_name_linenum = 12; % number of the line of the station name
station_name_line = split(lines(station_name_linenum),' = ');
station_name = convertStringsToChars(station_name_line(2));

% Retrieve measurements
az_linenum_start = 22; % number of the line of the first azimuth measurement
el_linenum_start = 23; % number of the line of the first elevation measurement
range_linenum_start = 24; % number of the line of the first range measurement

az_linenum_end = tot_num_lines-5; % number of the line of the last azimuth measurement
el_linenum_end = tot_num_lines-4; % number of the line of the last elevation measurement
range_linenum_end = tot_num_lines-3; % number of the line of the last range measurement

step = 4; % number of lines after which the same measurement appears again

% Split angle string from measurement time and measurement value string
angle_time_meas_az = [];
angle_time_meas_el = [];
angle_time_meas_range = [];
% Azimuth
for i = az_linenum_start:step:az_linenum_end  
    angle_time_meas_az = [angle_time_meas_az; split(lines(i),' = ')];
                                                  % cannot be preallocated 
                                                  % because its size
                                                  % changes for every TDM
                                                  % file considered                                                 
end
% Elevation
for i = el_linenum_start:step:el_linenum_end
    angle_time_meas_el = [angle_time_meas_el; split(lines(i),' = ')];
end
% Split range string from measurement time and measurement value string
for i = range_linenum_start:step:range_linenum_end
    angle_time_meas_range = [angle_time_meas_range; split(lines(i),' = ')];
end

% Split measurement time string from measurement value string
time_meas_az = [];
time_meas_el = [];
time_meas_range = [];
% Azimuth
for j = 2:2:length(angle_time_meas_az)
    time_meas_az = [time_meas_az; split(angle_time_meas_az(j),' ')];
end
% Elevation
for j = 2:2:length(angle_time_meas_el)
    time_meas_el = [time_meas_el; split(angle_time_meas_el(j),' ')];
end
% Range
for j = 2:2:length(angle_time_meas_range)
    time_meas_range = [time_meas_range; split(angle_time_meas_range(j),' ')];
end

% Get the azimuth measurement string and convert it to double
az = [];
for l = 2:2:length(time_meas_az)
    az = [az; str2double(time_meas_az(l))];
end

% Get the elevation measurement string and convert it to double 
el = [];
for l = 2:2:length(time_meas_el)
    el = [el; str2double(time_meas_el(l))];
end

% Get the range measurement string and convert it to double
range = [];
for l = 2:2:length(time_meas_range)
    range = [range; str2double(time_meas_range(l))];
end

% Get measurement time string and convert it to ephemeris time (it can be 
% retrieved for one measurement only because the time will be the same for 
% all the measurements)
time = [];
for k = 1:2:length(time_meas_az)
    time = [time; cspice_str2et(convertStringsToChars(time_meas_az(k)))];
end

% Store the station name in a vector with length equal to the number of
% measurements
Station_Name = string;
for i = 1:length(time)
    Station_Name(i,:) = station_name;
end

% Assemble output structure
Measurements = struct('Station_Name',Station_Name,...
                     'Start_Time',(struct('ET',start_time_et,'Date',start_time_datestr)),...
                     'Stop_Time',(struct('ET',stop_time_et,'Date',stop_time_datestr)),...
                     'Data',[time, az, el, range],'Measurement_ET',time,...
                     'Az',az,'El',el,'Range',range);
                 
end

%-------------------------------------------------------------------------%

%-------------------------------------------------------------------------%

% UTKalmanFilter.m

function [xk,Pk] = UTKalmanFilter(xk_1,Pk_1,tk_1,tk,yk,R,station_name,mu_sun)
%-------------------------------------------------------------------------%
%
% UTKalmanFilter.m updates sequentially the spacecraft state in terms of
% mean and covariance using an Unscented Kalman Filter.
%
%-------------------------------------------------------------------------%
% PROTOTYPE:
%  [xk,Pk] = UTKalmanFilter(xk_1,Pk_1,tk_1,tk,yk,R,station_name,mu_sun)
%
%-------------------------------------------------------------------------%
% INPUT ARGUMENTS:
%  xk_1          [6x1]     Initial mean state:               [km, km/s]
%                          Position and velocity     
%  Pk_1          [6x6]     Initial covariance matrix         [km^2, km^2/s, km^2/s^2]  
%  tk_1          [1]       Initial epoch (Ephemeris Time)    [s]
%  tk            [1]       Final epoch (Ephemeris Time)      [s]
%  yk            [3x1]     Measurements vector (Az, El, Rng) [deg, deg, km]
%  R             [3x3]     Measurements noise matrix         [deg^2, deg^2*km, km^2]
%  station_name  [char]    Ground Station ID                 [-]
%  mu_sun        [1]       Standard gravitational parameter  [km^3/s^2]
%                          of the Sun  
%                   
%-------------------------------------------------------------------------%
% OUTPUT ARGUMENTS: 
%  xk            [6x1]     Final mean state:                [km, km/s]
%                          Position and velocity     
%  Pk            [6x6]     Final covariance matrix          [km^2, km^2/s, km^2/s^2]
%
%-------------------------------------------------------------------------%
% AUTHOR:
%  Gian Marco Paldino
%
%-------------------------------------------------------------------------%
% VERSIONS: 
%  11/12/2021: First version
%
%-------------------------------------------------------------------------%

% Apply Unscented Transform

% Compute sigma points and weights

% Compute the scaling parameter lambda
alpha = 10^-3; % determines the spread of sigma points
beta = 2; % usually set to 2 for optimality
n = numel(xk_1); % number of states
k = 0; % usually set to 0
lambda = alpha^2 * (n+k) - n;

% Compute the square root matrix
sqrt_matrix = sqrtm((n+lambda)*Pk_1);
sqrt_matrix = nearestSPD((sqrt_matrix+sqrt_matrix')/2);

% Compose the sigma points
X0 = xk_1';
X = zeros(2*n,n);
for i = 1:(2*n)
    if i > n
        X(i,:) = xk_1' - real(sqrt_matrix(i-n,:));
    else
        X(i,:) = xk_1' + real(sqrt_matrix(i,:));
    end
end

% Compute the weights
W0_m = lambda/(n+lambda);
W0_c = W0_m + (1 + alpha^2 + beta);
Wi = 1/(2*(n+lambda));

% Perform the Prediction step

% Propagate the sigma points
% Preallocate variables for computational speed
Xi = zeros(6,2*n);
rv_target_eci_i = zeros(6,2*n);
azimuth = zeros(1,2*n);
elevation = zeros(1,2*n);
range = zeros(1,2*n);
gamma_i = zeros(3,2*n);
% Propagate sigma point (i=0)
Xi_0 = flow_2BP(X0(1:3)',X0(4:6)',tk_1,tk,mu_sun,1);
% Propagate sigma points (i=1:12)
% Compute the state of the Earth in Sun-Centric J2000
rv_earth = cspice_spkezr('EARTH',tk,'J2000','NONE','SUN'); 
for i = 1:(2*n)
        Xi(:,i) = flow_2BP(X(i,1:3)',X(i,4:6)',tk_1,tk,mu_sun,1);
        
        % Derive the predicted measurements
        % Compute spacecraft state wrt Earth in ECI
        rv_target_eci_0 = Xi_0 - rv_earth; 
        rv_target_eci_i(:,i) = Xi(:,i) - rv_earth;
        % Compute Az [deg], El [deg], Range [km]
        [azimuth_0, elevation_0, range_0, ~] = antenna_pointing(station_name, tk, rv_target_eci_0);
        [azimuth(:,i), elevation(:,i), range(:,i), ~] = antenna_pointing(station_name, tk, rv_target_eci_i(:,i));
        % Retrieve predicted measurements vector
        gamma_0 = [wrapTo360(rad2deg(azimuth_0)); rad2deg(elevation_0); range_0];
        gamma_i(:,i) = [wrapTo360(rad2deg(azimuth(:,i))); rad2deg(elevation(:,i)); range(:,i)];        
end

% Compute a priori mean state 
% (slightly better accuracy if all components are summed one by one rather 
%  than using sum command)
xk_hat_ = W0_m*Xi_0 + Wi*Xi(:,1) + Wi*Xi(:,2) + Wi*Xi(:,3) + Wi*Xi(:,4) + ...
          Wi*Xi(:,5) + Wi*Xi(:,6) + Wi*Xi(:,7) + Wi*Xi(:,8) + Wi*Xi(:,9) + ...
          Wi*Xi(:,10) + Wi*Xi(:,11) + Wi*Xi(:,12);
  
% Compute a priori measurement mean
yk_hat_ = W0_m*gamma_0 + Wi*gamma_i(:,1) + Wi*gamma_i(:,2) + Wi*gamma_i(:,3) + ...
          Wi*gamma_i(:,4) + Wi*gamma_i(:,5) + Wi*gamma_i(:,6) + Wi*gamma_i(:,7) + ...
          Wi*gamma_i(:,8) + Wi*gamma_i(:,9) + Wi*gamma_i(:,10) + Wi*gamma_i(:,11) + ...
          Wi*gamma_i(:,12);

% Compute a priori covariance
Pk_0 = W0_c*(Xi_0-xk_hat_)*(Xi_0-xk_hat_)';
Pk_i = zeros(6,6); % preallocate variable
for i = 1:(2*n)
    summation_Pk = Wi*(Xi(:,i)-xk_hat_)*(Xi(:,i)-xk_hat_)';
    Pk_i = Pk_i + summation_Pk;
end
Pk_ = Pk_0 + Pk_i;

% Perform the Update step    

% Compute measurement covariance    
Pyy_k_0 = W0_c*(gamma_0-yk_hat_)*(gamma_0-yk_hat_)';
Pyy_k_i = zeros(3,3); % preallocate variable
for i = 1:(2*n)
    summation_Pyy_k = Wi*(gamma_i(:,i)-yk_hat_)*(gamma_i(:,i)-yk_hat_)';
    Pyy_k_i = Pyy_k_i + summation_Pyy_k;
end
Pyy_k = Pyy_k_0 + Pyy_k_i;
Pyy_k = nearestSPD((Pyy_k+Pyy_k')/2); % make matrix symmetric positive definite
Pyy_k = Pyy_k + R; % Add measurement noise

% Compute cross covariance    
Pxy_k_0 = W0_c*(Xi_0-xk_hat_)*(gamma_0-yk_hat_)';
Pxy_k_i = zeros(6,3); % preallocate variable
for i = 1:(2*n)
    summation_Pxy_k = Wi*(Xi(:,i)-xk_hat_)*(gamma_i(:,i)-yk_hat_)';
    Pxy_k_i = Pxy_k_i + summation_Pxy_k;
end
Pxy_k = Pxy_k_0 + Pxy_k_i;

% Compute the Kalman Gain
Kk = Pxy_k/Pyy_k;
    
% Compute the a posteriori mean and covariance    
xk = xk_hat_ + (Kk*(yk - yk_hat_));
Pk = Pk_ - Kk*Pyy_k*Kk';
Pk = nearestSPD((Pk+Pk')/2);   

end

%-------------------------------------------------------------------------%

%-------------------------------------------------------------------------%

% antenna_pointing.m 

function [sat_azimuth, sat_elevation, sat_range, sat_range_rate] = antenna_pointing(stationName, et, rv_target_eci)

%station_visibility Summary of this function goes here
%Detailed explanation goes here

topoFrame = [stationName, '_TOPO'];

% Get from SPK kernel the station state in ECI reference frame
rv_station_eci = cspice_spkezr(stationName, et, 'J2000', 'NONE', 'EARTH');

% Compute state vector of satellite as seen from the station in J2000
rv_station_sat_eci = rv_target_eci - rv_station_eci;

% Get state rotation matrix from ECI (J2000) to TOPOCENTRIC
ROT_ECI2TOPO = cspice_sxform('J2000',topoFrame, et);

% Convert target ECI state into topocentric
rv_station_sat_topo = zeros(size(rv_station_sat_eci));
for i = 1:size(rv_station_sat_eci,2)
    rv_station_sat_topo(:,i) = ROT_ECI2TOPO(:,:,i)*rv_station_sat_eci(:,i);
end

% Compute range and range-rate
% Compute euclidean norm along direction 1 (column-wise)
sat_range = vecnorm(rv_station_sat_topo(1:3,:),2,1);

% Range rate is the scalar product of the velocity and a unit vector in the
% range direction
sat_range_rate = dot(rv_station_sat_topo(4:6,:), rv_station_sat_topo(1:3,:)./sat_range);

% Compute azimuth and elevation
rll_station_sat = cspice_xfmsta(rv_station_sat_topo,'RECTANGULAR','LATITUDINAL','EARTH');

%fprintf('%.5f %.5f\n',sat_range-rll_station_sat(1,:))
%fprintf('%.5f %.5f\n',sat_range_rate-rll_station_sat(4,:))

sat_azimuth = rll_station_sat(2,:);
sat_elevation = rll_station_sat(3,:);

end

%-------------------------------------------------------------------------%

%-------------------------------------------------------------------------%

% nearestSPD.m

function Ahat = nearestSPD(A)
% nearestSPD - the nearest (in Frobenius norm) Symmetric Positive Definite matrix to A
% usage: Ahat = nearestSPD(A)
%
% From Higham: "The nearest symmetric positive semidefinite matrix in the
% Frobenius norm to an arbitrary real matrix A is shown to be (B + H)/2,
% where H is the symmetric polar factor of B=(A + A')/2."
%
% http://www.sciencedirect.com/science/article/pii/0024379588902236
%
% arguments: (input)
%  A - square matrix, which will be converted to the nearest Symmetric
%    Positive Definite Matrix.
%
% Arguments: (output)
%  Ahat - The matrix chosen as the nearest SPD matrix to A.

if nargin ~= 1
  error('Exactly one argument must be provided.')
end

% test for a square matrix A
[r,c] = size(A);
if r ~= c
  error('A must be a square matrix.')
elseif (r == 1) && (A <= 0)
  % A was scalar and non-positive, so just return eps
  Ahat = eps;
  return
end

% symmetrize A into B
B = (A + A')/2;

% Compute the symmetric polar factor of B. Call it H.
% Clearly H is itself SPD.
[~,Sigma,V] = svd(B);
H = V*Sigma*V';

% get Ahat in the above formula
Ahat = (B+H)/2;

% ensure symmetry
Ahat = (Ahat + Ahat')/2;

% test that Ahat is in fact PD. if it is not so, then tweak it just a bit.
p = 1;
k = 0;
while p ~= 0
  [~,p] = chol(Ahat);
  k = k + 1;
  if p ~= 0
    % Ahat failed the chol test. It must have been just a hair off,
    % due to floating point trash, so it is simplest now just to
    % tweak by adding a tiny multiple of an identity matrix.
    mineig = min(eig(Ahat));
    Ahat = Ahat + (-mineig*k.^2 + eps(mineig))*eye(size(A));
  end
end
end