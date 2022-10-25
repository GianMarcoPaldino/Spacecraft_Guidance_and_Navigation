%-------------------------------------------------------------------------%
%                                                                         %                             
%                    SPACECRAFT GUIDANCE AND NAVIGATION                   %
%                              A.Y. 2021/2022                             %
%                                                                         %
%                        GIAN MARCO PALDINO - 968731                      %
%                                                                         %
%                               ASSIGNMENT 1                              %
%                                                                         %
%-------------------------------------------------------------------------%
% DISCLAIMERS:                                                            %
%          - This code was developed and tested using a Mac/Intel OSX     %
%            64bit machine.                                               %
%          - The folder 'textures' must be added to the path in order to  %
%            correctly visualize planets/Sun/Moon in the plots.           %
%-------------------------------------------------------------------------%
%% EXERCISE 1

% Clear memory workspace and set path
clearvars; close all; clc
addpath(genpath(pwd)); 

disp('Exercise 1')

% Load SPICE kernels: 
% de425s.bsp, gm_de431.tpc, naif0012.tls, pck00010.tpc
cspice_furnsh('assignment01.tm');

% Count total kernels number
fprintf('Total kernels number: %d\n', cspice_ktotal('ALL'));

% REQUEST 1

% Define a generic elliptical orbit
a = 3.4105e+4;               %  Semi-major axis        [km]
e = 0.3226;                  %  Eccentricity           [-]
inc = deg2rad(79.9808);      %  Inclination            [rad]
OM = deg2rad(300);           %  RAAN                   [rad]
om = deg2rad(10);            %  Argument of perigee    [rad]
th = deg2rad(0);             %  True anomaly           [rad]

% Earth gravitational constant [km^3/s^2]
mu_earth = cspice_bodvrd('EARTH','GM',1);

% Compute orbital period and set time span for the integration
T_orb = 2*pi*(sqrt(a^3/mu_earth)); 
tspan = linspace(0,T_orb,10000);

% Compute orbit's initial state (position [km] and velocity [km/s]) 
% from keplerian parameters
[rr0,vv0] = kep2car(a,e,inc,OM,om,th,mu_earth);
X0 = [rr0; vv0];

% Compute the flow at final time
X_end = flow_2BP(rr0,vv0,0,T_orb,mu_earth,1);
% Propagate the orbit
[~,X] = twobody_propagator(rr0,vv0,mu_earth,tspan);


% Retrieve position and velocity at each integration step
rr = X(:,1:3)';
vv = X(:,4:6)';

% VALIDATE THE PROPAGATOR
% Compute errors w.r.t initial condition to validate the propagator 
% (since tspan = T_orb, we expect the given initial state and the final 
% state computed through the integration to be equal)
err_r = norm(rr0(:,1)-rr(:,end));
err_v = norm(vv0(:,1)-vv(:,end));
err_state = norm(X0'-X(end,:));
err_end = norm(X0-X_end);

disp(' ')
disp('Request 1: 2BP propagator validation:')
fprintf('--> Absolute error between initial and final state: %g\n',err_end)
disp(' ')

% Orbit plot
figure(1)
plot3(rr(1,1),rr(2,1),rr(3,1),'ok')
hold on
plot3(rr(1,:),rr(2,:),rr(3,:),'r','linewidth',1.5)
hold on
PlotEarth
grid on 
box on
xlabel('$r_x$ [km]','Interpreter','latex')
ylabel('$r_y$ [km]','Interpreter','latex')
zlabel('$r_z$ [km]','Interpreter','latex')
legend('Initial state (periapsis)')
title('Ex 1 point 1: propagated orbit plot')


%-------------------------------------------------------------------------%

% REQUEST 2 
disp('Request 2: Implement a Lambert solver')
% In order to test and validate the developed Lambert solver, Mars Express 
% trajectory is used as reference.

% Define departure and arrival data
dep_body = 'EARTH';
arr_body = 'MARS';
dep_epoch_str = '2003-Jun-2 17:45:00 UTC'; 
arr_epoch_str = '2003-Dec-25 03:00:00 UTC';

% Convert string date to ephemeris time (to be used with SPICE)
t_dep = cspice_str2et(dep_epoch_str);
t_arr = cspice_str2et(arr_epoch_str);
% Compute the time of flight [s]
ToF_lambert = t_arr - t_dep;

% Select reference frame (using SPICE naming convention)
ref_frame = 'ECLIPJ2000';

% Define primary body parameters
primary_body = 'SUN';
mu_sun = cspice_bodvrd(primary_body,'GM',1);

% Retrieve initial and final state of the planets:

% Initial state: position [km] and velocity [km/s] of the Earth at departure
State_E = cspice_spkezr(dep_body,t_dep,ref_frame,'NONE',primary_body); 
r_E = State_E(1:3,:);
v_E = State_E(4:6,:);

% Final state: position [km] and velocity [km/s] of Mars at arrival
State_M = cspice_spkezr(arr_body,t_arr,ref_frame,'NONE',primary_body);   
r_M = State_M(1:3,:);
v_M = State_M(4:6,:);

% Solve Lambert problem using the developed solver 
v1_guess = v_E;

tic
lambert.sol = shootingfun_lambert(r_E,r_M,t_dep,t_arr,mu_sun,v1_guess);
time_lambert = toc; 

disp(' ')
fprintf('Initial velocity found using the implemented Lambert solver: [');
fprintf('%f, ',lambert.sol(1:2))
fprintf('%f]'' km/s\n',lambert.sol(end))
fprintf('Elapsed time using the implemented Lambert solver: %f sec \n',time_lambert)

% Solve Lambert problem again using the classic solver
tic
[~,~,~,ERROR,VI,VF,~,~]=lambertMR(r_E,r_M,ToF_lambert,mu_sun,0, 0, 0);
time_lambertMR = toc;

disp(' ')
fprintf('Initial velocity found using the classic Lambert solver: [');
fprintf('%f, ',VI(1:2))
fprintf('%f]'' km/s\n',VI(end))
fprintf('Elapsed time using the classic Lambert solver: %f sec \n',time_lambertMR)

% Compute the error between the two solutions of the Lambert solvers
err_lambert = norm(lambert.sol-VI');

disp(' ')
fprintf('Absolute error between the solutions of the two Lambert solvers: %d\n',err_lambert);

% Propagate and plot the trajectory:
% Compute the semi-major axis [km] of Earth and Mars heliocentric orbits
a_E = -(mu_sun/(2*((1/2*norm(v_E)^2)-(mu_sun/norm(r_E)))));
a_M = -(mu_sun/(2*((1/2*norm(v_M)^2)-(mu_sun/norm(r_M)))));

% Compute Earth and Mars revolution periods [s]
T_E = 2*pi*sqrt(a_E^3/mu_sun);
T_M = 2*pi*sqrt(a_M^3/mu_sun);

% Compute the timespan for the propagation of the Lambert arc and
% Earth and Mars heliocentric orbits
tsp_L = linspace(t_dep,t_arr,10000);
tsp_E = linspace(0,T_E,10000);
tsp_M = linspace(0,T_M,10000);

% Propagate the orbits
[~,L] = twobody_propagator(r_E,lambert.sol,mu_sun,tsp_L);
[~,E] = twobody_propagator(r_E,v_E,mu_sun,tsp_E);
[~,M] = twobody_propagator(r_M,v_M,mu_sun,tsp_M);

% Plot Lambert arc
fig2 = figure(2); 
plot3(E(:,1),E(:,2),E(:,3),'linewidth',1.5)
hold on
plot3(M(:,1),M(:,2),M(:,3),'linewidth',1.5)
hold on
plot3(L(:,1),L(:,2),L(:,3),'linewidth',2)
hold on
sunGlobe = drawBody('Sun', [0,0,0], fig2, 30);
depGlobe = drawBody('Earth', r_E, fig2, 14);
arrGlobe = drawBody('Mars', r_M, fig2, 11);
grid on
xlabel('$r_x$ [km]','Interpreter','latex')
ylabel('$r_y$ [km]','Interpreter','latex')
zlabel('$r_z$ [km]','Interpreter','latex')
legend('Departure orbit (Earth)','Arrival orbit (Mars)',...
    'Lambert arc')
title('Ex 1 point 2: Lambert arc')

%-------------------------------------------------------------------------%

% REQUEST 3

% Set departure time lower boundary
departure_LB = cspice_str2et('2024-Jan-1 00:00:00.0000 TDB');

% Convert Earth-Mars synodic period from days to seconds
synodic_period = 780*cspice_spd();

% Set departure time upper boundary 
departure_UB = departure_LB + synodic_period;

% Compute mean ephemeris time of the departure window 
et_mean = (departure_LB + departure_UB)/2;

% Compute the ToF for a parabolic transfer and an Homann transfer between
% Earth and Mars to define the minimum and maximum ToF
[TPAR, T_hohmann] = E2M_ToF_boundaries(dep_body,arr_body,primary_body,mu_sun,ref_frame,et_mean);

ToF_min = 1.2*TPAR;
ToF_max = 1.6*T_hohmann;

% Set arrival time lower and upper boundaries
arrival_LB = departure_LB + ToF_min;
arrival_UB = departure_UB + ToF_max;

% Set date format for SPICE
t_format = 'DD-Mon-YYYY HR:MN:SC::TDB';

% Compute datestring for arrival and departure UB and LB dates
departure_UB_datestring = cspice_timout(departure_UB,t_format);
arrival_LB_datestring = cspice_timout(arrival_LB,t_format);
arrival_UB_datestring = cspice_timout(arrival_UB,t_format);

% Perform a preliminary grid search to provide fmincon with an educated 
% initial guess 
grid_points = 500; % Number of grid points
departure_window = linspace(departure_LB,departure_UB,grid_points);
arrival_window = linspace(arrival_LB,arrival_UB,grid_points);

% Preallocate variables for computational speed
DV_tot_mat = zeros(length(departure_window),length(arrival_window));
DV_T1 = zeros(length(departure_window),length(arrival_window));
DV_T2 = zeros(length(departure_window),length(arrival_window));

% Grid search
for i=1:length(departure_window)
    
    % Compute Keplerian and cartesian coordinates of departure planet at
    % each step of the departure window
    State_dep = cspice_spkezr(dep_body,departure_window(i),ref_frame,'NONE',primary_body);
    
    for j=1:length(arrival_window)
    % Compute Keplerian and cartesian coordinates of arrival planet at
    % each step of the arrival window
    State_arr = cspice_spkezr(arr_body,arrival_window(j),ref_frame,'NONE',primary_body);
    
    % Compute Time of Flight
    ToF=(arrival_window(j)-departure_window(i)); % ToF must be in [s]
    
    % Check that arrival date is actually after the departure date, if this
    % is true then call lambertMR and get initial and final velocity on the
    % transfer orbit
    if arrival_window(j) < departure_window(i) || ToF > ToF_max
            DV_tot_mat(i,j) = NaN;
            DV_T1(i,j) = NaN;
            DV_T2(i,j)= NaN;
    else
    [~,~,~,~,vvT1,vvT2,~,~] = lambertMR(State_dep(1:3), State_arr(1:3) , ToF, mu_sun, 0, 0, 2);
     
    % Compute DeltaV of the transfer at each time step
    DV_T1(i,j) = norm(vvT1'-State_dep(4:6));
    DV_T2(i,j) = norm(State_arr(4:6)-vvT2');
    DV_tot_mat(i,j) = DV_T1(i,j) + DV_T2(i,j);
    end
    end
    
end

% Find minimum DeltaV tot from the grid search
DVmin_grid = min(min(DV_tot_mat));
% Find the index of the minimum DeltaV tot
[i_min,j_min] = find(DV_tot_mat==DVmin_grid);

% Retrieve optimal arrival and departure ephemeris time associated to the 
% minimum DeltaV obtained from the grid search
opt_dep_grid = departure_window(i_min);
opt_arr_grid = arrival_window(j_min);

% Convert optimal arrival and departure ephemeris time to date string
opt_dep_grid_date = cspice_timout(opt_dep_grid,t_format);
opt_arr_grid_date = cspice_timout(opt_arr_grid,t_format);

disp(' ')
disp('Request 3: Earth to Mars transfer')
disp(' ')
disp('Perform a preliminary grid search over the whole interval:')
fprintf('Min DeltaV obtained from preliminary grid search: %f km/s \n',DVmin_grid)
fprintf('Optimal departure date: %s \n',opt_dep_grid_date)
fprintf('Optimal arrival date: %s \n', opt_arr_grid_date)

% Porkchop plot
% Compute datenum for each date of the departure and arrival windows 
dep_dates = zeros(1,length(departure_window)); % preallocate variable for cpu speed
arr_dates = zeros(1,length(arrival_window)); % preallocate variable for cpu speed
for k = 1:length(departure_window) 
    dep_dates(k) = datenum(cspice_timout(departure_window(k),t_format));
    arr_dates(k) = datenum(cspice_timout(arrival_window(k),t_format));
end
figure(3)
% Define levels and plot contour lines
levels = linspace(floor(DVmin_grid),floor(DVmin_grid+20),40);
PorkChop = contour(dep_dates,arr_dates,DV_tot_mat',levels);
grid on
xlabel('Departure date')
ylabel('Arrival date')
xtickangle(30);
ytickangle(30);
datetick('x','yyyy mmm dd','keeplimits')
datetick('y','yyyy mmm dd','keeplimits')
hold on
% Adjust font size
ha = gca;
ha.FontSize = 12;
% Add colorbar
hbc = colorbar;
hbc.Title.String = '$\Delta v$ [km/s]';
hbc.Title.Interpreter = 'latex';
hbc.Title.FontSize = 15;
% Save the current color range
cRange = caxis; 
% Plot minimum DeltaV location on the grid
DVmin_marker = plot(dep_dates(i_min),arr_dates(j_min),'o','MarkerEdgeColor','k'...
    ,'MarkerFaceColor','r','MarkerSize',7); 
% Plot constant ToF lines (in days)
[C,h] = contour(dep_dates,arr_dates, arr_dates' - dep_dates,80:60:440,'k');
clabel(C,h);
title('porkchop plot')
legend(DVmin_marker,'$\Delta v_{min}$','interpreter','latex','location','southeast');
% Set the color range to the previous one (to avoid colorbar overlapping)
caxis(cRange);


% Define departure and arrival guesses (first guess is the optimal solution
% found by the grid search) for the fmincon optimization
dep_guess.et = [opt_dep_grid, ...
                cspice_str2et('20-May-2024 00:00:00.0000 TDB'), ...
                cspice_str2et('10-Nov-2024 00:00:00.0000 TDB')];
arr_guess.et = [opt_arr_grid, ...  
                cspice_str2et('20-Apr-2025 00:00:00.0000 TDB'),...
                cspice_str2et('10-Jul-2025 00:00:00.0000 TDB')];

% Number of variables of the problem
nvar = 2;
% Seconds per day
spd = cspice_spd;
% Preallocate variable before for loop to improve computational efficiency 
Y0 = zeros(nvar,length(dep_guess.et)); 

for n = 1:length(dep_guess.et)
% Convert guessed arrival and departure ephemeris times to date strings
dep_guess.datestring(n,:) = cspice_timout(dep_guess.et(n),t_format);
arr_guess.datestring(n,:) = cspice_timout(arr_guess.et(n),t_format);
% Compute datenum for each guessed date of departure and arrival
dep_guess.datenum(n) = datenum(dep_guess.datestring(n,:));
arr_guess.datenum(n) = datenum(arr_guess.datestring(n,:));
% Define fmincon initial guess for each guessed date of departure and
% arrival
Y0(:,n) = [dep_guess.et(:,n)/spd; arr_guess.et(:,n)/spd];
end

% DEFINE SHOOTING PROBLEM STRUCTURE
% Define objective function
E2M.objective = @(Y) E2M_objectivefun(Y,dep_body,arr_body,primary_body,mu_sun,ref_frame);

% Set lower and upper boundaries
E2M.lb = [departure_LB/spd; arrival_LB/spd];
E2M.ub = [departure_UB/spd; arrival_UB/spd];
E2M.A =[1 -1];
E2M.B = 0;

% Set solver options
E2M.solver='fmincon';
options = optimoptions('fmincon');
options.ConstraintTolerance = 1;
options.Display = 'iter';
options.PlotFcn = 'optimplotfval';
E2M.options = options;

% Set initial guess
E2M.x0 = Y0;

% Preallocate variables before for loop to improve computational efficiency 
t = zeros(nvar,length(Y0));
DV_fmincon = zeros(1,length(Y0));

% Solve the shooting problem for each guessed departure and arrival date
disp('Solve the shooting problem using fmincon for different initial guesses')
for l = 1:length(Y0)   
E2M.x0 = Y0(:,l);
[t(:,l),DV_fmincon(:,l)] = fmincon(E2M);
end

% Convert the solution (departure and arrival ephemeris time) from days to sec
t = t.*spd;

% Reorganize the solution of the fmincon into a structure containing 
% departure and arrival dates  
dep_sol.et = t(1,:);
arr_sol.et = t(2,:);

% Preallocate variables before for loop to improve computational efficiency
vv_T1 = zeros(3,length(dep_sol.et));
vv_T2 = zeros(3,length(dep_sol.et));
t_span = zeros(10000,length(dep_sol.et));

for n = 1:length(dep_sol.et)
    
% Convert obtained arrival and departure ephemeris times to date strings
dep_sol.datestring(n,:) = cspice_timout(dep_sol.et(n),t_format);
arr_sol.datestring(n,:) = cspice_timout(arr_sol.et(n),t_format);

% Compute datenum for each obtained date of departure and arrival 
dep_sol.datenum(n) = datenum(dep_sol.datestring(n,:));
arr_sol.datenum(n) = datenum(arr_sol.datestring(n,:));

% Compute the state of Earth and Mars at departure and arrival respectively
State.dep(:,n) = cspice_spkezr(dep_body,dep_sol.et(n),ref_frame,'NONE',primary_body);
State.arr(:,n) = cspice_spkezr(arr_body,arr_sol.et(n),ref_frame,'NONE',primary_body);

% Compute ToF associated to each solution found 
ToF(:,n) = arr_sol.et(n) - dep_sol.et(n);

% Characterize and propagate transfer trajectory associated to each 
% solution found
[~,~,~,~,vv_T1(:,n),vv_T2(:,n),~,~] = lambertMR(State.dep(1:3,n),State.arr(1:3,n),ToF(:,n),mu_sun,0,0,2);
t_span(:,n) = linspace(dep_sol.et(n),arr_sol.et(n),10000);
[~,State.transf_traj(:,:,n)] = twobody_propagator(State.dep(1:3,n),vv_T1(:,n),mu_sun,t_span(:,n)');

end

disp('First fmincon initial guess:')
fprintf('Guessed departure date: %s\n',dep_guess.datestring(1,:))
fprintf('Guessed arrival date: %s\n',arr_guess.datestring(1,:))
disp('Solution associated to the first initial guess:')
fprintf('Optimal departure date: %s\n',dep_sol.datestring(1,:))
fprintf('Optimal arrival date: %s\n',arr_sol.datestring(1,:))
fprintf('Time of flight: %d days\n',round(ToF(1)/spd))
fprintf('DeltaV associated to the first initial guess: %f km/s\n',DV_fmincon(1))
disp(' ')
disp('Second fmincon initial guess:')
fprintf('Guessed departure date: %s\n',dep_guess.datestring(2,:))
fprintf('Guessed arrival date: %s\n',arr_guess.datestring(2,:))
disp('Solution associated to the second initial guess:')
fprintf('Optimal departure date: %s\n',dep_sol.datestring(2,:))
fprintf('Optimal arrival date: %s\n',arr_sol.datestring(2,:))
fprintf('Time of flight: %d days\n',round(ToF(2)/spd))
fprintf('DeltaV associated to the second initial guess: %f km/s\n',DV_fmincon(2))
disp(' ')
disp('Third fmincon initial guess:')
fprintf('Guessed departure date: %s\n',dep_guess.datestring(3,:))
fprintf('Guessed arrival date: %s\n',arr_guess.datestring(3,:))
disp('Solution associated to the third initial guess:')
fprintf('Optimal departure date: %s\n',dep_sol.datestring(3,:))
fprintf('Optimal arrival date: %s\n',arr_sol.datestring(3,:))
fprintf('Time of flight: %d days\n',round(ToF(3)/spd))
fprintf('DeltaV associated to the third initial guess: %f km/s\n',DV_fmincon(3))

% Plot initial guesses and corresponding fmincon solutions
figure(5)
PorkChop2 = contour(dep_dates,arr_dates,DV_tot_mat',levels);
grid on
xlabel('Departure date')
ylabel('Arrival date')
xtickangle(30);
ytickangle(30);
datetick('x','yyyy mmm dd','keeplimits')
datetick('y','yyyy mmm dd','keeplimits')
hold on
% Adjust font size
ha = gca;
ha.FontSize = 12;
% Add colorbar
hbc = colorbar;
hbc.Title.String = '$\Delta v$ [km/s]';
hbc.Title.Interpreter = 'latex';
hbc.Title.FontSize = 15;
% Save the current color range
cRange = caxis; 
% Plot constant ToF lines (in days)
[C,h] = contour(dep_dates,arr_dates, arr_dates' - dep_dates,80:60:440,'k');
clabel(C,h);
% Set the color range to the previous one (to avoid colorbar overlapping)
caxis(cRange);
% Plot guesses and corresponding solutions
second_guess_marker = plot(dep_guess.datenum(2),arr_guess.datenum(2),...
    'o','MarkerEdgeColor','k','MarkerFaceColor','k','MarkerSize',7);
DV_secondsol_marker = plot(dep_sol.datenum(2),arr_sol.datenum(2),...
    'o','MarkerEdgeColor','k','MarkerFaceColor','r','MarkerSize',7); 
third_guess_marker = plot(dep_guess.datenum(3),arr_guess.datenum(3),...
    'o','MarkerEdgeColor','k','MarkerFaceColor','g','MarkerSize',7);
DV_thirdsol_marker = plot(dep_sol.datenum(3),arr_sol.datenum(3),...
    'o','MarkerEdgeColor','k','MarkerFaceColor','b','MarkerSize',7); 
legend([second_guess_marker DV_secondsol_marker third_guess_marker DV_thirdsol_marker],...
    'Second initial guess', 'Second solution','Third initial guess',...
    'Third solution','location','southeast');

% Plot transfer trajectories found by fmincon
titles = {'Ex 1 point 3: First solution trajectory',...
          'Ex 1 point 3: Second solution trajectory',...
          'Ex 1 point 3: Third solution trajectory'};
for p = 1:length(dep_sol.et)
figure(p+5); 
plot3(E(:,1),E(:,2),E(:,3),'linewidth',1.5)
hold on
plot3(M(:,1),M(:,2),M(:,3),'linewidth',1.5)
hold on
plot3(State.transf_traj(:,1,p),State.transf_traj(:,2,p),State.transf_traj(:,3,p),'linewidth',2)
hold on
sunGlobe = drawBody('Sun', [0,0,0], figure(p+5), 30);
depGlobe = drawBody('Earth', State.dep(1:3,p), figure(p+5), 14);
arrGlobe = drawBody('Mars', State.arr(1:3,p), figure(p+5), 11);
grid on
xlabel('$r_x$ [km]','Interpreter','latex')
ylabel('$r_y$ [km]','Interpreter','latex')
zlabel('$r_z$ [km]','Interpreter','latex')
legend('Departure orbit (Earth)','Arrival orbit (Mars)',...
    'Transfer trajectory') 
title(titles{p})
end

% Clear kernel pool
cspice_kclear

% Check that kernel pool is clear
fprintf('\nTotal kernels number after kclear: %d\n', cspice_ktotal('ALL'));

%-------------------------------------------------------------------------%
%% EXERCISE 2

% Clear memory workspace and set path
clearvars; close all; clc
addpath(genpath(pwd)); 

disp('Exercise 2')

% Load SPICE kernels: 
% de425s.bsp, gm_de431.tpc, naif0012.tls, pck00010.tpc
cspice_furnsh('assignment01.tm');

% Count total kernels number
fprintf('Total kernels number: %d\n', cspice_ktotal('ALL'));

% Define cell-array with object labels
labels = {'Earth'; 'Sun'; 'Moon'}; 

% Initialize planetary data for n-body propagation
bodies = nbody_init(labels);

% Earth standard gravitational parameter [km^3/s^2]
mu_earth = bodies{1}.GM;

% Set reference frame (using SPICE naming convention)
ref_frame = 'J2000';

% Define departure date guess
dep_epoch_str = '2021-Dec-24 00:00:00.0000 TDB';
t1 = cspice_str2et(dep_epoch_str);

% Set date format for SPICE conversion
t_format = 'DD-Mon-YYYY HR:MN:SC::TDB';

% Seconds per day
spd = cspice_spd;

% Initial guess generation (bielliptic 2BP transfer starting from a 200 km 
% altitude equatorial LEO)
% Guess initial orbit: perigee radius 6578.1 [km] apogee radius 750000 [km]
% Set perigee and apogee radius [km]
rr1 = [6578.1;0;0];
rr2 = [-750000;0;0];
% Compute first orbit semi-major axis [km], eccentricity, semi-latus rectum
% [km], specific angular momentum [km^2/s] and orbital period [s]
a1 = (norm(rr1)+norm(rr2))/2;
e1 = (norm(rr2)-norm(rr1))/(norm(rr2)+norm(rr1));
p1 = a1*(1-e1^2);
h1 = sqrt(mu_earth*p1);
T_orb1 = 2*pi*sqrt(a1^3/mu_earth);

% Compute 1st arc reference ToF
tof1_ideal = T_orb1/2;

% Compute velocity at pericenter of the first orbit
v_th_p1 = (mu_earth/h1)*(1+e1*cos(0));
vv1 = [0; v_th_p1; 0];

% Define second orbit
r_moon = 384400;
rp2 = r_moon + 60000;
ra2 = norm(rr2);

% Compute second orbit semi-major axis [km], eccentricity, semi-latus rectum
% [km], specific angular momentum [km^2/s] and orbital period [s]
a2 = (ra2+rp2)/2;
e2 = (ra2-rp2)/(ra2+rp2);
p2 = a2*(1-e2^2);
h2 = sqrt(mu_earth*p2);
T_orb2 = 2*pi*sqrt(a2^3/mu_earth);

% Compute 2nd arc reference ToF
tof2_ideal = T_orb2/2;

% Compute velocity at apocenter of the second orbit
v_th_a2 = (mu_earth/h2)*(1+e2*cos(pi));
vv2 = [0; v_th_a2; 0];

% Astronomical Unit in [km]
AU = 149597870.691;

% Initial state guess
X1i = [rr1/AU; vv1];
X2i = [rr2/AU; vv2];

% Define second burn time and arrival time guesses
t2 = tof1_ideal + t1;
t3 = tof2_ideal + t2;

% Initial guess array
Y0 = [X1i; X2i; t1/spd; t2/spd; t3/spd]; 

% DEFINE MULTIPLE SHOOTING PROBLEM STRUCTURE: E2EML2
% Define objective function
E2EML2.objective = @(y) E2EML2_DeltaV(y,bodies,ref_frame);

% Set initial guess
E2EML2.x0 = Y0;

% Set problem's constraints 
E2EML2.A = [zeros(1,12) 1 0 -1;
            zeros(1,12) 1 -1 0;
            zeros(1,12) 0 1 -1];
E2EML2.B = [0;0;0];
E2EML2.Aeq = [0 0 1 zeros(1,12)];
E2EML2.Beq = 0;
E2EML2.nonlcon = @(y) E2EML2_constraints(y,bodies,ref_frame);

% Set solver options
E2EML2.solver = 'fmincon';
options = optimoptions('fmincon');
options.Algorithm = 'active-set';
options.ConstraintTolerance = 1;
options.Display = 'iter';
options.PlotFcn = 'optimplotfval';
options.MaxFunctionEvaluations = 20000;
options.MaxIterations = 100;
E2EML2.options = options;

% Find optimal solution
[Y,DeltaV,~,output,lambda,gradient,H] = fmincon(E2EML2);

% Propagate and plot optimal solution
ode_options = odeset( 'RelTol', 1e-12, 'AbsTol', 1e-12);
[timevec1,State1] = ode113(@(t,x) nbody_rhs(t,x,bodies,ref_frame),[Y(13)*spd Y(14)*spd],[Y(1:3)*AU;Y(4:6)],ode_options);
[timevec2,State2] = ode113(@(t,x) nbody_rhs(t,x,bodies,ref_frame),[Y(14)*spd Y(15)*spd],[Y(7:9)*AU;Y(10:12)],ode_options);

% Retrieve trajectory and time
YY = [State1;State2];
timevec = [timevec1;timevec2];

% Plot 
figure(2)
plot3(State1(:,1),State1(:,2),State1(:,3),'linewidth',1)
hold on
plot3(State2(:,1),State2(:,2),State2(:,3),'linewidth',1)
hold on
r_moon = cspice_spkpos('MOON',timevec(end),ref_frame,'NONE','EARTH');
r_L2 = r_moon/norm(r_moon)*(norm(r_moon)+60e3);
plot3(r_L2(1),r_L2(2),r_L2(3),'or')
PlotEarth
hold on
grid on
axis equal
legend('First arc', 'Second arc','EML2 at final time')
xlabel('$r_x$ [km]','Interpreter','latex')
ylabel('$r_y$ [km]','Interpreter','latex')
zlabel('$r_z$ [km]','Interpreter','latex')

% Compute ToF (in days)
TOF1 = Y(14)-Y(13);
TOF2 = Y(15)-Y(14);

% Convert optimal arrival and departure ephemeris time to date string
optimal_departure = cspice_timout(Y(13)*spd,t_format);
optimal_arrival = cspice_timout(Y(15)*spd,t_format);

% Evaluate DeltaVs [km/s]
[DeltaV_tot,DeltaV_1,DeltaV_2,DeltaV_3] = E2EML2_DeltaV(Y,bodies,ref_frame);

% Evaluate constraints
[c, ceq] = E2EML2_constraints(Y,bodies,ref_frame);

% Display results of second request
disp(' ')
disp('fmincon optimal solution:')
fprintf('Guessed departure date: %s\n',dep_epoch_str)
fprintf('Optimal departure date: %s\n', optimal_departure);
fprintf('Optimal arrival date: %s\n', optimal_arrival);
fprintf('First arc ToF: %f days\n',TOF1)
fprintf('Second arc ToF: %f days\n',TOF2)
fprintf('Total ToF: %f days\n',TOF1+TOF2)
fprintf('DeltaV_1: %f km/s\n',DeltaV_1)
fprintf('DeltaV_2: %f km/s\n',DeltaV_2)
fprintf('DeltaV_3: %f km/s\n',DeltaV_3)
fprintf('DeltaV_tot: %f km/s\n',DeltaV_tot)
disp(' ')

%-------------------------------------------------------------------------%

% REQUEST 3
disp('Request 3: evaluate DeltaVtot as function of departure epoch')
disp(' ')
disp('Wait for computation')

% Define departure epochs spanning an entire year (2 points per month
% discretization
dep_epoch_str_1year_guess = '2022-Jan-01 00:00:00.0000 TDB';
arr_epoch_str_1year_guess = '2023-Jan-01 00:00:00.0000 TDB';
t_dep_1year_guess_in = cspice_str2et(dep_epoch_str_1year_guess);
t_dep_1year_guess_fin = cspice_str2et(arr_epoch_str_1year_guess);
% Departure guesses (2 per month, along 1 year timespan)
t1_1year_guess.vect = linspace(t_dep_1year_guess_in,t_dep_1year_guess_fin,26);

% Preallocate variables before for loop to improve computational speed
Y0_1year = zeros(15,length(t1_1year_guess.vect));

% Compute departure guesses, second burn time guesses, arrival date guesses
% and convert them from ephemeris time to date string and datenum 
for n = 1:length(t1_1year_guess.vect)
    
    t1_1year_guess.datestring(n,:) = cspice_timout(t1_1year_guess.vect(n),t_format);
    t1_1year_guess.datenum(n) = datenum(t1_1year_guess.datestring(n,:));
    
    % Adjust initial guesses for second burn and arrival times for each date 
    if n==3 || n==7 || n==12 || n==17 || n==20 || n==24 || n==26
        t2_1year_guess.vect(n) = floor(tof1_ideal/spd)*spd + t1_1year_guess.vect(n);
        t3_1year_guess.vect(n) = floor(tof2_ideal/spd)*spd + t2_1year_guess.vect(n);
    elseif n==21 
        t2_1year_guess.vect(n) = ceil(tof1_ideal/spd)*spd + t1_1year_guess.vect(n);
        t3_1year_guess.vect(n) = floor(tof2_ideal/spd)*spd + t2_1year_guess.vect(n);
    elseif n==5 || n==23
        t2_1year_guess.vect(n) = ceil(tof1_ideal/spd)*spd + t1_1year_guess.vect(n);
        t3_1year_guess.vect(n) = ceil(tof2_ideal/spd)*spd + t2_1year_guess.vect(n);
    elseif n==6 || n==19
        t2_1year_guess.vect(n) = floor((tof1_ideal/spd)-1)*spd + t1_1year_guess.vect(n);
        t3_1year_guess.vect(n) = floor(tof2_ideal/spd)*spd + t2_1year_guess.vect(n);
    elseif n==8 || n==10 || n==18
        t2_1year_guess.vect(n) = floor((tof1_ideal/spd)+2)*spd + t1_1year_guess.vect(n);
        t3_1year_guess.vect(n) = floor((tof2_ideal/spd)+5)*spd + t2_1year_guess.vect(n);
    else 
        t2_1year_guess.vect(n) = 13.399*spd + t1_1year_guess.vect(n); %tof1ideal+t1
        t3_1year_guess.vect(n) = 26.579*spd + t2_1year_guess.vect(n); %tof2ideal+t2
    end
    
    % Convert ephemeris time to date string and datenum
    t2_1year_guess.datestring(n,:) = cspice_timout(t2_1year_guess.vect(n),t_format);
    t2_1year_guess.datenum(n) = datenum(t2_1year_guess.datestring(n,:));
    
    t3_1year_guess.datestring(n,:) = cspice_timout(t3_1year_guess.vect(n),t_format);
    t3_1year_guess.datenum(n) = datenum(t3_1year_guess.datestring(n,:));
    
    % Define initial guesses matrix
    Y0_1year(:,n) = [rr1/AU; vv1; rr2/AU; vv2; t1_1year_guess.vect(n)/spd;...
                        t2_1year_guess.vect(n)/spd; t3_1year_guess.vect(n)/spd];
end

% Preallocate variables before for loop to improve computational speed
Y_1year = zeros(15,length(t1_1year_guess.vect));
DeltaV_1year = zeros(1,length(t1_1year_guess.vect));

% Solve the multiple shooting problem for several departure guesses
for l = 1:length(t1_1year_guess.vect)   
    
fprintf('Iteration %d\n',l)
E2EML2.x0 = Y0_1year(:,l);
[Y_1year(:,l),DeltaV_1year(:,l)] = fmincon(E2EML2);

% Retrieve the solution for the optimal departure date, optimal second burn
% time and optimal arrival date for each guess
t1_1year.vect(l) = Y_1year(13,l)*spd;
t2_1year.vect(l) = Y_1year(14,l)*spd;
t3_1year.vect(l) = Y_1year(15,l)*spd;

% Convert ephemeris time to date string and datenum for better results presentation
t1_1year.datestring(l,:) = cspice_timout(t1_1year.vect(l),t_format);
t1_1year.datenum(l) = datenum(t1_1year.datestring(l,:));
t2_1year.datestring(l,:) = cspice_timout(t2_1year.vect(l),t_format);
t2_1year.datenum(l) = datenum(t2_1year.datestring(l,:));
t3_1year.datestring(l,:) = cspice_timout(t3_1year.vect(l),t_format);
t3_1year.datenum(l) = datenum(t3_1year.datestring(l,:));

% Retrieve the ToF (in days) from each solution
TOF1_1year.vect(l) = Y_1year(14,l) - Y_1year(13,l);
TOF2_1year.vect(l) = Y_1year(15,l) - Y_1year(14,l);

end

% Plot DeltaVtot as function of departure epoch
figure(4)
plot(t1_1year.datenum,DeltaV_1year,'r')
hold on
plot(t1_1year.datenum,DeltaV_1year,'ob','MarkerSize',7,'linewidth',1.5)
datetick('x','yyyy mmm dd','keeplimits')
xlabel('Departure epoch','fontsize',13)
ylabel('\Deltav_{tot} [km/s]','fontsize',13)
ylim([3,6])
grid on

% Clear kernel pool
cspice_kclear

% Check that kernel pool is clear
fprintf('\nTotal kernels number after kclear: %d\n', cspice_ktotal('ALL'));

%-------------------------------------------------------------------------%
%% EXERCISE 3

% Clear memory workspace and set path
clearvars; close all; clc
addpath(genpath(pwd)); 

disp('Exercise 3')

% Load SPICE kernels: 
% de425s.bsp, gm_de431.tpc, naif0012.tls, pck00010.tpc
cspice_furnsh('assignment01.tm');

% Count total kernels number
fprintf('Total kernels number: %d\n', cspice_ktotal('ALL'));

% Define initial and final state (position [km] and velocity [km/s])
rr1_0 = [0; -29597.43; 0];
vv1_0 = [1.8349; 0.0002; 3.1783];

rr2_f = [0; -29617.43; 0];
vv2_f = [1.8371; 0.0002; 3.1755];

% S/C mass [kg]
m0 = 735;
% Maximum thrust [kN]
Tmax = 100*10^-6;  
% Specific impulse [s]
Isp = 3000;
% Standard acceleration of free fall [km/s^2]
g0 = 9.80665*10^-3;
% Earth standard gravitational parameter [km^3/s^2]
mu = cspice_bodvrd('Earth','GM',1);

% Compute semi-major axis [km] of initial and final orbit
a_1 = -(mu/(2*((1/2*norm(vv1_0)^2)-(mu/norm(rr1_0)))));
a_2 = -(mu/(2*((1/2*norm(vv2_f)^2)-(mu/norm(rr2_f)))));
% Compute orbital period [s] of initial and final orbit

T_orb1 = 2*pi*sqrt(a_1^3/mu);
T_orb2 = 2*pi*sqrt(a_2^3/mu);

% Set initial time instant and define time span for the orbit propagation
t0 = 0;
tspan_1 = linspace(t0,T_orb1,10000);
tspan_2 = linspace(t0,T_orb2,10000);

% Propagate the initial and final orbits (in the geocentric 2BP)
[~,X1] = twobody_propagator(rr1_0,vv1_0,mu,tspan_1);
[~,X2] = twobody_propagator(rr2_f,vv2_f,mu,tspan_2);

% Seconds per day
spd = cspice_spd;

% Initial state of the s/c
x0 = [rr1_0; vv1_0; m0];

% Initial unknown guess 
z0_guess = [5; 5; 5; 1e4; 1e4; 1e4; 1; T_orb1/spd];

% Set fsolve options
options = optimoptions('fsolve');
options.Display = 'iter';
options.MaxFunctionEvaluations = 2000;
options.MaxIterations = 500;

% Solve the TPBVP 
[SOL,fval] = fsolve(@(z) shootingfun_TPBVP(z,x0,t0,rr2_f,vv2_f,mu,Tmax,Isp,g0),z0_guess,options);

% Retrieve costate and final time from the solution of the TPBVP
lambda_r0 = SOL(1:3);
lambda_v0 = SOL(4:6);
lambda_m0 = SOL(7);
tf = SOL(8)*spd; % [s]

% Convert transfer time from [s] to [hh:mm:ss] 
transfer_time = datestr(seconds(tf),'HH:MM:SS');

% Display results
disp(' ')
fprintf('Solution of the TPBVP is: [');
fprintf('%f, ',SOL(1:7))
fprintf('%f]''\n',SOL(end))
fprintf('Minimum transfer time is: %s hrs',transfer_time(:,1:2))
fprintf(' %s min',transfer_time(:,4:5))
fprintf(' %s sec\n',transfer_time(:,7:8))

% Compute the flow of the problem starting from the solution found
[constraints,time,dynamics] = shootingfun_TPBVP(SOL,x0,t0,rr2_f,vv2_f,mu,Tmax,Isp,g0);

% Retrieve the state and the costate flow 
r_flow = dynamics(:,1:3);
v_flow = dynamics(:,4:6);
m_flow = dynamics(:,7);
lambda_r_flow = dynamics(:,8:10);
lambda_v_flow = dynamics(:,11:13);
lambda_m_flow = dynamics(:,14);

% Plot transfer orbit
figure(1)
plot3(X1(1,1),X1(1,2),X1(1,3),'o','MarkerFaceColor','g','MarkerEdgeColor','k','MarkerSize',8)
hold on
plot3(X2(1,1),X2(1,2),X2(1,3),'o','MarkerFaceColor','y','MarkerEdgeColor','k','MarkerSize',8)
hold on 
plot3(X1(:,1),X1(:,2),X1(:,3),'r','linewidth',1)
hold on
plot3(X2(:,1),X2(:,2),X2(:,3),'b','linewidth',1)
hold on
plot3(r_flow(:,1),r_flow(:,2),r_flow(:,3),'c','linewidth',1)
hold on
PlotEarth
axis equal
grid on
legend('Initial point','Target point','Initial orbit','Target orbit',...
       'Thrust arc')
xlabel('$r_x$ [km]','Interpreter','latex')
ylabel('$r_y$ [km]','Interpreter','latex')
zlabel('$r_z$ [km]','Interpreter','latex')

figure(2)
plot3(r_flow(:,1),r_flow(:,2),r_flow(:,3),'linewidth',1.3)
hold on
plot3(X1(1,1),X1(1,2),X1(1,3),'o','MarkerFaceColor','g','MarkerEdgeColor','k','MarkerSize',8)
hold on
plot3(X2(1,1),X2(1,2),X2(1,3),'o','MarkerFaceColor','y','MarkerEdgeColor','k','MarkerSize',8)
hold on
PlotEarth
grid on
axis equal 
legend('Thrust arc','Initial point','Target point')
xlabel('$r_x$ [km]','Interpreter','latex')
ylabel('$r_y$ [km]','Interpreter','latex')
zlabel('$r_z$ [km]','Interpreter','latex')

% Compute thrust throttle factor (u) and switching function (St) for the
% solution found
u = zeros(length(time),1);
St = zeros(length(time),1);

for k = 1:length(time)
[u(k,:),St(k,:)] = thrust_switch(lambda_v_flow(k,:),lambda_m_flow(k,:),Isp,m_flow(k,:),g0);
end

% Plot lambda_m, switching function and thrust throttle factor
figure(3)
% lambda_m
subplot(3,1,1) 
plot(time,lambda_m_flow,'linewidth',1.5)
xlim([0,time(end)])
ylim([0,0.12])
xlabel('Time [s]')
ylabel('\lambda_m')
grid on
% St
subplot(3,1,2) 
plot(time,St,'linewidth',1.5)
xlim([0,time(end)])
xlabel('Time [s]')
ylabel('S_t')
grid on
% u
subplot(3,1,3) 
plot(time,u,'linewidth',1.5)
xlim([0,time(end)])
xlabel('Time [s]')
ylabel('u')
grid on

%-------------------------------------------------------------------------%

%% OPTIONAL REQUEST: Compute tf for several values of Tmax
disp(' ')
disp('Optional request:')
disp('Transfer time as function of the maximum thrust value:')

% Define an array of Tmax
Tmax_vect = (100:100:1000)*10^-6;

% Solve the TPBVP for each value of Tmax finding the transfer time tf
tf_vect = zeros(1,length(Tmax_vect));
for j = 1:length(Tmax_vect)
    
    if Tmax_vect(j) < 400*10^-6
        z0_guess = [5; 5; 5; 1e4; 1e4; 1e4; 1; T_orb1/spd];
    else % use a different lambda0_guess to grant convergence
        z0_guess = [0; 5; 0; 4e3; 4e3; 4e3; 0; T_orb1/spd];
    end    
    
    [SOL(:,j),fval(:,j)] = fsolve(@(z) shootingfun_TPBVP(z,x0,t0,rr2_f,vv2_f,mu,Tmax_vect(j),Isp,g0),z0_guess,options);
    tf_vect(:,j) = SOL(8,j)*spd;
    
end

% Plot tf as function of Tmax
figure(5)
plot(Tmax_vect/10^-6,tf_vect/3600,'r','linewidth',1)
hold on
plot(Tmax_vect/10^-6,tf_vect/3600,'ob')
xlim([100,1000])
xlabel('Thrust [mN]')
ylabel('Transfer time [h]')
% legend('Transfer time','Thrust value')
grid on

% Display results of the optional requests
tf_vect_string = datestr(seconds(tf_vect),'HH:MM:SS');
for n = 1:length(Tmax_vect)
    fprintf('Thrust: %g mN',Tmax_vect(n)/10^-6)
    fprintf(' --> Transfer time: %s hrs',tf_vect_string(n,1:2))
    fprintf(' %s min',tf_vect_string(n,4:5))
    fprintf(' %s sec\n',tf_vect_string(n,7:8))
end

% Clear kernel pool
cspice_kclear

% Check that kernel pool is clear
fprintf('\nTotal kernels number after kclear: %d\n', cspice_ktotal('ALL'));

%-------------------------------------------------------------------------%

%-------------------------------------------------------------------------%
%% FUNCTIONS

% function: kep2car.m

function [rr, vv] = kep2car(a,e,i,OM,om,th,mu)
%-------------------------------------------------------------------------%
%
% kep2car.m trasforms Keplerian parameters to cartesian coordinates.
%
%-------------------------------------------------------------------------%
% PROTOTYPE:
%  [rr, vv] = kep2car(a, e, i, OM, om, th, mu)
% 
%-------------------------------------------------------------------------%
% INPUT ARGUMENTS:
%  a             [1]     Semi-major axis                      [km]
%  e             [1]     Eccentricity                         [-]
%  i             [1]     Inclination                          [rad]
%  OM            [1]     RAAN                                 [rad]
%  om            [1]     Argument of periapsis                [rad]
%  th            [1]     True anomaly                         [rad]
%  mu            [1]     Standard gravitational parameter     [km^3/s^2]
%                        of the primary body
%
%-------------------------------------------------------------------------%
% OUTPUT ARGUMENTS:
%  rr            [3x1]   Position vector                      [km]
%  vv            [3x1]   Velocity vector                      [km/s]
%
%-------------------------------------------------------------------------%
% AUTHOR:
%  Gian Marco Paldino
%
%-------------------------------------------------------------------------%
% VERSIONS: 
%  23-10-2021: First version
%
%-------------------------------------------------------------------------%

    % Semi-latus rectum
    p = a*(1-e^2); 
    
    % Position vector norm
    r = p/(1+e*cos(th));
        
    % Position and velocity vectors in Perifocal Coordinates [3x1]
    rr_pf = r.*[cos(th); sin(th); 0];
    
    vv_pf = (sqrt(mu/p))*[-sin(th); e+cos(th); 0];
    
    % First rotation matrix (around k versor)
    R_OM = [cos(OM) , sin(OM) , 0;...
           -sin(OM) , cos(OM) , 0;...
           0        , 0       , 1];

    % Second rotation matrix (around i' versor)
    R_i  = [1      , 0       ,     0 ;...
            0      , cos(i)  , sin(i);...
            0      , -sin(i) , cos(i)];   
        
    % Third rotation matrix (around k'' versor)
    R_om = [cos(om) , sin(om) ,     0 ;...
            -sin(om), cos(om) ,     0 ;...
            0       , 0       ,     1 ];
    
    % Rotation matrix from ECI to PF
    T_eci2pf = R_om * R_i * R_OM;
    
    % Rotation matrix from PF to ECI
    T_pf2eci = T_eci2pf.';
   
    % Position and velocity vectors rotated in ECI 
    rr = T_pf2eci*rr_pf;
    vv = T_pf2eci*vv_pf;
       
end
   
%-------------------------------------------------------------------------%    

%-------------------------------------------------------------------------%

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

% Normalize position vector
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

% function: twobody_propagator.m

function [T,X] = twobody_propagator(rr0,vv0,mu,tspan)
%-------------------------------------------------------------------------%
%
% orbit_propagation.m performs the numerical integration of the equations 
% of motion for the 2BP. 
%
%-------------------------------------------------------------------------%
% PROTOTYPE;
%  [T,X] = twobody_propagator(rr0,vv0,mu,tspan)
%
%-------------------------------------------------------------------------%
% INPUT ARGUMENTS:
%  rr0              [3x1]  Initial position                  [km]
%  vv0              [3x1]  Initial velocity                  [km/s]
%  mu               [1]    Standard gravitational parameter  [km^3/s^2]
%                          of the primary body
%  t_span           [1xn]  Integration time steps vector     [s]  
%
%-------------------------------------------------------------------------%
% OUTPUT ARGUMENTS: 
%  T               [1xn]   Integrated orbit time steps       [s]
%  X               [6x1]   State of the system:              [km, km/s]
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
ode_options = odeset( 'RelTol', 1e-13, 'AbsTol', 1e-14 );

% Perform the Integration 
[T, X] = ode113( @(t,x) odetwobody(t,x,mu), tspan, X0, ode_options);

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
flow_options = odeset( 'RelTol', 1e-13, 'AbsTol', 1e-14 );

% Set initial conditions for the integration of the ODE  
X0 = [rr0; vv0];

% Perform the integration
[~, X] = ode113( @(t,x) odetwobody(t,x,mu), tspan, X0, flow_options );

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

% function: shootingfun_lambert.m

function v1 = shootingfun_lambert(r1,r2,ti,tf,mu,v1_guess)
% -------------------------------------------------------------------------%
% 
% shootingfun_lambert.m provides the initial velocity solving the Lambert's
% problem using a finite differences approach.
% 
% -------------------------------------------------------------------------%
% PROTOTYPE:
%  Delta_r = shootingfun_lambert(vi, ti, tf,mu,r1,r2)
%
%-------------------------------------------------------------------------%
% INPUT ARGUMENTS:
%  vi              [3x1]  Initial velocity                   [km/s]
%  ti              [1]    Initial time                       [s]
%  tf              [1]    Final time                         [s]
%  mu              [1]    Standard gravitational parameter   [km^3/s^2]
%                         of the primary body
%  r1              [3x1]  Initial position                   [km]
%  r2              [3x1]  Final position                     [km]
%                   
%-------------------------------------------------------------------------%
% OUTPUT ARGUMENTS: 
%  Delta_r         [3x1] Difference between final position   [km]
%
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

% Retrieve velocity components from the initial guess
v_x = v1_guess(1);
v_y = v1_guess(2);
v_z = v1_guess(3);

% Set while loop parameters (iterations, error, and tolerance)
iter = 1;
iter_max = 30;
err = 1; 
tol = 1e-12;

% Compute while loop
while iter < iter_max && err > tol

% Define initial condition
x0 = [r1; v1_guess];

% Compute the flow of the 2BP starting from the initial guessed velocity
flow = flow_2BP(r1,v1_guess,ti,tf,mu,1);

% Compute perturbations along vx vy vz
eps_vx0 = [0; 0; 0; sqrt(eps)*max(1,abs(v_x)); 0; 0];
eps_vy0 = [0; 0; 0; 0; sqrt(eps)*max(1,abs(v_y)); 0];
eps_vz0 = [0; 0; 0; 0; 0; sqrt(eps)*max(1,abs(v_z))];

% Define perturbed initial conditions
X0_pert_vx = x0 + eps_vx0;
X0_pert_vy = x0 + eps_vy0;
X0_pert_vz = x0 + eps_vz0;


% Compute the flow of the 2BP starting from the perturbed initial
% conditions
phi_eps_vx0 = flow_2BP(X0_pert_vx(1:3),X0_pert_vx(4:6),ti,tf,mu,1);
phi_eps_vy0 = flow_2BP(X0_pert_vy(1:3),X0_pert_vy(4:6),ti,tf,mu,1);
phi_eps_vz0 = flow_2BP(X0_pert_vz(1:3),X0_pert_vz(4:6),ti,tf,mu,1);

% Initialize the Jacobian 
J=zeros(3,3);

% Compute the Jacobian components
for k = 1:3
    J(k,1) = (phi_eps_vx0(k)-flow(k))/(sqrt(eps)*max(1,abs(v_x)));
    J(k,2) = (phi_eps_vy0(k)-flow(k))/(sqrt(eps)*max(1,abs(v_y)));
    J(k,3) = (phi_eps_vz0(k)-flow(k))/(sqrt(eps)*max(1,abs(v_z)));
end

% Compute delta_r and delta_v
delta_r2 = flow(1:3) - r2;
delta_v1 = (J^-1) * delta_r2;

% Retrieve initial velocity
v1 = v1_guess - delta_v1;

v1_guess = v1;

% Iterate
iter = iter+1;
err = max(abs(flow(1:3) - r2));

end

end

%-------------------------------------------------------------------------%

%-------------------------------------------------------------------------%

% function: lambertMR.m

function [A,P,E,ERROR,VI,VF,TPAR,THETA] = lambertMR(RI,RF,TOF,MU,orbitType,Nrev,Ncase,optionsLMR)

% lambertMR.m - Lambert's problem solver for all possible transfers
%   (multi-revolution transfer included).
%
% PROTOTYPE:
%   [A,P,E,ERROR,VI,VF,TPAR,THETA] = lambertMR(RI,RF,TOF,MU,orbitType,Nrev,Ncase,optionsLMR)
%
% DESCRIPTION:
%   Lambert's problem solver for all possible transfers:
%       1- zero-revolution (for all possible types of orbits: circles, ellipses,
%       	parabolas and hyperbolas)
%       2- multirevolution case
%       3- inversion of the motion
%
%   1- ZERO-REVOLUTION LAMBERT'S PROBLEM
%
%   For the solution of Lambert's problem with number of revolution = 0 the
%   subroutine by Chris D'Souza is included here.
%   This subroutine is a Lambert algorithm which given two radius vectors
%   and the time to get from one to the other, it finds the orbit
%   connecting the two. It solves the problem using a new algorithm
%   developed by R. Battin. It solves the Lambert problem for all possible
%   types of orbits (circles, ellipses, parabolas and hyperbolas).
%   The only singularity is for the case of a transfer angle of 360 degrees,
%   which is a rather obscure case.
%   It computes the velocity vectors corresponding to the given radius
%   vectors except for the case when the transfer angle is 180 degrees
%   in which case the orbit plane is ambiguous (an infinite number of
%   transfer orbits exist).
% 
%   2- MULTIREVOLUTION LAMBERT'S PROBLEM
%
%   For the solution of Lambert's problem with Nrev>0 number of revolution,
%   Battin's formulation has been extended to accomodate N-revolution
%   transfer orbits, by following the paper: "Using Battin Mathod to obtain 
%   Multiple-revolution Lambert's Solutions" by Shen and Tsiotras.
%
%   When Nrev>0 the possible orbits are just ellipses.
%   If 0<=Nrev<=Nmax, there are two Nrev-revolution transfer orbits.
%   These two transfer orbits have different semi-major axis and they may 
%   be all combinations of large-e and small-e transfer orbits.
%   The Original Successive Substitution Method by Battin converges to one
%   of the two possible solution with a viable initial guest, however it
%   diverges from the other one. Then a Reversed Successive Substitution is
%   used to converge to the second solution.
%   A procedure is implemented in order to guarantee to provide initial
%   guesses in the convergence region. If Nrew exceeds the maximum number
%   of revolution an ERROR is given:
%   warning('off','lambertMR:SuccessiveSubstitutionDiverged') to take out
%   the warnings or use optionsLMR(1) = 0.
% 
%   3- INVERSION OF THE MOTION
% 
%   Direct or retrograde option can be selected for the transfer.
%   
%   The algorithm computes the semi-major axis, the parameter (semi-latus 
%   rectum), the eccentricity and the velocity vectors.
% 
%   NOTE: If ERROR occurs or the 360 or 180 degree transfer case is 
%   encountered. 
%
% INPUT:
%	RI[3]           Vector containing the initial position in Cartesian
%                   coordinates [L].
%	RF[3]           Vector containing the final position vector in
%                   Cartesian coordinates [L].
%	TOF[1]          Transfer time, time of flight [T].
%  	MU[1]           Planetary constant of the planet (mu = mass * G) [L^3/T^2]
%	orbitType[1]    Logical variable defining whether transfer is
%                       0: direct transfer from R1 to R2 (counterclockwise)
%                       1: retrograde transfer from R1 to R2 (clockwise)
%	Nrev[1]         Number of revolutions.
%                   if Nrev = 0 ZERO-REVOLUTION transfer is calculated
%                   if Nrev > 0 two transfers are possible. Ncase should be
%                          defined to select one of the two.
%	Ncase[1]        Logical variable defining the small-a or large-a option
%                   in case of Nrev>0:
%                       0: small-a option
%                       1: large-a option
%	optionsLMR[1]	lambertMR options:
%                    optionsLMR(1) = display options:
%                                    0: no display
%                                    1: warnings are displayed only when
%                                       the algorithm does not converge
%                                    2: full warnings displayed
%
% OUTPUT:
%	A[1]        Semi-major axis of the transfer orbit [L].
% 	P[1]        Semi-latus rectum of the transfer orbit [L].
%  	E[1]        Eccentricity of the transfer orbit.
%	ERROR[1]	Error flag
%                   0:	No error
%                   1:	Error, routine failed to converge
%                   -1:	180 degrees transfer
%                   2:  360 degrees transfer
%                   3:  the algorithm doesn't converge because the number 
%                       of revolutions is bigger than Nrevmax for that TOF
%                   4:  Routine failed to converge, maximum number of
%                       iterations exceeded.
%	VI[3]       Vector containing the initial velocity vector in Cartesian
%               coordinates [L/T].
%	VT[1]		Vector containing the final velocity vector in Cartesian
%               coordinates [L/T].
%	TPAR[1] 	Parabolic flight time between RI and RF [T].
%	THETA[1]	Transfer angle [radians].
%
% NOTE: The semi-major axis, positions, times, and gravitational parameter
%       must be in compatible units.
%
% CALLED FUNCTIONS:
%   qck, h_E (added at the bottom of this file)
%
% REFERENCES:
%   - Shen and Tsiotras, "Using Battin method to obtain Multiple-Revolution
%       Lambert's solutions".
%   - Battin R., "An Introduction to the Mathematics and Methods of
%       Astrodynamics, Revised Edition", 1999.
%
% FUTURE DEVELOPMENT:
%   - 180 degrees transfer indetermination
%   - 360 degrees transfer singularity
%   - Nmax number of max revolution for a given TOF:
%     work in progress - Camilla Colombo
%
% ORIGINAL VERSION:
%   Chris D'Souza, 20/01/1989, MATLAB, lambert.m
%       verified by Darrel Monroe, 10/25/90
%       - Labert.m solved only direct transfer, without multi-revolution
%         option
%
% AUTHOR:
%   Camilla Colombo, 10/11/2006, MATLAB, lambertMR.m
%
% CHANGELOG:
%   13/11/2006, Camilla Colombo: added ERROR = 3 if Nrev > NrevMAX
%	21/11/2006, Camilla Colombo: added another case of ERROR = 3 (index
%   	N3) corresponding to the limit case when small-a solution = large-a
%       solution. No solution is given in this case.
%	06/08/2007, Camilla Colombo: optionsLMR added as an input
%	28/11/2007, Camilla Colombo: minor changes
%   29/01/2009, Matteo Ceriotti:
%       - Introduced variable for maximum number of iterations nitermax.
%       - Corrected final check on maximum number of iterations exceeded, from
%           "==" to ">=" (if N1 >= nitermax || N >= nitermax).
%       - Increased maxumum number of iterations to 2000, not to lose some
%           solutions.
%       - In OSS loop, added check for maximum number of iterations exceeded,
%           which then sets checkNconvOSS = 0.
%       - Changed the way of coumputing X given Y1 in RSS. Now the
%           Newton-Raphson method with initial guess suggested by Shen,
%           Tsiotras is used. This should guarantee convergence without the
%           need of an external zero finder (fsolve).
%       - Changed absolute tolerance into relative tolerance in all loops X0-X.
%           Now the condition is: while "abs(X0-X) >= abs(X)*TOL+TOL".
%       - Added return immediately when any error is detected.
%       - Moved check on 4*TOF*LAMBDA==0 after computing LAMBDA.
%       - Moved check on THETA==0 || THETA==2*PI after computing THETA.
%       - Added error code 4 (number of iterations exceeded).
%       - Removed variable Nwhile, as not strictly needed.
%       - Removed variable PIE=pi.
%   29/01/2009, REVISION: Matteo Ceriotti
%   21/07/2009, Matteo Ceriotti, Camilla Colombo:
%       added condition to detect case 180 degrees transfer indetermination
%   30/01/2010, Camilla Colombo: Header and function name in accordance
%       with guidlines.
%
% Note: Please if you have got any changes that you would like to be done,
%   do not change the function, please contact the author.
%
% -------------------------------------------------------------------------

% Check inputs
if nargin < 8
    optionsLMR = 0;
    if nargin < 6
        Nrev = 0;
        if nargin < 5
            orbitType = 0;
            if nargin < 4
                error('Not enough input arguments. See lambertMR.');
            end
        end
    end
end

nitermax = 2000; % Maximum number of iterations for loops
TOL = 1e-14;

TWOPI=2*pi;

% Reset
% A=0;P=0;E=0;
VI=[0,0,0];VF=[0,0,0];

% ----------------------------------
% Compute the vector magnitudes and various cross and dot products

RIM2   = dot(RI,RI);
RIM    = sqrt(RIM2);
RFM2   = dot(RF,RF);
RFM    = sqrt(RFM2);
CTH    = dot(RI,RF)/(RIM*RFM);
CR     = cross(RI,RF);
STH    = norm(CR)/(RIM*RFM);

% Choose angle for up angular momentum
switch orbitType
    case 0 % direct transfer
        if CR(3) < 0 
            STH = -STH;
        end
    case 1 % retrograde transfer
        if CR(3) > 0 
            STH = -STH;
        end
    otherwise
		error('%d is not an allowed orbitType',orbitType);
end
        
THETA  = qck(atan2(STH,CTH));
% if abs(THETA - pi) >= 0.01
if THETA == TWOPI || THETA==0
    ERROR = 2;
    A=0; P=0; E=0; VI=[0,0,0]; VF=[0,0,0]; TPAR=0; THETA=0;
    return
end

B1     = sign(STH); if STH == 0; B1 = 1; end

% ----------------------------------
% Compute the chord and the semi-perimeter

C= sqrt(RIM2 + RFM2 - 2*RIM*RFM*CTH);
S= (RIM + RFM + C)/2;
%BETA   = 2*asin(sqrt((S-C)/S));
%PMIN   = TWOPI*sqrt(S^3/(8*MU));
% TMIN   = PMIN*(pi-BETA+sin(BETA))/(TWOPI);
LAMBDA = B1*sqrt((S-C)/S);

if 4*TOF*LAMBDA == 0 || abs((S-C)/S) < TOL
    ERROR = -1;
    A=0; P=0; E=0; VI=[0,0,0]; VF=[0,0,0]; TPAR=0; THETA=0;
    return
end

% ----------------------------------
% Compute L carefully for transfer angles less than 5 degrees

if THETA*180/pi <= 5
   W   = atan((RFM/RIM)^.25) - pi/4;
   R1  = (sin(THETA/4))^2;
   S1  = (tan(2*W))^2;
   L   = (R1+S1)/(R1+S1+cos(THETA/2));
else
   L   = ((1-LAMBDA)/(1+LAMBDA))^2;
end

M= 8*MU*TOF^2/(S^3*(1+LAMBDA)^6);
TPAR   = (sqrt(2/MU)/3)*(S^1.5-B1*(S-C)^1.5);
L1     = (1 - L)/2;

CHECKFEAS = 0;
N1 = 0;
N = 0;

if Nrev == 0
    % ----------------------------------
    % Initialize values of y, n, and x

    Y= 1;
    N= 0;
    N1=0;
    ERROR  = 0;
    % CHECKFEAS=0;

    if (TOF-TPAR) <= 1e-3
        X0  = 0;
    else
        X0  = L;
    end

    X= -1.e8;

    % ----------------------------------
    % Begin iteration
    
    % ---> CL: 26/01/2009, Matteo Ceriotti: 
    %       Changed absolute tolerance into relative tolerance here below.
    while (abs(X0-X) >= abs(X)*TOL+TOL) && (N <= nitermax)
        N   = N+1;
        X   = X0;
        ETA = X/(sqrt(1+X) + 1)^2;
        CHECKFEAS=1;

        % ----------------------------------
        % Compute x by means of an algorithm devised by
        % Gauticci for evaluating continued fractions by the
        % 'Top Down' method
        
        DELTA = 1;
        U     = 1;
        SIGMA = 1;
        M1    = 0;

        while abs(U) > TOL && M1 <= nitermax
            M1    = M1+1;
            GAMMA = (M1 + 3)^2/(4*(M1+3)^2 - 1);
            DELTA = 1/(1 + GAMMA*ETA*DELTA);
            U     = U*(DELTA - 1);
            SIGMA = SIGMA + U;
        end

        C1 = 8*(sqrt(1+X)+1)/(3+1/(5 + ETA + (9*ETA/7)*SIGMA));

        % ----------------------------------
        % Compute H1 and H2
        
        if N == 1
            DENOM = (1 + 2*X + L)*(3*C1 + X*C1 +4*X);
            H1 = (L+X)^2*(C1 + 1 + 3*X)/DENOM;
            H2 = M*(C1+X-L)/DENOM;
        else
            QR = sqrt(L1^2 + M/Y^2);
            XPLL = QR - L1;
            LP2XP1 = 2*QR;
            DENOM = LP2XP1*(3*C1 + X*C1+4*X);
            H1 = ((XPLL^2)*(C1 + 1 + 3*X))/DENOM;
            H2 = M*(C1+X-L)/DENOM;
        end
        
        B = 27*H2/(4*(1+H1)^3);
        U = -B/(2*(sqrt(B+1)+1));

        % ----------------------------------
        % Compute the continued fraction expansion K(u)
        % by means of the 'Top Down' method
        
        % Y can be computed finding the roots of the formula and selecting
        % the real one:
        % y^3 - (1+h1)*y^2 - h2 = 0     (7.113) Battin
        %
        % Ycami_ = roots([1 -1-H1 0 -H2])
        % kcami = find( abs(imag(Ycami_)) < eps );
        % Ycami = Ycami_(kcami)

        DELTA = 1;
        U0 = 1;
        SIGMA = 1;
        N1 = 0;

        while N1 < nitermax && abs(U0) >= TOL
            if N1 == 0
                GAMMA = 4/27;
                DELTA = 1/(1-GAMMA*U*DELTA);
                U0 = U0*(DELTA - 1);
                SIGMA = SIGMA + U0;
            else
                for I8 = 1:2
                    if I8 == 1
                        GAMMA = 2*(3*N1+1)*(6*N1-1)/(9*(4*N1 - 1)*(4*N1+1));
                    else
                        GAMMA = 2*(3*N1+2)*(6*N1+1)/(9*(4*N1 + 1)*(4*N1+3));
                    end
                    DELTA = 1/(1-GAMMA*U*DELTA);
                    U0 = U0*(DELTA-1);
                    SIGMA = SIGMA + U0;
                end
            end

            N1 = N1 + 1;
        end

        KU = (SIGMA/3)^2;
        Y = ((1+H1)/3)*(2+sqrt(B+1)/(1-2*U*KU));    % Y = Ycami
        
        X0 = sqrt(((1-L)/2)^2+M/Y^2)-(1+L)/2;
        % fprintf('n= %d, x0=%.14f\n',N,X0);
    end
    
% MULTIREVOLUTION
elseif (Nrev > 0) && (4*TOF*LAMBDA~=0) %(abs(THETA)-pi > 0.5*pi/180)

    checkNconvRSS = 1;
    checkNconvOSS = 1;
    N3 = 1;
    
    while N3 < 3
        
        if Ncase == 0 || checkNconvRSS == 0

            % - Original Successive Substitution -
            % always converges to xL - small a

            % ----------------------------------
            % Initialize values of y, n, and x
            
            Y= 1;
            N= 0;
            N1=0;
            ERROR = 0;
            % CHECKFEAS = 0;
%             if (TOF-TPAR) <= 1e-3
%                 X0 = 0;
%             else
            if checkNconvOSS == 0
                X0 = 2*X0;
                checkNconvOSS = 1;
                % see p. 11 USING BATTIN METHOD TO OBTAIN 
                % MULTIPLE-REVOLUTION LAMBERT'S SOLUTIONS - Shen, Tsiotras
            elseif checkNconvRSS == 0
                % X0 is taken from the RSS
            else
                X0 = L;
            end

            X = -1.e8;

            % ----------------------------------
            % Begin iteration
            			
            % ---> CL: 26/01/2009,Matteo Ceriotti 
            %   Changed absolute tolerance into relative tolerance here
            %   below.
            while (abs(X0-X) >= abs(X)*TOL+TOL) && (N <= nitermax)
                N   = N+1;
                X   = X0;
                ETA = X/(sqrt(1+X) + 1)^2;
                CHECKFEAS = 1;

                % ----------------------------------
                % Compute x by means of an algorithm devised by
                % Gauticci for evaluating continued fractions by the
                % 'Top Down' method
                

                DELTA = 1;
                U     = 1;
                SIGMA = 1;
                M1    = 0;

                while abs(U) > TOL && M1 <= nitermax
                    M1    = M1+1;
                    GAMMA = (M1 + 3)^2/(4*(M1+3)^2 - 1);
                    DELTA = 1/(1 + GAMMA*ETA*DELTA);
                    U     = U*(DELTA - 1);
                    SIGMA = SIGMA + U;
                end

                C1 = 8*(sqrt(1+X)+1)/(3+1/(5 + ETA + (9*ETA/7)*SIGMA));

                % ----------------------------------
                % Compute H1 and H2
                
                if N == 1
                    DENOM = (1 + 2*X + L)*(3*C1 + X*C1 +4*X);
                    H1 = (L+X)^2*(C1 + 1 + 3*X)/DENOM;
                    H2 = M*(C1+X-L)/DENOM;
                else
                    QR = sqrt(L1^2 + M/Y^2);
                    XPLL = QR - L1;
                    LP2XP1 = 2*QR;
                    DENOM = LP2XP1*(3*C1 + X*C1+4*X);
                    H1 = ((XPLL^2)*(C1 + 1 + 3*X))/DENOM;
                    H2 = M*(C1+X-L)/DENOM;
                end

                H3 = M*Nrev*pi/(4*X*sqrt(X));
                H2 = H3+H2;

                B = 27*H2/(4*(1+H1)^3);
                U = -B/(2*(sqrt(B+1)+1));

                % ----------------------------------
                % Compute the continued fraction expansion K(u)
                % by means of the 'Top Down' method
                
                % Y can be computed finding the roots of the formula and selecting
                % the real one:
                % y^3 - (1+h1)*y^2 - h2 = 0     (7.113) Battin
                %
                % Ycami_ = roots([1 -1-H1 0 -H2])
                % kcami = find( abs(imag(Ycami_)) < eps );
                % Ycami = Ycami_(kcami)

                DELTA = 1;
                U0 = 1;
                SIGMA = 1;
                N1 = 0;

                while N1 < nitermax && abs(U0) >= TOL
                    if N1 == 0
                        GAMMA = 4/27;
                        DELTA = 1/(1-GAMMA*U*DELTA);
                        U0 = U0*(DELTA - 1);
                        SIGMA = SIGMA + U0;
                    else
                        for I8 = 1:2
                            if I8 == 1
                                GAMMA = 2*(3*N1+1)*(6*N1-1)/(9*(4*N1 - 1)*(4*N1+1));
                            else
                                GAMMA = 2*(3*N1+2)*(6*N1+1)/(9*(4*N1 + 1)*(4*N1+3));
                            end
                            DELTA = 1/(1-GAMMA*U*DELTA);
                            U0 = U0*(DELTA-1);
                            SIGMA = SIGMA + U0;
                        end
                    end

                    N1 = N1 + 1;
                end

                KU = (SIGMA/3)^2;
                Y = ((1+H1)/3)*(2+sqrt(B+1)/(1-2*U*KU));	% Y = Ycami
                if Y > sqrt(M/L)
                    if optionsLMR(1) == 2
                        warning('lambertMR:SuccessiveSubstitutionDiverged',...
                                ['Original Successive Substitution is diverging\n'...
                                '-> Reverse Successive Substitution used to find the proper XO.\n']);
                    end
                    checkNconvOSS = 0;
                    break
                end
                
                X0 = sqrt(((1-L)/2)^2+M/Y^2)-(1+L)/2;
                % fprintf('N: %d X0: %.14f\n',N,X0);
            end
            
            % When 2 solutions exist (small and big a), the previous loop
            % must either converge or diverge because Y > sqrt(M/L) at some
            % point. Thus, the upper bound on the number of iterations
            % should not be necessary. Though, nothing can be said in the
            % case tof<tofmin and so no solution exist. In this case, an
            % upper bound on number of iterations could be needed.
            
            if N >= nitermax % Checks if previous loop ended due to maximum number of iterations
                if optionsLMR(1) == 2
                    warning('lambertMR:SuccessiveSubstitutionExceedMaxIter',...
                            ['Original Successive Substitution exceeded max number of iteration\n'...
                            '-> Reverse Successive Substitution used to find the proper XO.\n']);
                end
                checkNconvOSS = 0;
            end
        end
        if (Ncase == 1 || checkNconvOSS == 0) && ~(checkNconvRSS == 0 && checkNconvOSS == 0)

            % - Reverse Successive Substitution -
            % always converges to xR - large a

            % ----------------------------------
            % Initialize values of y, n, and x
            
            N = 0;
            N1 = 0;
            ERROR  = 0;
            % CHECKFEAS=0;
            if checkNconvRSS == 0
                X0 = X0/2; % XL/2
                checkNconvRSS = 1;
                % see p. 11 USING BATTIN METHOD TO OBTAIN 
                % MULTIPLE-REVOLUTION LAMBERT'S SOLUTIONS - Shen, Tsiotras
            elseif checkNconvOSS == 0
                % X0 is taken from the OSS
            else
                X0 = L;
            end

            X = -1.e8;

            % ----------------------------------
            % Begin iteration
            
            % ---> CL: 26/01/2009, Matteo Ceriotti
            %   Changed absolute tolerance into relative tolerance here
            %   below.
            while (abs(X0-X) >= abs(X)*TOL+TOL) && (N <= nitermax)
                N = N+1;
                X = X0;
                CHECKFEAS=1;

                Y = sqrt(M/((L+X)*(1+X))); % y1 in eq. (8a) in Shen, Tsiotras

                if Y < 1
                    if optionsLMR(1) == 2
                        warning('lambertMR:SuccessiveSubstitutionDiverged',...
                                ['Reverse Successive Substitution is diverging\n' ...
                                '-> Original Successive Substitution used to find the proper XO.\n']);
                    end
                    checkNconvRSS = 0;
                    break
                end
                
                % ---> CL: 27/01/2009, Matteo Ceriotti
                %   This is the Newton-Raphson method suggested by USING
                %   BATTIN METHOD TO OBTAIN MULTIPLE-REVOLUTION LAMBERT'S
                %   SOLUTIONS - Shen, Tsiotras
                
                % To assure the Newton-Raphson method to be convergent
                Erss = 2*atan(sqrt(X));
                while h_E(Erss,Y,M,Nrev) < 0
                    Erss = Erss/2;
                end
                
                Nnew = 1;
                Erss_old = -1.e8;
                
                % The following Newton-Raphson method should always
                % converge, given the previous first guess choice,
                % according to the paper. Therefore, the condition on
                % number of iterations should not be neccesary. It could be
                % necessary for the case tof < tofmin.
                while (abs(Erss-Erss_old) >= abs(Erss)*TOL+TOL) && Nnew < nitermax
                    Nnew = Nnew+1;
                    [h, dh] = h_E(Erss,Y,M,Nrev);
                    Erss_old = Erss;
                    Erss = Erss - h/dh;
                    % fprintf('Nnew: %d Erss: %.16f h_E: %.16f\n',Nnew,Erss,h);
                end
                if Nnew >= nitermax
                    if optionsLMR(1) ~= 0
                        warning('lambertMR:NewtonRaphsonIterExceeded', 'Newton-Raphson exceeded max iterations.\n');
                    end
                end
                X0 = tan(Erss/2)^2;
            end
        end
        if checkNconvOSS == 1 && checkNconvRSS == 1
            break
        end
        
        if checkNconvRSS == 0 && checkNconvOSS == 0
            if optionsLMR ~=0
                warning('lambertMR:SuccessiveSubstitutionDiverged',...
                        ['Both Original Successive Substitution and Reverse ' ...
                        'Successive Substitution diverge because Nrev > NrevMAX.\n' ...
                        'Work in progress to calculate NrevMAX.\n']);
            end
            ERROR = 3;
            A=0; P=0; E=0; VI=[0,0,0]; VF=[0,0,0]; TPAR=0; THETA=0;
            return
        end
        
        N3 = N3+1;
    end
    
    if N3 == 3
        if optionsLMR ~=0
            warning('lambertMR:SuccessiveSubstitutionDiverged',...
                    ['Either Original Successive Substitution or Reverse ' ...
                    'Successive Substitution is always diverging\n' ...
                    'because Nrev > NrevMAX or because large-a solution = small-a solution (limit case).\n' ...
                    'Work in progress to calculate NrevMAX.\n']);
        end
        ERROR = 3;
        A=0; P=0; E=0; VI=[0,0,0]; VF=[0,0,0]; TPAR=0; THETA=0;
        return
    end
end

% ----------------------------------
% Compute the velocity vectors

if CHECKFEAS == 0
    ERROR = 1;
    A=0; P=0; E=0; VI=[0,0,0]; VF=[0,0,0]; TPAR=0; THETA=0;
    return
end

if N1 >= nitermax || N >= nitermax
    ERROR = 4;
    if optionsLMR ~=0
        disp('Lambert algorithm has not converged, maximum number of iterations exceeded.');
    end
    A=0; P=0; E=0; VI=[0,0,0]; VF=[0,0,0]; TPAR=0; THETA=0;
    return
end

CONST = M*S*(1+LAMBDA)^2;
A = CONST/(8*X0*Y^2);

R11 = (1 + LAMBDA)^2/(4*TOF*LAMBDA);
S11 = Y*(1 + X0);
T11 = (M*S*(1+LAMBDA)^2)/S11;

VI(1:3) = -R11*(S11*(RI(1:3)-RF(1:3))-T11*RI(1:3)/RIM);
VF(1:3) = -R11*(S11*(RI(1:3)-RF(1:3))+T11*RF(1:3)/RFM);

P = (2*RIM*RFM*Y^2*(1+X0)^2*sin(THETA/2)^2)/CONST;
E = sqrt(1 - P/A);

return

% -------------------------------------------------------------------------

function [angle] = qck(angle)

% qck.m - Reduce an angle between 0 and 2*pi
%
% PROTOTYPE:
%   [angle]=qck(angle)
%
% DESCRIPTION:
%   This function takes any angle and reduces it, if necessary,
% 	so that it lies in the range from 0 to 2 PI radians.
% 
% INPUTS:
%   ANGLE[1]    Angle to be reduced (in radians)
% 
% OUTPUTS:
%   QCK[1]      The angle reduced, if necessary, to the range
%               from 0 to 2 PI radians (in radians)
% 
% CALLED FUNCTIONS:
%   pi (from MATLAB)
%
% AUTHOR:
%   W.T. Fowler, July, 1978
%
% CHANGELOG:
%   8/20/90, REVISION: Darrel Monroe
%
% -------------------------------------------------------------------------

twopi = 2*pi;
 
diff = twopi * (fix(angle/twopi) + min([0,sign(angle)]));

angle = angle -diff;

return

% -------------------------------------------------------------------------

end
end

%-------------------------------------------------------------------------%    

%-------------------------------------------------------------------------%

% function: drawBody.m

function [globe] = drawBody(body_label, position, handle, scaleFactor)
%-------------------------------------------------------------------------%
% 
% drawBody.m draws a celestial body (Planet, Sun, Moon) of the Solar System
% in a 3D plot.
%
%-------------------------------------------------------------------------%
% PROTOTYPE:
%  [globe] = drawBody(body_label, position, handle, scaleFactor)
%
%-------------------------------------------------------------------------%
% INPUT ARGUMENTS:
%  body_label   [char]  String identifying the celestial body 
%                       1:  'Mercury'
%                       2:  'Venus'
%                       3:  'Earth'
%                       4:  'Mars'
%                       5:  'Jupiter'
%                       6:  'Saturn'
%                       7:  'Uranus'
%                       8:  'Neptune'
%                       9:  'Pluto'
%                       10: 'Sun'
%                       11: 'Moon'
%  position    [3x1]    Position of the celestial body
%  handle      [1]      Figure handle
%  scaleFactor [1]      Scale factor of the celestial body
%
%-------------------------------------------------------------------------%
% OUTPUT ARGUMENTS:
%  [globe]     [1x1 Surface] 3D surface of the celestial body
%
%-------------------------------------------------------------------------%
% CALLED FUNCTIONS:
%   (none)
%
%-------------------------------------------------------------------------%
% CONTRIBUTORS:
%   Unknown original author
%   Gian Marco Paldino
% 
%-------------------------------------------------------------------------%
% VERSIONS: 
%   Unknown date, First version: plotPlanet.m
%   13/12/2020, Gian Marco Paldino: code cleanup and revision, added Moon,
%               function name changed, header created in accordance to
%               guidelines.
%   29/10/2021, Gian Marco Paldino: use SPICE function cspice_bodvrd to
%               retrieve constants.
%
%-------------------------------------------------------------------------%
   
% Preparing Figure Object
    if nargin<3
        HAXIS = gca;
    elseif ishandle(handle)==0
            msg = ('The figure handle is not valid');
            error(msg)
    else
        try
            HAXIS=gca(handle);
        catch
            HAXIS=handle;  
        end
        hold on
    end
    
    if nargin<4
        if body_label == Sun
            scaleFactor = 1;
        else
            radius = cspice_bodvrd(body_label,'RADII',3);
            radius_sun = cspice_bodvrd('Sun','RADII',3);
            scaleFactor = radius(1)/radius_sun(1);
        end
    end
    
    radius_sun = cspice_bodvrd('Sun','RADII',3);
    % Planet radius w.r.t. sun [km]
    Rplanet = radius_sun(1)*scaleFactor; 
    % Number of globe panels around the equator [deg/panel] = [360/npanels]
    npanels = 360;
    erad=Rplanet; % Equatorial radius [km]
    prad=Rplanet; % Polar radius [km]
    hold on;
    axis equal;
    axis vis3d;
    % Create a 3D meshgrid of the sphere points using the ellipsoid function
    [x,y,z] = ellipsoid(position(1), position(2), position(3), erad, erad, prad, npanels);
    globe = surf(HAXIS, x,y,z,'FaceColor','none','EdgeColor',0.5*[1 1 1], 'HandleVisibility','off');
    % RMK.: HandleVisibility=off removes unuseful legends for the plotted
    % globe
    % Load Earth image for texture map
    cdata = imread(sprintf('%s.jpg',body_label)); 
    % Set image as color data (cdata) property, and set face color to 
    % indicate a texturemap, which Matlab expects to be in cdata.

    globe.FaceColor = 'texturemap';
    globe.CData = flip(cdata); % W/o flip() the texture looks upside down
    globe.EdgeColor = 'none';    
    globe.FaceLighting = 'gouraud';
    globe.AmbientStrength = 0.5;
    
     if body_label ~= 10
        globe.FaceAlpha = 0.9;
     else
        globe.FaceAlpha = 1.0;
        globe.FaceLighting = 'none';
     end
      
    rotate(globe,[0 0 1],180, position);
    
end 

%-------------------------------------------------------------------------%    

%-------------------------------------------------------------------------%

% function: E2M_ToF_boundaries.m

function [TPAR, T_hohmann] = E2M_ToF_boundaries(dep_body,arr_body,primary_body,mu_primary,ref_frame,et)
%-------------------------------------------------------------------------%
%
% E2M_ToF_boundaries.m computes time of flight lower and upper boundary
% to be used as constraints for the optimization algorithm for a Earth to 
% Mars transfer.
%
%-------------------------------------------------------------------------%
% PROTOTYPE:
%  [TPAR, T_hohmann] = E2M_ToF_boundaries(dep_body,arr_body,primary_body,mu_primary,ref_frame,et)
%
%-------------------------------------------------------------------------%
% INPUT ARGUMENTS:
%  dep_body         [char]  Departure planet label (Earth)       [-]
%  arr_body         [char]  Arrival planet label   (Mars)        [-]
%  primary_body     [char]  Primary body label     (Sun)         [-]
%  mu_primary       [1]     Primary body standard gravitational  [km^3/s^2]  
%                           constant
%  ref_frame        [char]  Reference frame string               [-]
%  et               [1]     Ephemeris time                       [s]
%
%-------------------------------------------------------------------------%
% OUTPUT ARGUMENTS:
%  TPAR             [1]  Parabolic ToF between departure planet  [s]
%                        and arrival planet
%  T_hohmann        [1]  Hohmann ToF between departure planet    [s]
%                        and arrival planet
%-------------------------------------------------------------------------%
% AUTHOR:
%  Gian Marco Paldino
%
%-------------------------------------------------------------------------%
% VERSIONS: 
%  29/10/2021: first version
%
%-------------------------------------------------------------------------%

% Get state of departure planet
State_dep = cspice_spkezr(dep_body,et,ref_frame,'NONE',primary_body);
rr_dep = State_dep(1:3);
vv_dep = State_dep(4:6);

% Convert the state to keplerian elements
[kep_dep(1),~,~,~,~,~] = car2kep(rr_dep,vv_dep,mu_primary);

% Get state of arrival planet
State_arr = cspice_spkezr(arr_body,et,ref_frame,'NONE',primary_body);
rr_arr = State_arr(1:3);
vv_arr = State_arr(4:6);

% Convert the state to keplerian elements
[kep_arr(1),~,~,~,~,~] = car2kep(rr_arr,vv_arr,mu_primary);

% Define the apocenter radius and pericenter radius of the Hohmann transfer
% orbit
r_p = kep_dep(1);
r_a = kep_arr(1);

% Compute the semi-major axis [km] and the half of the period of the 
% Hohmann transfer orbit [s]
a_hohmann = (r_a+r_p)*0.5;
T_hohmann = pi*sqrt((a_hohmann^3)/mu_primary);

% Compute the parabolic ToF using the classic Lambert solver
[~,~,~,~,~,~,TPAR,~] = lambertMR(rr_dep, rr_arr, T_hohmann, mu_primary, 0, 0, 0, 1);

end

%-------------------------------------------------------------------------%    

%-------------------------------------------------------------------------%

% function: E2M_objectivefun.m

function DeltaV_tot = E2M_objectivefun(Y,dep_body,arr_body,primary_body,mu_primary,ref_frame)
%-------------------------------------------------------------------------%
%
% E2M_objectivefun.m provides the objective function for a fuel optimal 
% Earth to Mars transfer, computing the total DeltaV required by the
% transfer.
%
%-------------------------------------------------------------------------%
% PROTOTYPE:
%  DeltaV_tot = E2M_objectivefun(Y,dep_body,arr_body,primary_body,... 
%                                mu_primary,ref_frame)
%
%-------------------------------------------------------------------------%
% INPUT ARGUMENTS:
%  Y                [2x1]   Array with departure and arrival     [d]
%                           ephemeris times 
%  dep_body         [char]  Departure planet label (Earth)       [-]
%  arr_body         [char]  Arrival planet label   (Mars)        [-]
%  primary_body     [char]  Primary body label     (Sun)         [-]
%  mu_primary       [1]     Primary body standard gravitational  [km^3/s^2]  
%                           constant
%  ref_frame        [char]  Reference frame string               [-]
%
%-------------------------------------------------------------------------%
% OUTPUT ARGUMENT:
%  DeltaV_tot       [1]     Total DeltaV of the transfer         [km/s]
%
%-------------------------------------------------------------------------%
% AUTHOR:
%  Gian Marco Paldino
%
%-------------------------------------------------------------------------%
% VERSIONS: 
%  29/10/2021: first version
%
%-------------------------------------------------------------------------%

% Retrieve departure and arrival ephemeris times from the input vector
et_dep = Y(1);
et_arr = Y(2);

% Compute the state of Earth and Mars at departure and arrival time
% respectively
State_dep = cspice_spkezr(dep_body,et_dep*24*3600,ref_frame,'NONE',primary_body);
State_arr = cspice_spkezr(arr_body,et_arr*24*3600,ref_frame,'NONE',primary_body);

% Characterize transfer trajectory
[~,~,~,~,vv_T1,vv_T2,~,~] = lambertMR(State_dep(1:3)', State_arr(1:3)', (et_arr-et_dep)*24*3600, mu_primary, 0, 0, 2);

% Compute DeltaV of the transfer
Deltav_1 = norm(vv_T1' - State_dep(4:6));
Deltav_2 = norm(vv_T2' - State_arr(4:6));
DeltaV_tot = Deltav_1 + Deltav_2;

end

%-------------------------------------------------------------------------%    

%-------------------------------------------------------------------------%

% function nbody_init.m

function [bodies] = nbody_init(labels)
%NBODY_INIT Initialize planetary data for n-body propagation
%   Given a set of labels of planets and/or barycentres, returns a
%   cell array populated with structures containing the body label and the
%   associated gravitational constant.
%
%
% Author
%   Name: ALESSANDRO 
%   Surname: MORSELLI
%   Research group: DART
%   Department: DAER
%   University: Politecnico di Milano 
%   Creation: 26/09/2021
%   Contact: alessandro.morselli@polimi.it
%   Copyright: (c) 2021 A. Morselli, Politecnico di Milano. 
%                  All rights reserved.
%
%
% Notes:
%   This material was prepared to support the course 'Satellite Guidance
%   and Navigation', AY 2021/2022.
%
%
% Inputs:
%   labels : [1,n] cell-array with object labels
%
% Outputs:
%   bodies : [1,n] cell-array with struct elements containing the following
%                  fields
%                  |
%                  |--bodies{i}.name -> body label
%                  |--bodies{i}.GM   -> gravitational constant [km**3/s**2]
%
%
% Prerequisites:
%   - MICE (Matlab SPICE)
%   - Populated kernel pool (PCK kernels)
%

% Initialize output
bodies = cell(size(labels));

% Loop over labels
for i = 1:length(labels)
    % Store body label
    bodies{i}.name = labels{i};
    % Store body gravitational constant
    bodies{i}.GM   = cspice_bodvrd(labels{i}, 'GM', 1);
end

end

%-------------------------------------------------------------------------%

%-------------------------------------------------------------------------%

% function: nbody_rhs.m

function [dxdt] = nbody_rhs(t, x, bodies, frame)
%NBODY_RHS Evaluates the right-hand-side of a N-body propagator
%   Evaluates the right-hand-side of a newtonian N-body propagator.
%   The integration centre is the Solar-System-Barycentre (SSB) and only
%   Newtonian gravitational accelerations are considered.
%
%
% Author
%   Name: ALESSANDRO 
%   Surname: MORSELLI
%   Research group: DART
%   Department: DAER
%   University: Politecnico di Milano 
%   Creation: 26/09/2021
%   Contact: alessandro.morselli@polimi.it
%   Copyright: (c) 2021 A. Morselli, Politecnico di Milano. 
%                  All rights reserved.
%
%
% Notes:
%   This material was prepared to support the course 'Satellite Guidance
%   and Navigation', AY 2021/2022.
%
%
% Inputs:
%   t      : [1,1] ephemeris time (ET SPICE), seconds past J2000 (TDB)
%   x      : [6,1] cartesian state vector wrt Earth-Barycentre
%   bodies : [1,6] cell-array created with function nbody_init
%
% Outputs:
%   dxdt   : [6,1] RHS, newtonian gravitational acceleration only
%
% Prerequisites:
%   - MICE (Matlab SPICE)
%   - Populated kernel pool (SPK, LSK, PCK kernels)
%
%-------------------------------------------------------------------------%
% VERSIONS: 
%  29/09/2021: Alessandro Morselli, first version
%  05/11/2021: Gian Marco Paldino, second version:
%              Modified function in order to work w.r.t. Earth 
%              instead of Solar System Barycenter.
%
%-------------------------------------------------------------------------%

if not( strcmpi(frame, 'ECLIPJ2000') || strcmpi(frame, 'J2000') )
    msg = 'Invalid integration reference frame, select either J2000 or ECLIPJ2000';
    error(msg);
end

% Initialize right-hand-side
dxdt = zeros(6,1);

% Position derivative is object's velocity
dxdt(1:3) = x(4:6);

% Extract the object position at state x
rr_earth_obj = x(1:3);

% Compute square distance and distance
norm_rr2 = dot(rr_earth_obj,rr_earth_obj);
norm_rr = sqrt(norm_rr2);

% Compute primary body acceleration
dxdt(4:6) = -bodies{1}.GM*rr_earth_obj/(norm_rr*norm_rr2);

% Compute perturbing accelerations from the other bodies
for i=2:length(bodies)

    % Retrieve position of i-th celestial body w.r.t. Earth
    % in inertial frame
    rho_body_earth = cspice_spkpos(bodies{i}.name, t, frame, 'NONE', 'EARTH');
    
    % Compute square distance and distance
    rho2 = dot(rho_body_earth,rho_body_earth);
    rho = sqrt(rho2);
    
    % Extract object position wrt. i-th celestial body
    dd_body_obj = rr_earth_obj - rho_body_earth;
    
    % Compute square distance and distance
    dd2 = dot(dd_body_obj, dd_body_obj);
    dd = sqrt(dd2);
    
    % Compute the gravitational acceleration using Newton's law
    aa_grav =  - bodies{i}.GM * (dd_body_obj /(dd*dd2) + rho_body_earth/(rho*rho2));

    % Sum up acceleration to right-hand-side
    dxdt(4:6) = dxdt(4:6) + aa_grav;

end

end

%-------------------------------------------------------------------------%

%-------------------------------------------------------------------------%

% function flow_NBP.m

function Xt = flow_NBP(rr0,vv0,ti,tf,bodies,ref_frame,outflag)
%-------------------------------------------------------------------------%
%
% flow_NBP.m computes the flow of the NBP.
%
%-------------------------------------------------------------------------%
% PROTOTYPE:
%  Xt = flow_NBP(rr0,vv0,ti,tf,bodies,ref_frame,0)
%  Xt = flow_NBP(rr0,vv0,ti,tf,bodies,ref_frame,1)
%
%-------------------------------------------------------------------------%
% INPUT ARGUMENTS:
%  rr0             [3x1]  Initial position                   [km]
%  vv0             [3x1]  Initial velocity                   [km/s]
%  ti              [1]    Initial time                       [s]
%  tf              [1]    Final time                         [s]
%  bodies          [1xn]  Cell-array created with function   [-]
%                         nbody_init
%  ref_frame       [char] Reference frame string             [-]
%  outflag         [1]    Output flag: position only (0)     [-]
%                         or full state (1)     
%                   
%-------------------------------------------------------------------------%
% OUTPUT ARGUMENTS: 
%  Xt              [6x1]  Flow of the NBP at time tf        [km, km/s]
% 
%-------------------------------------------------------------------------%
% AUTHOR:
%  Gian Marco Paldino
%
%-------------------------------------------------------------------------%
% VERSIONS: 
%  05/11/2021: First version
%
%-------------------------------------------------------------------------%

% Set integration timespan 
tspan = [ti tf];

% Set integration options
options = odeset('RelTol', 1e-13, 'AbsTol', 1e-13);

% Set initial conditions for the integration of the ODE  
X0 = [rr0; vv0];

% Perform the integration
[~, X] = ode113(@(t,x) nbody_rhs(t,x,bodies,ref_frame), tspan, X0, options);

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

% function E2EML2_DeltaV.m

function [DeltaV_tot,DeltaV_1,DeltaV_2,DeltaV_3] = E2EML2_DeltaV(Y,bodies,ref_frame)
%-------------------------------------------------------------------------%
%
% E2EML2_DeltaV.m computes the DeltaV for a Earth to Earth-Moon L2 transfer
% with three impulses.
%
%-------------------------------------------------------------------------%
% PROTOTYPE:
%  [DeltaV_tot,DeltaV_1,DeltaV_2,DeltaV_3] = E2EML2_DeltaV(Y,bodies,ref_frame)
%
%-------------------------------------------------------------------------%
% INPUT ARGUMENTS:
%  Y                [15x1]  Optimization variable:
%                             r1       [3x1] Initial position  [AU]   
%                             v1       [3x1] Initial velocity  [km/s]  
%                             r2       [3x1] Position at t2    [AU] 
%                             v2       [3x1] Velocity at t2    [km/s]  
%                             t1       [1]   Departure time    [d]
%                             t2       [1]   Second burn time  [d]
%                             t3       [1]   Arrival time      [d]
%  bodies           [1,6]   Cell-array created with function   [-]
%                           nbody_init
%  ref_frame        [char]  Reference frame string             [-] 
%
%-------------------------------------------------------------------------%
% OUTPUT ARGUMENT:
%  DeltaV_tot      [1]      Total cost of the transfer         [km/s]
%  DeltaV_1        [1]      Cost of the first impulse          [km/s]
%  DeltaV_2        [1]      Cost of the second impulse         [km/s]
%  DeltaV_3        [1]      Cost of the third impulse          [km/s]
%
%-------------------------------------------------------------------------%
% AUTHOR:
%  Gian Marco Paldino
%
%-------------------------------------------------------------------------%
% VERSIONS: 
%  05/11/2021: first version
%
%-------------------------------------------------------------------------%

% Seconds per day
spd = cspice_spd;

% Astronomical Unit in [km]
AU = 149597870.691;

% Retrieve optimization variables
X1 = Y(1:6);
X2 = Y(7:12);
t1 = Y(13)*spd;
t2 = Y(14)*spd;
t3 = Y(15)*spd;

% Retrieve state at time t1 (position [km] and velocity [km/s])
rr1 = X1(1:3)*AU;
vv1 = X1(4:6);

% Retrieve state at time t2 (position [km] and velocity [km/s])
rr2 = X2(1:3)*AU;
vv2 = X2(4:6);

% Set Earth standard gravitational parameter [km^3/s^2]
mu_earth = bodies{1}.GM;

% Compute Moon state vector at final time
X_moon = cspice_spkezr(bodies{3}.name,t3,ref_frame,'NONE',bodies{1}.name);

% EML2 velocity is the same as the Moon velocity
vv_L2 = X_moon(4:6);

% Compute the flows from t1 to t2 and from t2 to t3
X2_flow = flow_NBP(rr1,vv1,t1,t2,bodies,ref_frame,1);
X3_flow = flow_NBP(rr2,vv2,t2,t3,bodies,ref_frame,1);

% Compute velocity vector along the initial parking orbit
vv_park = cross([0;0;1],rr1/norm(rr1))*sqrt(mu_earth/norm(rr1));

% Compute DeltaV of the the transfer
DeltaV_1 = norm(vv1 - vv_park);
DeltaV_2 = norm(vv2 - X2_flow(4:6));
DeltaV_3 = norm(vv_L2 - X3_flow(4:6));

DeltaV_tot = DeltaV_1 + DeltaV_2 + DeltaV_3;

end

%-------------------------------------------------------------------------%

%-------------------------------------------------------------------------%

% function E2EML2_constraints.m

function [c, ceq] = E2EML2_constraints(Y,bodies,ref_frame)
%-------------------------------------------------------------------------%
%
% E2EML2_constraints.m computes the nonlinear constraints for a Earth to 
% Earth-Moon L2 transfer with three impulses.
%
%-------------------------------------------------------------------------%
% PROTOTYPE:
%  [c, ceq] = E2EML2_constraints(Y,bodies,ref_frame)
%
%-------------------------------------------------------------------------%
% INPUT ARGUMENTS:
%  Y                [15x1]  Optimization variable:
%                             r1       [3x1] Initial position  [AU]   
%                             v1       [3x1] Initial velocity  [km/s]  
%                             r2       [3x1] Position at t2    [AU] 
%                             v2       [3x1] Velocity at t2    [km/s]  
%                             t1       [1]   Departure time    [d]
%                             t2       [1]   Second burn time  [d]
%                             t3       [1]   Arrival time      [d]
%  bodies           [1,6]   Cell-array created with function   [-]
%                           nbody_init
%  ref_frame        [char]  Reference frame string             [-] 
%
%-------------------------------------------------------------------------%
% OUTPUT ARGUMENT:
%  c               [1]      Nonlinear inequality constraints 
%                           r1 >= Earth radius
%                           r(t2) >= 500000 [km]
%                           v1 <= v_escape
%                           v(t1) >= v_circ
%  ceq             [1]      Nonlinear equality constraints         
%                           r(t1) = r_park
%                           r(t2) = r2
%                           r(t3) = rL2
%-------------------------------------------------------------------------%
% AUTHOR:
%  Gian Marco Paldino
%
%-------------------------------------------------------------------------%
% VERSIONS: 
%  05/11/2021: first version
%
%-------------------------------------------------------------------------%
% Seconds per day
spd = cspice_spd;

% Astronomical Unit in [km]
AU = 149597870.691;

% Retrieve optimization variables
X1 = Y(1:6);
X2 = Y(7:12);
t1 = Y(13)*spd;
t2 = Y(14)*spd;
t3 = Y(15)*spd;

% Retrieve state at time t1 (position [km] and velocity [km/s])
rr1 = X1(1:3)*AU;
vv1 = X1(4:6);

% Retrieve state at time t2 (position [km] and velocity [km/s])
rr2 = X2(1:3)*AU;
vv2 = X2(4:6);

% Set Earth standard gravitational parameter [km^3/s^2]
mu_earth = bodies{1}.GM;

% Set Earth radius 
Radii_earth = cspice_bodvrd(bodies{1}.name, 'RADII', 3);
R_earth = Radii_earth(1);

% Compute initial parking orbit radius
r_park = R_earth + 200;

% Compute Moon state vector at final time
X_moon = cspice_spkezr(bodies{3}.name,t3,ref_frame,'NONE',bodies{1}.name);

% Moon position vector at final time [km]
rr_moon = X_moon(1:3);

% EML2 Position vector at final time [km] (same as Moon + 60000 km)
rr_L2 = (norm(rr_moon) + 60000) * rr_moon/norm(rr_moon);

% Compute nonlinear inequality constraints
c = [R_earth - norm(rr1);
     500000 - norm(rr2);
     norm(vv1) - sqrt(2*mu_earth/norm(rr1))
     sqrt(mu_earth/norm(rr1)) - norm(vv1)];
 
% Compute nonlinear equality constraints
ceq = [norm(rr1) - r_park;
       flow_NBP(rr1,vv1,t1,t2,bodies,ref_frame,0) - rr2;
       flow_NBP(rr2,vv2,t2,t3,bodies,ref_frame,0) - rr_L2];
   
end

%-------------------------------------------------------------------------%

%-------------------------------------------------------------------------%

% function TPBVP_dynamics

function dY = TPBVP_dynamics(~,Y,mu,Tmax,Isp,g0)
%-------------------------------------------------------------------------%
%
% TPBVP_dynamics.m provides the odefun for the low thrust controlled s/c time 
% optimal TPBVP.
%
%-------------------------------------------------------------------------%
% PROTOTYPE:
%  dY = TPBVP_dynamics(t,Y,mu,Tmax,Isp,g0)
%
%-------------------------------------------------------------------------%
% INPUT ARGUMENTS:
%  t                [1]     Time instant
%  Y                [14x1]  State and costate array:
%                           - State:
%                             r          [3x1]    Position
%                             v          [3x1]    Velocity
%                             m          [1]      Mass
%                           - Costate (Lagrange multipliers):
%                             lambda_r   [3x1]    Position costate
%                             lambda_v   [3x1]    Velocity costate
%                             lambda_m   [1]      Mass costate
%  Tmax             [1]     Maximum thrust                         [kN]
%  Isp              [1]     Specific impulse of the thruster       [s]
%  g0               [1]     Standard acceleration of free fall     [km/s^2]  
%
%-------------------------------------------------------------------------%
% OUTPUT ARGUMENT:
%  dY               [14x1]  Derivative of state and costate
%
%-------------------------------------------------------------------------%
% AUTHOR:
%  Gian Marco Paldino
%
%-------------------------------------------------------------------------%
% VERSIONS: 
%  01/11/2021: first version
%
%-------------------------------------------------------------------------%

% Extract the state and costate vectors
x = Y(1:7);
lambda = Y(8:14);

% Extract position [km] velocity [km/s] and mass [kg] from the state vector
r = x(1:3);     
v = x(4:6);     
m = x(7);       

% Extract lambda_r, lambda_v, and lambda_m from the costate vector
lambda_r = lambda(1:3);
lambda_v = lambda(4:6);
lambda_m = lambda(7);

% Compute control action
[u,~] = thrust_switch(lambda_v,lambda_m,Isp,m,g0);

% Compute otimal thrust pointing directions
alfa = - lambda_v/norm(lambda_v);

% Define equations of motion for the controlled 2BP
dr = v;
dv = -(mu/norm(r)^3)*r + u*(Tmax/m)*alfa;     
dm = -u * (Tmax/(Isp*g0));

% Derivative of the state
dx = [dr; dv; dm];

% Derivative of the costate
dlambda_r = -3*mu/norm(r)^5*dot(r,lambda_v)*r + mu/norm(r)^3*lambda_v;
dlambda_v = -lambda_r;
dlambda_m = -u*norm(lambda_v)*Tmax/m^2;

dlambda = [dlambda_r; dlambda_v; dlambda_m];

% Full dynamics
dY = [dx; dlambda];

end

%-------------------------------------------------------------------------%    

%-------------------------------------------------------------------------%

% function: thrust_switch.m

function [u,St] = thrust_switch(lambda_v,lambda_m,Isp,m,g0)
%-------------------------------------------------------------------------%    
%
% thrust_switch.m computes the switching function for the time optimal 
% TPBVP for continous low-thrust guidance in the 2BP.
%
%-------------------------------------------------------------------------%
% PROTOTYPE:
%  [u,St] = thrust_switch(lambda_v,lambda_m,Isp,m,g0)
%
%-------------------------------------------------------------------------%
% INPUT ARGUMENTS:
%  lambda_v         [3x1]   Velocity costate
%  lambda_m         [3x1]   Mass costate
%  Isp              [1]     Specific impulse of the thruster      [s]
%  g0               [1]     Standard acceleration of free fall    [km/s^2]  
%
%-------------------------------------------------------------------------%
% OUTPUT ARGUMENTS:
%  u                [1]     Thrust throttle factor (∈ [0,1])      [-]
%  St               [1]     Switching function
%
%-------------------------------------------------------------------------%
% AUTHOR:
%  Gian Marco Paldino
%
%-------------------------------------------------------------------------%
% VERSIONS: 
%  01/11/2021: first version
%
%-------------------------------------------------------------------------%

% Compute the switching function
St = - norm(lambda_v)*Isp*g0/m - norm(lambda_m);

% Compute thrust throttle action depending on the sign of the switching 
% function 
if St >= 0
    u = 0;
else
    u = 1;
end

end

%-------------------------------------------------------------------------%    

%-------------------------------------------------------------------------%

% function: shootingfun_TPBVP.m

function [constraints,time,S] = shootingfun_TPBVP(z,x0,t0,rr_f,vv_f,mu,Tmax,Isp,g0)
%-------------------------------------------------------------------------%    
%
% shootingfun_TPBVP.m provides the shooting function for the time optimal 
% TPBVP for continous low-thrust guidance in the 2BP.
%
%-------------------------------------------------------------------------%
% PROTOTYPE:
%  [constraints,time,S] = shootingfun_TPBVP(z,x0,t0,rr_f,vv_f,mu,Tmax,Isp,g0)
%
%-------------------------------------------------------------------------%
% INPUT ARGUMENTS:
%  z                [8x1]   Array containing the unknowns 
%                           of the problem:
%                           - lambda_r   [3x1]   Position costate
%                           - lambda_v   [3x1]   Velocity costate
%                           - lambda_m   [1]     Mass costate
%                           - tf         [1]     Final time       [d]
%  x0               [7x1]   Initial state of the s/c:
%                           - rr_0       [3x1]   Initial Position [km]
%                           - vv_0       [3x1]   Initial Velocity [km/s]
%                           - m0         [1]     Initial Mass     [kg]
%  t0               [1]     Initial time                          [d]
%  rr_f             [3x1]   Final position                        [km]
%  vv_f             [3x1]   Final veclocity                       [km/s]
%  mu_primary       [1]     Earth standard gravitational          [km^3/s^2]  
%                           constant
%  Tmax             [1]     Maximum thrust                        [kN]
%  Isp              [1]     Specific impulse of the thruster      [s]
%  g0               [1]     Standard acceleration of free fall    [km/s^2]  
%
%-------------------------------------------------------------------------%
% OUTPUT ARGUMENTS:
%  constraints      [8x1]   Constraints of the TPBVP:
%                           - r(tf) - rf = 0        [3x1]
%                           - v(tf) - vf = 0        [3x1]
%                           - lambda_m = 0          [1]
%                           - H(tf) = 0             [1]
%  time             [1xn]   Integration time steps                [d]
%                           (where n is the number of steps)
%  S                [14x1]  Integrated dynamics 
%                           - State:
%                             r          [3x1]    Position        [km]
%                             v          [3x1]    Velocity        [km/s]
%                             m          [1]      Mass            [kg]
%                           - Costate (Lagrange multipliers):
%                             lambda_r   [3x1]    Position costate
%                             lambda_v   [3x1]    Velocity costate
%                             lambda_m   [1]      Mass costate
%
%-------------------------------------------------------------------------%
% AUTHOR:
%  Gian Marco Paldino
%
%-------------------------------------------------------------------------%
% VERSIONS: 
%  01/11/2021: first version
%
%-------------------------------------------------------------------------%

% Retrieve lambda0 and final time
lambda0 = z(1:7);
tf = z(8);    

% Seconds per day
spd = cspice_spd;
% Convert time to seconds
tf = tf*spd;

% Set initial conditions
Y0 = [x0; lambda0];

% Integrate the dynamics
options = odeset('RelTol', 3e-14, 'AbsTol', 1e-14);
[time,S] = ode113(@(t,Y) TPBVP_dynamics(t,Y,mu,Tmax,Isp,g0),[t0 tf],Y0,options);

% Retrieve the solution at final time
YF = S(end,:);
XF = YF(1:7)';
Lambda_f = YF(8:14)';
r = XF(1:3);     
v = XF(4:6);     
m = XF(7);       
lambda_r = Lambda_f(1:3);
lambda_v = Lambda_f(4:6);
lambda_m = Lambda_f(7);

% Compute switching function and thrust throttle factor at final time
[u,St] = thrust_switch(lambda_v,lambda_m,Isp,m,g0);

% Compute the Hamiltonian at final time
Hf = 1 + dot(lambda_r,v) - (mu/(norm(r))^3)*dot(r,lambda_v) + (Tmax/(Isp*g0))*u*St;
       
% Compute constraints 
error_pos = XF(1:3) - rr_f;
error_vel = XF(4:6) - vv_f;
error_lambda_m = lambda_m*(1000);
error_H = Hf;
constraints = [error_pos; error_vel; error_lambda_m; error_H];

end
