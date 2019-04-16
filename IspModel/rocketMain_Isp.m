%% ASEN 2004 - Rocket Bottle Lab - Rocket Equation Model
%{
    Authors: Jashan Chopra (107689146)
    Adapted from 2012 Project 2
    Date Created: March 20th, 2018

[ Rocket Equation Model - Model 2 ]

Script Purposes and goals:
    1) Import static test stand data
    2) Integrate test stand data to get Isp
    3) Solve for DeltaV
    4) Solve for final distance

Outputs:
    1) Final Distance
    2) Max Height
    3) Total Time
    4) Angle Deviation from launch
%}

% General Housekeeping
clc; clear; close all;

%% Define Initial Variables

g = 9.80665;                % gravitational constant
mProp = 1;                  % [kg] [1 L = 1000g = 1kg]
mBottle=130/1000;           %[kg] % Mass of Empty Bottle
Patm=83427;                 %[Pa] % Atmospheric Pressure
Pgage=275790;               %[Pa] % Intial Gauge Pressure of air in bottle

volBottle=.002;             %[m^3] % Empty bottle Volume
volWater=.001;              %[m^3] % Initial Volume of Water
volAir=volBottle-volWater;  % Volume of air in Bottle Initial

R=287;          % Universal Gas Constant Air
Tair=290.15;    %[K] % Initital Temp of Air

mAir=((Pgage+Patm)*volAir)/(R*Tair);    % Mass of Air Initital
mAirFinal=((Patm)*volAir)/(R*Tair);     % Mass of Air Final


%% Import Test Stand Data
% Import data from the rocket test stand
    f = 'LA8am_test3';
    data = fileLoad(f);                                                 % Load data

    % Shift up negative values
    negData = data < 0;                                                 % negative data values
    low = mean(data(negData));                                          % mean of negative values
    low = 4*low;                                                        % 4x correction factor provides best estimate
    data = data + abs(low);                                             % add that mean to shift values

    % Remove extraneous negative
    indicies = find(data <= 0);                                         % negative indices
    data(indicies) = [];                                                % remove negative indices

    % curve fit
    frequency = 1.652 * 1000;                                           % [Hz] Sampling Rate
    time = (1 / frequency) * linspace(0,length(data),length(data))';    % time vector
    fitobject = fit(time,data,'smoothingspline');                       % cubic interp

    % extraneous value removal
    fx = abs(differentiate(fitobject, time));                           % calculate slope of various points
    deletion = find(fx <= 600);                                         % deletion parameter
    data(deletion) = []; time(deletion) = [];                           % remove values from data
    time = time - time(1);                                              % reset time to 0
    fitobject = fit(time,data,'cubicinterp');                           % refit data
%% Ploting

    % Plot thrust
    figure(1)
    plot(fitobject,time,data)
    title('Force over time')
    xlabel('Time [S]')
    ylabel('Force [N]')
    legend('8am Test 3','Fitted Interpolation')

        %% BOTH INTEGRATION METHODS RETURN SAME VALUE --
        % Integrate to find area under curve
        isp = trapz(time,data);
        isp = isp / (mProp*g);

        % alternative method to verify ISP
        isp2 = integrate(fitobject,time(end),time(1));
        isp2 = isp2 / (mProp*g);

%% Calculations

% Use Isp to calculate deltaV
    mInitial = mBottle + mProp;         % Initial Mass
    mFinal = nBottle;                   % Final Mass

    deltaV = isp*g*log(mInitial / mFinal); % [m/s] % DeltaV

%% Call ode45 using initial instantaneous deltaV

% Wind angles
launchTheta = 204; launchTheta = launchTheta * (pi/180);    % launch angle from North
sufWindTh = 67.5; surWind = sufWindTh * (pi/180);           % surface wind from North
alofWindTh = 22.5; alofWindTh = alofWindTh * (pi/180);      % aloft wind from North

% Wind magnitudes
sufWind = 3; sufWind = sufWind * 0.44704;                   % convert to m/s wind
alofWind = 3; alofWind = alofWind * 0.44704;                % convert to m/s wind

% calculate relative downrange and crossrange wind components
theta3 = launchTheta - alofWindTh - (90*pi/180);            % angles to crossrange & downrange
theta4 = launchTheta - alofWindTh;
alofY = alofWind*cos(theta3);                               % dot products
alofX = -alofWind*cos(theta4);

theta3 = launchTheta + sufWindTh - (90*pi/180);             % angles to crossrange & downrange
theta4 = launchTheta + sufWindTh;
sufY = sufWind*cos(theta3);                                 % dot products
sufX = -sufWind*cos(theta4);

vwx = sufX + alofX; vwy = sufY + alofY;                     % wind added

vw = [vwx;vwy;0]; % Initial velocity of the wind

theta = pi/4; % Initial Angle

    vx = deltaV*cos(theta); % Initial Velocity x
    vy = 0;                 % Initial Velocity y
    vz = deltaV*sin(theta); % Initial Velocity z

    x0 = 0;     %[m] % Initial Horizontal Distance
    z0 = .25;   %[m] % Initial Vertical Distance
    y0 = 0;     %[m] % Initial Lateral Distance

% Initialize the state vector
stateVec = [vx;vz;vy;x0;z0;y0];

% Manual timestep helps ode45 accuracy
tspan = 0:.00001:5; % Create a time span, 0 to 2 seconds

% Odeset triggers a cancel event in ode45
optstage1 = odeset('Events',@opt1); % opt1.m outputs the cancel signal

% Call to ode45 with state1func,state1vec,tspan, and optstage1
[t, ds] = ode45(@(t,s) state2funcROCKET(t,s,vw), tspan, stateVec, optstage1);

%% Plotting

figure(2) % Plot the distance

% Remove values where z has gone below the surface of the Earth
indices = find(ds(:,5) < 0);
ds(indices,:) = [];

% Display Max Height
disp('Max Height [m]')
disp(max(ds(:,5)));

% Display Distance
disp('Final Distance [m]')
disp(ds(end,4));

% Display Final Time
disp('Final Time [s]')
disp(t(indices(1)))

% Display Final Beta Value
betaFinal = atan((ds(length(ds),6) / ds(length(ds),4)));
betaFinal = (betaFinal * 180) / pi;
disp('Lateral Angle Displament [deg]')
disp(abs(betaFinal))

% Display Summed Distance
sum = sqrt(ds(end,4)^2 + abs(ds(end,6))^2);
disp('Total Distance [m]')
disp(sum)

% Point of Max Height and point of landing
[maxZ, I] = max(ds(:,5));      % Max height
thisX = ds(I,4);               % Indicies at max height
thisY = ds(I,6);

maxX = max(ds(:,4));           % Final distances
maxY = max(ds(:,6));

p = plot3(ds(:,4),ds(:,6),ds(:,5)); % Plot height v. distance
p.LineWidth = 2;
hold on
plot3(thisX,thisY,maxZ,'o')     % Point at max height
plot3(maxX,maxY,0,'o')          % Point at landing
title('Height v. Distance')     % Title
xlabel('Downrange [m]')         % X label
ylabel('Crossrange [m]')        % Y label
zlabel('Height')                % Z Label
ylim([-5 5]); xlim auto; zlim auto;                 % scaling
legend('Trajectory','Max Height','Final Location')  % legend
set(gca,'Ydir','reverse')                           % flip y axis

%% Checks & Comments
    % higher initial drag due to larger velocity, less inertia later
    % greater initial angle of trajectory due to large initial velocity
    % Overall, these two factors lead to less distance
