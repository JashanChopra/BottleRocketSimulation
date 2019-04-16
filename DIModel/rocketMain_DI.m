%% ASEN 2004 - Rocket Bottle Lab - Direct Interpolation Model
%{
    Authors: Jashan Chopra (107689146)
    Adapted from 2012 Project 2
    Date Created: March 21st, 2018

[ Direct Interpolation - Model 3 ]

Script Purposes and goals:
    1) Import static test stand data
    2) Adjust ode45 state functions to use interpolated force
    3) Rearrange thermodynamic models using interpolated force
    4) Solve for final distance


Outputs:
    1) Final Distance
    2) Max Height
    3) Total Time
    4) Angle Deviation from launch
%}

% General Housekeeping
clc; clear; close all;

%% Initial Variables

% Given Constants
g=9.80665;              %[m/s] % Gravity Constant
Cd=.8;                  % Discharge Coefficient
rhoAtm=.961;            %[kg/m3] % Ambient Air Density
rhoWater=1000;          %[kg/m3] % Density of Water
volBottle=.002;         %[m^3] % Empty bottle Volume
Patm=83427;             %[Pa]   % Atmospheric Pressure
gamma=1.4;              % Ratio of Specific Heats for Air
Dthroat=.021;           %[m] % Diameter of Throat
Dbottle=.105;           %[m] % Diameter of Bottle
R=287;                  % Universal Gas Constant Air
mBottle=117/1000;       %[kg] % Mass of Empty Bottle
CDrag=.3361;            % Drag Coefficient
Pgage=275790;           %[Pa] % Intial Gauge Pressure of air in bottle
volWater=0.000962;      %[m^3] % Initial Volume of Water
Tair=275.372;           %[K] % Initital Temp of Air
v0=0;                   % [m/s] Initial Velocity
theta=pi/4;             %[rad] % Initial Angle
ls=.5;                  %[m] % Length of Test Stand

% Calculated Constants
volAir=volBottle-volWater;              % Volume of air in Bottle Initial
mAir=((Pgage+Patm)*volAir)/(R*Tair);    % Mass of Air Initital
Athroat=pi*(Dthroat/2)^2;               % Area of Throat
Abottle=pi*(Dbottle/2)^2;               % Area of Bottle
mWater = rhoWater * volWater;           % Initial mass of water

x0=0;   %[m] % Initial Horizontal Distance
z0=.25; %[m] % Initial Vertical Distance
y0 = 0; %[m] % Initial Lateral Distance

% Wind angles
launchTheta = 204; launchTheta = launchTheta * (pi/180);    % launch angle from North
sufWindTh = 67.5; surWind = sufWindTh * (pi/180);           % surface wind from North
alofWindTh = 22.5; alofWindTh = alofWindTh * (pi/180);      % aloft wind from North

% Wind magnitudes
sufWind = 1; sufWind = sufWind * 0.44704;                   % convert to m/s wind
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

vx = 0; % Initial Velocity x
vy = 0; % Initial Velocity y
vz = 0; % Initial Velocity z

%% Import Test Stand Data
% Import data from the rocket test stand
    f = 'LA8am_test3';
    data = fileLoad(f);                                                 % Load data

    % Shift up negative values
    negData = data < 0;                                                 % negative data values
    low = mean(data(negData));                                          % mean of negative values
    data = data + abs(low);                                             % add that mean to shift values

    % Remove extraneous negative
    indicies = data <= 0;                                               % negative indices
    data(indicies) = [];                                                % remove negative indices

    % curve fit
    frequency = 1.652 * 1000;                                           % [Hz] Sampling Rate
    time = (1 / frequency) * linspace(0,length(data),length(data))';    % time vector
    fitobject = fit(time,data,'smoothingspline');                       % cubic interp

    % extraneous value removal
    fx = abs(differentiate(fitobject, time));                           % calculate slope of various points
    deletion = find(fx <= 1200);                                        % deletion parameter
    data(deletion) = []; time(deletion) = [];                           % remove values from data
    time = time - time(1);                                              % reset time to 0
    fitobject = fit(time,data,'cubicinterp');                           % refit data

%% Plot static test stand data

    % Plot thrust
    figure(1)
    plot(fitobject,time,data) % Plot the fitted line with actual data points
    hold on
    title('Force over time')
    xlabel('Time [S]')
    ylabel('Force [N]')
    legend('8am Test 3','Fitted Interpolation')

%% Phase 1 - Intepolation of the thrust as the static test data

% First we create a state vector to put into ode45
state1Vec = [volAir vx vz vy x0 z0 y0];

% Manual timestep helps ode45 accuracy
tspan = 0:.00001:.5; % Create a time span, 0 to .5 seconds

% Odeset triggers a cancel event in ode45
optstage1 = odeset('Events',@opt1); % opt1.m outputs the cancel signal

% Call to ode45 with state1func,state1vec,tspan, and optstage1
[t, ds] = ode45(@(t,s) state1funcDI(t,s,Cd,Athroat,gamma,...
    rhoAtm,Patm,Pgage,volAir,Abottle,mBottle,mAir,CDrag,volBottle,...
    rhoWater,vw,fitobject), tspan, state1Vec, optstage1);

%% Phase 2 - Intepolation of the thrust as the static test data

m2 = mAir;
Pair = Pgage + Patm; % Defines pressure of air inside bottle

% Initialize the second state vector based off previous conditions
% [ mass_air x-vel z-vel y-vel dist-x dist-z dist-y]
state2Vec = [m2;ds(end,2);ds(end,3);ds(end,4);ds(end,5);ds(end,6);ds(end,7)];

% Time span for ode45
tspan2 = t(end):.000001:6;

% Call to ode45 for the second stage, similar to before
optstage2 = odeset('Events',@opt2); % This function modifies ode45 to stop
[t2,ds2] = ode45(@(t2,s) state2funcDI(t2,s,mAir,gamma,Pair,volBottle,R,Patm,...
    Athroat,Cd,Abottle,rhoAtm,volAir,CDrag,mBottle,vw,fitobject),tspan2, state2Vec,optstage2);

%% Plot Distances
figure(2) % Opens second figure

% Remove values where z has gone below the surface of the Earth
indices = find(ds2(:,6) < 0);
ds2(indices,:) = [];

% Concatenate two stages
xpos=[ds(:,5);ds2(:,5)];
zpos=[ds(:,6);ds2(:,6)];
ypos=[ds(:,7);ds2(:,7)];

% Display Max Height
disp('Max Height [m]')
disp(max(ds2(:,6)));

% Display Distance
disp('Final Distance [m]')
disp(ds2(end,5));

% Display Final Time
disp('Final Time [s]')
if ~isempty(indices)
    disp(t2(indices(1)))
else
    disp('Model is messed up')
end

% Display Final Beta Value
betaFinal = atan((ds2(length(ds2),7) / ds2(length(ds2),5)));
betaFinal = (betaFinal * 180) / pi;
disp('Lateral Angle Displament [deg]')
disp(abs(betaFinal))

% Display Summed Distance
sum = sqrt(ds2(end,5)^2 + abs(ds2(end,7))^2);
disp('Total Distance')
disp(sum)

% Point of Max Height and point of landing
[maxZ, I] = max(ds2(:,6));      % Max height
thisX = ds2(I,5);               % Indicies at max height
thisY = ds2(I,7);

maxX = max(ds2(:,5));           % Final distances
if ds2(end,7) > 0
    maxY = max(ds2(:,7));       % if statement controls max/min (+/-)
else
    maxY = min(ds2(:,7));
end

p = plot3(xpos,ypos,zpos);      % Plot height v. distance
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

%% Appendix 1 - Intepolating Thrust at any given timestep
% ! This method does not work for the given scenario ! %
%{
    This following code will be IN the ode45 state function

    % At any point we have the force vector 'force1'
        % Inside ode45 we have the integer time value, 't'
        % Also bring in 'tspan' as an input to ode45 along with 'force1'

    1) find the index of tspan that equals current t value
        index = find(tspan == t);

    2) Take the index value and divide it by the ratio of force values in tspan
        % divsor = length(tspan) / length(force1);

    3) Round up the divided index value --> this is the index of force1 required
        % index = ceil(index / divsor)

    4) Index 'force1' at that value for usage in rest of the equations
        forcethrust = force1(index)

   % Combine into one nice line
   index = ceil((find(tspan == t)) / (length(tspan) / length(force1)));

%}
