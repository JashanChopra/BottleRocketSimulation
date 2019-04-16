%% ASEN 2004 - Rocket Bottle Lab - Main Code
%{
    Authors: Jashan Chopra (107689146), Dom Dougherty (106073184)
    Adapted from 2012 Project 2
    Date Created: March 20th, 2018

[ THERMODYNAMIC MODEL - Model 1 ]

Script Purposes and goals:
    1) Combine best aspects of Jashan & Sasha's original 2012 code
        Specifically fancy plotting features
    2) Adjust code to incorporate the heading velocity vector
    3) Adjust code to incorporate affects of wind in heading velocity
    4) Sensitivity Analysis on multiple variables to determine improvments
    5) Convert the plots to 3D to predict exact trajectory of launch

Displays:
    1) Final Distance
    2) Max Height
    3) Total Time
    4) Angle Deviation from launch
    5) Total Distance
%}

% General Housekeeping
clc; clear; close all;

%% Identify Constants

% Given Constants
g=9.80665;              %[m/s] % Gravity Constant
Cd=.8;                  % Discharge Coefficient
rhoAtm=.961;            %[kg/m3] % Ambient Air Density
rhoWater=1000;          %[kg/m3] % Density of Water
volBottle=.002;         %[m^3] % Empty bottle Volume
Patm=83427;             %[Pa]   % Atmospheric Pressure
gamma=1.4;              % Ratio of Specific Heats for Air
Dthroat=21.32/1000;     %[m] % Diameter of Throat
Dbottle=107.5 / 1000;   %[m] % Diameter of Bottle
R=287;                  % Universal Gas Constant Air
mBottle=130/1000;     %[kg] % Mass of Empty Bottle
CDrag=.3773;            % Drag Coefficient
Pgage=275790;           %[Pa] % Intial Gauge Pressure of air in bottle
volWater=0.001;         %[m^3] % Initial Volume of Water
Tair=290.15;            %[K] % Initital Temp of Air
v0=0;                   % [m/s] Initial Velocity
theta= 40 * pi / 180;   %[rad] % Initial Angle
ls=.5;                  %[m] % Length of Test Stand

x0=0;   %[m] % Initial Horizontal Distance
z0=.25; %[m] % Initial Vertical Distance
y0 = 0; %[m] % Initial Lateral Distance

% Calculated Constants

volAir=volBottle-volWater;              % Volume of air in Bottle Initial
mAir=((Pgage+Patm)*volAir)/(R*Tair);    % Mass of Air Initital
Athroat=pi*(Dthroat/2)^2;               % Area of Throat
Abottle=pi*(Dbottle/2)^2;               % Area of Bottle

mWater = rhoWater * volWater; % Initial mass of water

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

vx = 0; % Initial Velocity x
vy = 0; % Initial Velocity y
vz = 0; % Initial Velocity z

%% Phase 1

% First we create a state vector to put into ode45
% [Volume Air, Velocities(x,y,z), Distances(x,y,z)
state1Vec = [volAir vx vz vy x0 z0 y0];

% Manual timestep helps ode45 accuracy
tspan = 0:.00001:.5; % Create a time span, 0 to .5 seconds

% Odeset triggers a cancel event in ode45
optstage1 = odeset('Events',@opt1); % opt1.m outputs the cancel signal

% Call to ode45 with state1func,state1vec,tspan, and optstage1
[t, ds] = ode45(@(t,s) state1func(t,s,Cd,Athroat,gamma,...
    rhoAtm,Patm,Pgage,volAir,Abottle,mBottle,mAir,CDrag,volBottle,...
    rhoWater,vw), tspan, state1Vec, optstage1);

% Initialize a vector of thrust values
Fthrust=zeros(1,length(t));

% For each value in the first time vector, we run the state func with our
% full state matrix to get the thrust
for i=1:length(t)
    in1=t(i);       % Time value at any time
    in2=ds(i,:);    % Row of the state matrix for time
    [~,Fthrust(i)] = state1func(in1,in2,Cd,Athroat,gamma,rhoAtm,Patm,Pgage,volAir,...
    Abottle,mBottle,mAir,CDrag,volBottle,rhoWater,vw);
end

%% Phase 2

% Define a new variable just to avoid confusion with initial mass of air
m2 = mAir;
Pair = Pgage + Patm; % Defines pressure of air inside bottle

% Initialize the second state vector based off previous conditions
% [ mass_air x-vel z-vel y-vel dist-x dist-z dist-y]
state2Vec = [m2;ds(end,2);ds(end,3);ds(end,4);ds(end,5);ds(end,6);ds(end,7)];

% Time span for ode45
tspan2 = t(end):.000001:6;

% Call to ode45 for the second stage, similar to before
optstage2 = odeset('Events',@opt2); % This function modifies ode45 to stop
[t2,ds2] = ode45(@(t2,s) state2func(t2,s,mAir,gamma,Pair,volBottle,R,Patm,...
    Athroat,Cd,Abottle,rhoAtm,volAir,CDrag,mBottle,vw),tspan2, state2Vec,optstage2);

% Exact same process as before, but for stage 2
Fthrust2=zeros(1,length(t2));
beta = zeros(1,length(t2)); % Also gets beta this time around

for i=1:length(t2)
    in3=t2(i);
    in4=ds2(i,:);
    [~,Fthrust2(i),beta(i)] = state2func(in3,in4,mAir,gamma,Pair,volBottle,R,Patm,...
    Athroat,Cd,Abottle,rhoAtm,volAir,CDrag,mBottle,vw);
end

%% Plot Thrust
figure(1) % Opens the first figure
% Concatenate two stages
timebot = [t;t2];
thrusttotal = [Fthrust Fthrust2];

% Plot the time v. thrust
p = plot(timebot,thrusttotal);
p.LineWidth = 2;
title('Thrust Diagram') % Title
xlabel('Time [Seconds]') % X label
ylabel('Thrust [N]') % Y label
xlim([0 .45]) % Shortens the x axis to avoid misc 0 values of thrust

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
disp(t2(indices(1)))

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

% Double check to make sure displacement laterally is parabolic
%{
    figure(6)
    plot(ds2(:,5),ds2(:,7))
    hold on
    plot([ds2(1,5) ds2(end,5)],[ds2(1,7) ds2(end,7)])
    title('Check')
    xlabel('X Displacement')
    ylabel('Y Displacement')
%}
