function deltaV = dataRead()
%% ASEN 2004 - Rocket Bottle Lab
%{

    Authors: Jashan Chopra (107689146)
    Date Created: April 4th, 2019

Script Purposes and goals:
  1) Read static test stand data
  2) edit data to appropriate ranges
  3) automate the data curtailing process
  4) output data to matrix for statistical usage

%}

% General Housekeeping
clc; clear; close all;

%% Import the test stand data, define initial conditions
g = 9.80665;        % gravitational constant
mProp = 1;          % [kg] [1 L = 1000g = 1kg]
mBottle=.15;        %[kg] % Mass of Empty Bottle
Patm=83427;         %[Pa] % Atmospheric Pressure
Pgage=345738;       %[Pa] % Intial Gauge Pressure of air in bottle

volBottle=.002;             %[m^3] % Empty bottle Volume
volWater=.001;              %[m^3] % Initial Volume of Water
volAir=volBottle-volWater;  % Volume of air in Bottle Initial

R=287;      % Universal Gas Constant Air
Tair=300;   %[K] % Initital Temp of Air

mAir=((Pgage+Patm)*volAir)/(R*Tair);    % Mass of Air Initital
mAirFinal=((Patm)*volAir)/(R*Tair);     % Mass of Air Final

frequency = 1.652 * 1000;               % [Hz] Sampling Rate

 filenames = ["Goup03_10am_statictest1","group04_10am_statictest1","Group01_01pm_statictest1",...
"Group01_08am_statictest1","Group01_10am_statictest1","Group02_08am_statictest1",...
"Group03_1pm_statictest1","Group03_8am_statictest1","Group3_8AM_statictest2",...
"Group04_01pm_statictest1","Group04_08am_statictest1","Group05_08am_statictest1",...
"Group06_8am_statictest1","Group06_10am_statictest1",...
"Group07_01pm_statictest500g","Group07_01pm_statictest750g","Group07_01pm_statictest1250g",...
"Group07_08am_statictest1","Group07_10am_statictest1","Group07_10am_statictest2",...
"Group07_10am_statictest500g","Group07_10am_statictest600g","Group07_10am_statictest700g",...
"Group08_08AM_statictest1","Group08_8AM_statictest500g","Group08_8AM_statictest750g",...
"Group08_10AM_statictest2","Group08_10AM_statictest4","Group09_1PM_statictest1",...
"Group09_08am_statictest1","Group09_10Am_statictest2","Group09_10AM_statictest3",...
"Group09_10AM_statictest4","Group10_01pm_statictest1","Group10_1PM_statictest600g",...
"Group10_8am_statictest1.csv","Group10_8am_statictest600","Group10_8am_statictest600mL2",...
"Group11_01PM_statictest640g","Group11_08AM_statictest1.txt"];  % Create a vector of all filenames

numFiles = length(filenames);                 % Number of files

%% Iteration of deltaV calculation
deltaV = zeros(numFiles,1);     % preallocate deltaV
for f = filenames               % loop iterates over all files

    % automated data trim
    data = fileLoad(f);                                                % Load data
    indicies = find(data <= 0);                                        % negative indices
    data(indicies) = [];                                               % remove negative indices

    % curve fit
    time = (1 / frequency) * linspace(0,length(data),length(data))';    % time vector
    fitobject = fit(time,data,'smoothingspline');                       % cubic interp

    % extraneous value removal
    fx = abs(differentiate(fitobject, time));                           % calculate slope of various points
    deletion = find(fx <= 900);                                         % deletion parameter
    data(deletion) = []; time(deletion) = [];                           % remove values from data
    time = time - time(1);                                              % reset time to 0
    fitobject = fit(time,data,'cubicinterp');                           % refit data

    % calculations
    area = integrate(fitobject,time(end),time(1));                      % integration
    isp = area / (mProp*g);                                             % isp calc

    mInitial = mBottle + mProp + mAir;                                  % initial mass
    mFinal = mInitial - mProp - mAir + mAirFinal;                       % final mass

    deltaV(find(filenames==f),:) = isp*g*log(mInitial / mFinal);        % [m/s] % DeltaV

    %% optional graphing block
    figure(1)
    plot(fitobject,time,data)
    hold on
    title('Force over time')
    xlabel('Time [S]')
    ylabel('Force [N]')
    legend('Data Set','Fitted Interpolation')

end
end
