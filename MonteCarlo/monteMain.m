%% ASEN 2004 - Rocket Bottle Lab - Monte Carlo Results & Plots
%{

    Authors: Jashan Chopra (107689146)
    Date Created: April 3rd, 2019

Script Purposes and goals:
    1) Call monte.m to get monte carlo simulation data
    2) Display mean and STD of each variable for each method
    3) Display conical error plots with x/y data values overlayed
    4) Combine data from two models to compare error bounds and overlay
    5) Display all plots and data neatly

acknowledgements:
    1) https://www.xarg.org/2018/04/how-to-plot-a-covariance-error-ellipse/

%}
tic
% General Housekeeping
clc; clear; close all;

%% Varying Parameters
numTrials = 250;        % number of trials
windMult = 1.5;           % Wind multiplier

% Wind Values [Randomized - normal]
vwx = windMult*abs(randn(numTrials,1));     % Normal distribution, only negative values
vwy = windMult*abs(randn(numTrials,1));
vw = [vwx,vwy,zeros(numTrials,1)];          % Randomized wind vector

% Initial Volume Water [Randomized - normal]
    % variation directly in monte.m

% Initial mass of bottle [Randomized - normal]
    % standard deviation is .01 kg, 10g -- Too much?

% Coefficient of drag [Randomized - normal]
    % normal value is .5, .01 is acceptable variation from dirt, etc?

% Tempature of Air [Randomized - normal]
    % normal value is 300K, 1 degree kelvin variation, maybe more?

%% Call to monte carlo simulation function
monteMat = monte(numTrials,vw);    % call monte.m

%% Plot Error Ellipses [Thermo Data]
prob = [.99,.95,.9];    % probability values at 99%, 95%, 90%
for p = prob            % Iterates for each error ellipse wanted

    % covariance is defined as [ sigmaX^2 , sigmaXY  ]
    %                          [ sigmaYX  , sigmaY^2 ]

    s = -2 * log(1 - p);                        % calculate the Mahalanobis radius
    sigma = cov(monteMat(:,1),monteMat(:,2));   % calculate the sigma covariance
    [V,D] = eig(s * sigma);                     % calculate eiganvals of weighted covariance
    theta = linspace(0,2*pi,360);               % theta linspace vector for ellipse

    radii = (V * sqrt(D)) * [cos(theta(:))' ; sin(theta(:))'];  % ellipse functions
    x = radii(1,:) + mean(monteMat(:,1));                       % x component
    y = radii(2,:) + mean(monteMat(:,2));                       % y component

    figure(1)
    e1 = plot(x,y);   % plotting
    hold on
end

%% Plot Error Ellipses [Rocket Data]
prob = [.99,.95,.9];    % probability values at 99%, 95%, 90%
for p = prob            % Iterates for each error ellipse wanted

    s = -2 * log(1 - p);                        % calculate the Mahalanobis radius
    sigma = cov(monteMat(:,4),monteMat(:,5));   % calculate the sigma covariance
    [V,D] = eig(s * sigma);                     % calculate eiganvals of weighted covariance
    theta = linspace(0,2*pi,360);               % theta linspace vector for ellipse

    radii = (V * sqrt(D)) * [cos(theta(:))' ; sin(theta(:))'];  % ellipse functions
    x = radii(1,:) + mean(monteMat(:,4));                       % x component
    y = radii(2,:) + mean(monteMat(:,5));                       % y component

    figure(1)
    e2 = plot(x,y);   % plotting
    hold on
end
%{
%% Plot Error Ellipses [DI Data]
prob = [.99,.95,.9];    % probability values at 99%, 95%, 90%
for p = prob            % Iterates for each error ellipse wanted

    s = -2 * log(1 - p);                        % calculate the Mahalanobis radius
    sigma = cov(monteMat(:,7),monteMat(:,8));   % calculate the sigma covariance
    [V,D] = eig(s * sigma);                     % calculate eiganvals of weighted covariance
    theta = linspace(0,2*pi,360);               % theta linspace vector for ellipse

    radii = (V * sqrt(D)) * [cos(theta(:))' ; sin(theta(:))'];  % ellipse functions
    x = radii(1,:) + mean(monteMat(:,7));                       % x component
    y = radii(2,:) + mean(monteMat(:,8));                       % y component

    figure(1)
    e3 = plot(x,y);   % plotting
    hold on
end
%}
%% Final Plots

figure(1)
plot(mean(monteMat(:,1)),mean(monteMat(:,2)),'*','MarkerSize',10);    % Thermo mean
plot(mean(monteMat(:,4)),mean(monteMat(:,5)),'*','MarkerSize',10);    % Rocket mean
% plot(mean(monteMat(:,7)),mean(monteMat(:,8)),'*','MarkerSize',10);    % DI mean
p3 = scatter(monteMat(:,1),monteMat(:,2));          % All thermo data
p4 = scatter(monteMat(:,4),monteMat(:,5));          % All rocket data
% p5 = scatter(monteMat(:,7),monteMat(:,8));        % All DI data
alpha(p3,.1); alpha(p4,.1); % alpha(p5,.1);         % set transparency
ylim([-10 20]); xlim([0 100]);                        % scaling

legend('99% T','95% T','90% T',...
    '99% Isp','95% Isp','90% Isp',...
    'Thermo Avg','Isp Avg',...
    'Thermo Data','Isp Data')

title('Thermodynamic & Rocket Model Error Ellipses')
xlabel('Downrange Distance [m]')
ylabel('Crossrange Distance [m]')

toc
