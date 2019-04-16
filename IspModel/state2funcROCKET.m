function [ds2] = state2funcROCKET(t,stateVec,vw)
% This function serves as the main input for ode45 for stage two and stage
% three
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Authors: Jason Balke (107025894), Jashan Chopra (107689146)
% Date Created: Novemeber 20th, 2018
% Last Date Modified: November 30th, 2018

%% Initial conditions

g=9.80665;              %[m/s] % Gravity Constant
rhoAtm=.961;            %[kg/m3] % Ambient Air Density
Dbottle= 107.5/1000;    %[m] % Diameter of Bottle
mBottle=130/1000;       %[kg] % Mass of Empty Bottle
CDrag=.3773;            % Drag Coefficient

% Calculated Constants
Abottle=pi*(Dbottle/2)^2;   % Area of Bottle

vx=stateVec(1);    % X Vel
vz=stateVec(2);    % Z Vel
vy=stateVec(3);    % Y Vel
x=stateVec(4);     % X Dist
z=stateVec(5);     % Z Dist
y=stateVec(6);     % Y Dist

%% Calculations

    vRel = [vx + vw(1), vy + vw(2), vz + vw(3)];        % add wind components
    velocity = norm(vRel);                              % Vector magnitude
    Fdrag = (1/2)*rhoAtm*(velocity^2)*CDrag*Abottle;    % Drag Force

    hypo = sqrt(vRel(1)^2+vRel(2)^2);   % hypotenuse of flight path xy proj
    theta = atan(vRel(3)/ hypo);        % The angle value at any point
    beta = atan(vRel(2)/ hypo);         % Lateral angle value at any point

    dragx = Fdrag*cos(theta)*cos(beta);          % x-component of drag
    dragz = Fdrag*sin(theta);                    % z-component of drag
    dragy = Fdrag*cos(theta)*sin(beta);          % y-component of drag

    mSystem = mBottle;      % Mass of the system for use in F=ma
    Fgravity = mSystem*g;   % Force from Gravity

    dvx= -dragx / mSystem;              % Acceleration in x --> integrates to vel x
    dvz= (-dragz-Fgravity) / mSystem;   % Same as above, in z direction
    dvy = dragy / mSystem;             % acceleration in y direction

    dx = vx; % velocity x --> integrates to distance x
    dy = vy; % same as above, in y direction
    dz = vz; % same as above, in z direction

    %% Output

    ds2=[dvx;dvz;dvy;dx;dz;dy]; % output call for ode45

end
