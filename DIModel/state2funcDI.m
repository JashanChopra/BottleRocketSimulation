function [ds2,Fthrust2,beta] = state2funcDI(t2,state2Vec,mAir,gamma,Pair,volBottle,R,Patm,Athroat,Cd,Abottle,rhoAtm,volAir,CDrag,mBottle,vw,fitobject)
% This function serves as the main input for ode45 for stage two and stage
% three
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Authors: Jason Balke (107025894), Jashan Chopra (107689146)
% Date Created: Novemeber 20th, 2018s
% Last Date Modified: November 30th, 2018

%% Initial conditions
g = 9.80665; %[m / s^2]
Pgage=275790;

m=state2Vec(1); % Mass of air
vx=state2Vec(2); % X Vel
vz=state2Vec(3); % Z Vel
vy=state2Vec(4); % Y Vel
x=state2Vec(5); % X Dist
z=state2Vec(6); % Z Dist
y=state2Vec(7); % Y Dist

%% Calculations
Pend = Pgage * ((volAir / volBottle)^gamma); % Pressure at the end

% m is the actual mass of air, mAir is the initial mass of air
P = ((m/mAir)^gamma)*Pend; % Calculates the pressure at any point

if P > Patm
    rhoBottle=m/volBottle; % Density of air
    TBottle=P/(rhoBottle*R); % Temperature inside bottle
    Pcrit=P*(2/(gamma+1))^(gamma/(gamma-1)); % Calculates the critical pressure

    % If the critical pressure is greater, we are not choked
    if Pcrit>Patm
        Texit=(2/(gamma+1))*TBottle; % Calculates the exit tempature
        rhoExit=Pcrit/(R*Texit); % Calculates the exit density
        Ve=sqrt(gamma*R*Texit); % Calculates the exit velocity
    % Else, we have choked flow and our equations change
    else
        M=sqrt((((P/Patm)^((gamma-1)/gamma))-1)/((gamma-1)/2)); % Mach Number
        Texit=TBottle/(1+((gamma-1)/2)*M^2); % Exit temp
        rhoExit=Patm/(R*Texit); % Exit density
        Ve=M*sqrt(gamma*R*Texit); % Exit velocity
    end

    mdot=-Cd*rhoExit*Athroat*Ve; % The mass flow, to be integrated
    Fthrust2=-mdot*Ve+(P-Patm)*Athroat; % Force of Thrust

    % Drag calculation from relative velocity
    vRel = [vx + vw(1), vy + vw(2), vz + vw(3)]; % add wind components
    velocity = norm(vRel); % Find the vector magnitude
    Fdrag = (1/2)*rhoAtm*(velocity^2)*CDrag*Abottle; % Drag Force Equation

    hypo = sqrt(vRel(1)^2+vRel(2)^2);   % hypotenuse of flight path xy proj
    theta = atan(vRel(3)/ hypo);        % The angle value at any point
    beta = atan(vRel(2)/ hypo);         % Lateral angle value at any point

    thrustx = Fthrust2*cos(theta)*cos(beta);      % x-component of Thrust
    thrustz = Fthrust2*sin(theta);                % z-component of Thrust
    thrusty = Fthrust2*cos(theta)*sin(beta);      % y-component of Thrust

    dragx = Fdrag*cos(theta)*cos(beta);          % x-component of drag
    dragz = Fdrag*sin(theta);                    % z-component of drag
    dragy = Fdrag*cos(theta)*sin(beta);          % y-component of drag

    mSystem = mBottle + m; % Mass of the system for use in F=ma
    Fgravity = mSystem*g; % Force from Gravity

    dvx=(thrustx-dragx)/mSystem; % Acceleration in x --> integrates to vel x
    dvz=(thrustz-dragz-Fgravity)/mSystem; % Same as above, in z direction
    dvy = (thrusty-dragy)/mSystem; % acceleration in y direction

    % If the velocities are negative make the distance 0 so it does not
    % appear to go through the test stand
    % Now that it's falling, negatives are fine

    dx = vx; % velocity x --> integrates to distance x
    dy = vy; % same as above, in y direction
    dz = vz; % same as above, in z direction

    %% Output

    ds2=[mdot;dvx;dvz;dvy;dx;dz;dy]; % output call for ode45
else
    Fthrust2 = 0; % Thrust is zero now
    mdot = 0; % Change in mass is zero

    % Drag calculation from relative velocity
    vRel = [vx + vw(1), vy + vw(2), vz + vw(3)]; % add wind components
    velocity = norm(vRel); % Find the vector magnitude
    Fdrag = (1/2)*rhoAtm*(velocity^2)*CDrag*Abottle; % Drag Force Equation

    hypo = sqrt(vRel(1)^2+vRel(2)^2);   % hypotenuse of flight path xy proj
    theta = atan(vRel(3)/ hypo);        % The angle value at any point
    beta = atan(vRel(2)/ hypo);         % Lateral angle value at any point

    thrustx = Fthrust2*cos(theta)*cos(beta);      % x-component of Thrust
    thrustz = Fthrust2*sin(theta);                % z-component of Thrust
    thrusty = Fthrust2*cos(theta)*sin(beta);      % y-component of Thrust

    dragx = Fdrag*cos(theta)*cos(beta);          % x-component of drag
    dragz = Fdrag*sin(theta);                    % z-component of drag
    dragy = Fdrag*cos(theta)*sin(beta);          % y-component of drag

    mSystem = mBottle + m; % Mass of the system for use in F=ma
    Fgravity = mSystem*g; % Force from Gravity

    dvx = (thrustx-dragx)/mSystem; % Acceleration in x --> integrates to vel x
    dvz = (thrustz-dragz-Fgravity)/mSystem; % Same as above, in z direction
    dvy = (thrusty-dragy)/mSystem; % acceleration in y direction

    % Now that it's falling, negatives are fine

    dx = vx; % velocity x --> integrates to distance x
    dy = vy; % same as above, in y direction
    dz = vz; % same as above, in z direction

    %% Output

    ds2=[mdot;dvx;dvz;dvy;dx;dz;dy]; % output call for ode45
end

end
