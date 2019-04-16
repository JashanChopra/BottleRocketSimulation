function [ds,Fthrust] = state1funcDI(t,state1Vec,Cd,Athroat,gamma,rhoAtm,Patm,Pgage,volAir,Abottle,mBottle,mAir,CDrag,volBottle,rhoWater,vw,fitobject)
% This function serves as the main input for ode45 for stage one, as well
% as when the rocket is on the rail
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Authors: Jason Balke (107025894), Jashan Chopra (107689146)
% Date Created: Novemeber 20th, 2018
% Last Date Modified: November 30th, 2018


% Initial conditions
g = 9.81; %[m/s^2]

vol = state1Vec(1); % Volume of the air
vx = state1Vec(2); % x-component of velocity
vz = state1Vec(3);% z-component of velocity
vy = state1Vec(4); % y-component of velocity
x = state1Vec(5);% x position
z = state1Vec(6);% z position
y = state1Vec(7); % y position

% If statement checks the distances and compares to the distance of the
% rail. Y distance should not change on the test stand.
if sqrt(x^2+z^2) < 0.531678
   dvol = Cd*Athroat*sqrt((2/rhoWater)*(((Pgage+Patm)*(volAir/vol)^gamma)-Patm)); % Volume Function
   P = ((Pgage+Patm)*volAir^gamma)/ (vol^gamma);                                  % Pressure Function

   % Thrust is intepolated from the styleSpline fit at given ode45 timestep
   Fthrust = fitobject(t);

   % Drag calculation from relative velocity
   vRel = [vx + vw(1), vy + vw(2), vz + vw(3)];         % add wind components
   velocity = norm(vRel);                               % Find the vector magnitude
   Fdrag = (1/2)*rhoAtm*(velocity^2)*CDrag*Abottle;     % Drag Force Equation

   theta = pi/4;    % Theta value constant for rails
   beta = 0;        % Beta value is constant for the rails

   thrustx = Fthrust*cos(theta)*cos(beta);      % x-component of Thrust
   thrustz = Fthrust*sin(theta);                % z-component of Thrust
   thrusty = Fthrust*cos(theta)*sin(beta);      % y-component of Thrust

   dragx = Fdrag*cos(theta)*cos(beta);          % x-component of drag
   dragz = Fdrag*sin(theta);                    % z-component of drag
   dragy = Fdrag*cos(theta)*sin(beta);          % y-component of drag

   volWater = volBottle - vol; % Volume of the water at a given point
   mWater = rhoWater * volWater; % Equation for the change in mass of water

   mSystem = mWater + mBottle + mAir; % Mass of the system for use in F=ma
   Fgravity=mSystem*g; % Force of gravity always acts in the z direction

   dvx = (thrustx-dragx)/mSystem; % acceleration in x direction
   dvz = (thrustz-dragz-Fgravity)/mSystem; % acceleration in zx direction
   dvy = (thrusty-dragy)/mSystem; % acceleration in y direction

    % If the velocities are negative make the distance 0 so it does not
    % appear to go through the test stand
    if sign(vx) == -1
        dx = 0;
    else
        dx = vx; % velocity x --> integrates to distance x
    end
    if sign(vy) == -1
        dy = 0;
    else
        dy=vy; % same as above, in y direction
    end
    if sign(vz) == -1
        dz = 0;
    else
        dz=vz; % same as above, in z direction
    end

   %% Output
   ds = [dvol;dvx;dvz;dvy;dx;dz;dy];
else % Off the rails
    dvol = Cd*Athroat*sqrt((2/rhoWater)*(((Pgage+Patm)*(volAir/vol)^gamma)-Patm)); % Volume Function

    % Thrust is intepolated from the styleSpline fit at given ode45 timestep
    Fthrust = fitobject(t);

    % Drag calculation from relative velocity
    vRel = [vx + vw(1), vy + vw(2), vz + vw(3)];        % add wind components
    velocity = norm(vRel);                              % Find the vector magnitude
    Fdrag = (1/2)*rhoAtm*(velocity^2)*CDrag*Abottle;    % Drag Force Equation

    theta = atan(vRel(3)/vRel(1));    % The angle value at any point
    beta = atan(vRel(2)/vRel(1));     % Lateral angle value at any point

    thrustx = Fthrust*cos(theta);   % x-component of Thrust
    thrustz = Fthrust*sin(theta);   % z-component of Thrust
    thrusty = Fthrust*sin(beta);    % y-component of Thrust

    dragx = Fdrag*cos(theta); % x-component of drag
    dragz = Fdrag*sin(theta); % z-component of drag
    dragy = Fdrag*sin(beta); % y-component of drag

    volWater = volBottle - vol;     % Volume of the water at a given point
    mWater = rhoWater * volWater;   % Equation for the change in mass of water

    mSystem = mWater + mBottle + mAir;  % Mass of the system for use in F=ma
    Fgravity=mSystem*g;                 % Force of gravity always acts in the z direction

    dvx = (thrustx-dragx)/mSystem;              % acceleration in x direction
    dvz = (thrustz-dragz-Fgravity)/mSystem;     % acceleration in y direction
    dvy = (thrusty-dragy)/mSystem;              % acceleration in y direction

    dx = vx; % velocity x --> integrates to distance x
    dy = vy;
    dz = vz;


    %% Output
    ds = [dvol;dvx;dvz;dvy;dx;dz;dy];
end

end
