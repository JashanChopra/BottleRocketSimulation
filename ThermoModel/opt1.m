function [value,isterminal,direction] = opt1(t,ds)
% This function modifies to stop ode45 when the phase 1 condition is met
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Authors: Jason Balke (107025894), Jashan Chopra (107689146)
% Date Created: Novemeber 20th, 2018
% Last Date Modified: November 30th, 2018

% Stopping condition is the pressure inside the bottle is equal to the
% atomspheric condition
value = (ds(1) > .0020000000) - 0.5; % the -0.5 avoids a negative error thrown

isterminal = 1; % This will stop the ode45 simulation
direction = 0; % Misc, required
end
