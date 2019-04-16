function [value,isterminal,direction] = opt2(t,ds2)
% This function modifies to stop ode45 when the phase 3 condition is met
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Authors: Jason Balke (107025894), Jashan Chopra (107689146)
% Date Created: Novemeber 20th, 2018
% Last Date Modified: November 30th, 2018
% Stopping condition is when it hits the ground

% If z = 0
value = (ds2(5) < 0) - 0.5; % the -0.5 avoids a negative error thrown

isterminal = 1; % This stops the ode45 simulation
direction = 0; % Misc, required
end
