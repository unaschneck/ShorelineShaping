function [Re,shields] = tidal_shields(flow_speed,flow_depth, grain_density,fluid_type,mannings_coef,d50,g)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Objective: 
%   This function will return the particle Reynolds number and the Shields number for unidirectional
%   flow (i.e. tides).
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% inputs:
%   flow_speed = velocity of liquid [m/s]
%   flow_depth = depth of liquid [m]
%   grain_density = density of non-cohesive grains [kg/m3]
%   fluid_type = information about fluid
%       fluid_type(1) = fluid density [kg/m3]
%       fluid_type(2) = kin_viscocity [m2/s]
%   manning_coef = Manning coefficient of bed 
%   d50 = diameter of grains [m]
%   g = gravity [m2/s]
% outputs:
%   Re = Particle Reynolds number
%   shields = Shields number
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% prepare variables in form to be read into equation
u = flow_speed;
H = flow_depth;
rho_s = grain_density;
rho = fluid_type(1);
kin_visc = fluid_type(2);
man_coef = mannings_coef;

% find Re and Shields number for each grain type and grain size
for i = 1:length(rho_s)
% PARTICLE REYNOLDS NUMBER
    Re(i,:) = (d50.*sqrt((rho_s(i)/rho - 1).*g.*d50))./kin_visc; % particle reynolds number
    
    % SHIELDS NUMBER
    shields(i,:) = (rho*man_coef*u^2)./((H^(1/3)).*(rho_s(i)-rho).*g.*d50);
end

end