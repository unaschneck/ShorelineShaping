function [h1, h2,h3] = make_shields_diagram()
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Objective:
%   This function creates a figure object with the entrainment threshold and suspension 
%   threshold with properties extended to Titan and Earth enviroments.
%   
%   Y-Axis: Shields Number
%   X-Axis: Particle Reynolds Number
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% inputs:
%   none
% outputs:
%   h1, h2, h3 = figure handles for thresholds
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% MAKE SHIELDS DIAGRAM AT EARTH:
rho_e = 1000; %water, kg/m3
rho_s_e = 2000; %quartz sand, kg/m3
spec_weight = (rho_s_e/rho_e) - 1; %dimensionless
g_e = 9.81; %m/s2
k = 0.4;
kin_visc_e = 1.05e-6; %m2/s;
D = [6.35e-6:1e-6:0.6]; %larger range of D for wider range of Re_p

%non-dimensional diameter
Re_particle = (D/kin_visc_e).*sqrt(spec_weight.*g_e.*D);
Re_p_1 = Re_particle;
D_star = Re_particle.^2;
% non-dimensional settling velocity
log_nondim_w_s = -3.76715 + 1.92944.*(log10(D_star)) - 0.09815.*(log10(D_star)).^2 - 0.00575.*(log10(D_star)).^3 + 0.00056.*(log10(D_star)).^4;
nondim_w_s = 10.^(log_nondim_w_s);
% dimensional settling velocity
w_s = (spec_weight.*g_e.*kin_visc_e.*nondim_w_s).^(1/3);

% finding the suspension threshold
numerator = (w_s./(1.2*k)).^2;
denominator = g_e.*spec_weight.*D;
suspension_threshold = numerator./denominator;
% finding entrainment threshold
Re_particle = (D/kin_visc_e).*sqrt(spec_weight.*g_e.*D);
Re_p_2 = Re_particle;
entrainment_threshold = 0.5.*(0.22.*(Re_particle.^(-0.6)) + 0.06.*(10.^(-7.7.*(Re_particle.^(-0.6)))));

% PLOTTING
h1 = plot(Re_particle(suspension_threshold>entrainment_threshold),suspension_threshold(suspension_threshold>entrainment_threshold),'--k','LineWidth',3);
hold on;
h2 = plot(Re_particle,entrainment_threshold,'-k','LineWidth',3);
set(gca, 'YScale', 'log')
set(gca, 'XScale', 'log')
grid on
xlabel('Re_p')
ylabel('\theta')
axis([10^-1 10^7 10^-3 10^3])
grid on;

% EXTENDING ENTRAINMENT THRESHOLD CURVE FOR TITAN
g = 1.352;
rho_s = 940; % ice grains
rho = 550; % liquid hydrocarbon (Vincent 2018)
kin_visc = 3e-7; % m2/s (Vincent 2018)

spec_weight = (rho_s/rho) - 1; %dimensionless

%non-dimensional diameter
Re_particle = (D/kin_visc).*sqrt(spec_weight.*g.*D);
D_star = Re_particle.^2;

% finding entrainment threshold
Re_particle = (D/kin_visc).*sqrt(spec_weight.*g.*D);
entrainment_threshold = 0.5.*(0.22.*(Re_particle.^(-0.6)) + 0.06.*(10.^(-7.7.*(Re_particle.^(-0.6)))));

h3 = plot(Re_particle,entrainment_threshold,'-k','LineWidth',3);

end