clc
clear
close all
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Author: Una Schneck (schneck.una@gmail.com)
% 
% This code computes the unidirectionally flow speed at different depth within the strait and the nearshore.
%
% u(z) = (u_star/k)*ln(z/z0) <- Law of the Wall (with 4/10 rule which says that the depth-averaged flow is aprox equal to the flow at 4/10 of the depth)
%   z = from z0 to H
%   k = 0.4 (Von Karman Constant, dimensionless)
%       z0 = D/30 <-- assume rough turbulence (constant z0)
%       z0 = kin_visc/(9*u_star) <-- smooth turbulence 
%--------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
% External function requirements:
%   (1) figure_settings.m -- makes figures pretty
%--------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
% Code produces:
%   (1)  Plot of flow speed with depth for smooth and rough hydraulically flow in strait
%   (2)  Plot of flow speed with depth for smooth and rough hydraulically flow in nearshore
%---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
% %Aprox time to run: < 1 minute
%---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
% Cite: The Shoreline Shaping Capability of Waves and Tides at Titan's Lakes (Schneck et al.)
% 
% --------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
tic
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% --------------------------------------------------- 1 . CHARACTER OF GRAINS AND FLOW --------------------------------------------------------- %% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% CONSTANTS
k = 0.4; % Von Karman Constant
d50 = [6.35e-5:1e-4:0.01]; % diameters for [finegrain sand gravel] [m]
D_mean = mean(d50); % [m]

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% --------------------------------------------------- 2 . PREVIOUS TIDAL PARAMETERS --------------------------------------------------------- %% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% PREVIOUS TIDAL MODEL PARAMETERS
u_max_strait = 0.64; % m/s (Vincent 2018) table 3
u_min_strait = 0.12; % m/s (Vincent 2018) table 3
u_max_lake = 0.046; % m/s (Vincent 2016)  
u_min_lake = 0.02; % m/s (Vincent 2016)

u_strait = u_max_strait;
u_coast = u_max_lake;

% PREVIOUS TIDAL MODEL PARAMETERS AND RESULTS
rho_strait = 550; % kg/m3 (Vincent 2018)
rho_lake = 662; % kg/m3 (Vincent 2016)

dyn_vis_lake = 1736e-6; % Pa.S (Lorenz 2010)
kin_vis_lake = dyn_vis_lake/rho_lake; % m2/s
kin_vis_strait = 3e-7; % m2/s (Vincent 2018)

fluid_strait = [rho_strait kin_vis_strait];
fluid_lake = [rho_lake kin_vis_lake];

max_depth_strait = 15; % m (V-shaped basin) (Vincent 2018)
H_strait = max_depth_strait/2; % use the average? the reported velocity would be a depth-averaged one
H_lake = 3; % ~3 m (Vincent 2016, emailed to ask for better precision but not response yet)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% --------------------------------------------------- 3 . LAW OF THE WALL --------------------------------------------------------- %% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% A. - NEARSHORE %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% (1) HYDRAULICALLY ROUGH 
z0_coast = D_mean/30;
u_star_coast = (u_coast*k)/log(0.37*H_lake/z0_coast); % using 4/10 rule to define u_star
z_coast = z0_coast:H_lake; % depth of flow to top of log layer
u_z_coast = (u_star_coast/k).*log(z_coast./z0_coast);
% (2) HYDRAULICALLY SMOOTH 
z0_coast_smooth = kin_vis_lake/(9*u_star_coast);
u_star_coast_smooth = (u_coast*k)/log(0.37*H_lake/z0_coast_smooth); % using 4/10 rule to define u_star
z_coast_smooth = z0_coast_smooth:H_lake;
u_z_coast_smooth = (u_star_coast_smooth/k).*log(z_coast_smooth./z0_coast_smooth);

% B. - STRAIT %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% (1) HYDRAULICALLY ROUGH 
z0_strait = D_mean/30;
u_star_strait = (u_strait*k)/log(0.37*H_strait/z0_strait); % using 4/10 rule to define u_star
z_strait = z0_strait:H_strait;
u_z_strait = (u_star_strait/k).*log(z_strait./z0_strait);
% (2) HYDRAULICALLY SMOOTH 
z0_strait_smooth = kin_vis_strait/(9*u_star_strait);
u_star_strait_smooth = (u_strait*k)/log(0.37*H_strait/z0_strait_smooth); % using 4/10 rule to define u_star
z_strait_smooth = z_strait:H_strait;
u_z_strait_smooth = (u_star_strait_smooth/k).*log(z_strait_smooth./z0_strait_smooth);


% percent differences between rough and smooth
difference_strait = (abs(u_z_strait - u_z_strait_smooth)./((u_z_strait + u_z_strait_smooth)./2)).*100
difference_lake = (abs(u_z_coast - u_z_coast_smooth)./((u_z_coast + u_z_coast_smooth)./2)).*100

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% --------------------------------------------------- 4 . PLOTTING --------------------------------------------------------- %% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% PLOTTING VARIABLES
vertical_line_ustrait = u_strait.*ones(1,5);
vertical_line_ucoast = u_coast.*ones(1,5);
z_fourtenths_strait = .37.*H_strait.*ones(1,5); 
z_fourtenths_lake = .37.*H_lake.*ones(1,5); 

figure;
set(gcf, 'Units', 'Normalized', 'OuterPosition', [0 0 1 1]);
h1 = plot(u_z_coast,z_coast,'--r','LineWidth',4);
hold on;
h2 = plot(u_z_coast_smooth,z_coast_smooth,'-r','LineWidth',4);
h3 = plot([0 0.4 0.5 0.6 0.8],z_fourtenths_lake,'--k','LineWidth',4);
h4 = plot(vertical_line_ucoast,[10^-6 10^-4 10^-2 10^0 10^1],'--k','LineWidth',4);
hold off;
legend([h1 h2],{'Rough','Smooth'})
xlabel('u [m/s]')
ylabel('z [m]')
grid on;
title('Vertical Velocity Profile (Nearshore)')
box on;
xlim([0 0.06])
ylim([0 H_lake])

figure;
set(gcf, 'Units', 'Normalized', 'OuterPosition', [0 0 1 1]);
hh1 = plot(u_z_strait,z_strait,'--b','LineWidth',4);
hold on;
hh2 = plot(u_z_strait_smooth,z_strait_smooth,'-b','LineWidth',4);
hh3 = plot([0 0.4 0.5 0.6 0.8],z_fourtenths_strait,'--k','LineWidth',4);
hh4 = plot(vertical_line_ustrait,[10^-6 10^-4 10^-2 10^0 10^1],'--k','LineWidth',4);
hold off;
legend([hh1 hh2],{'Rough','Smooth'})
xlabel('u [m/s]')
ylabel('z [m]')
grid on;
title('Vertical Velocity Profile (Strait)')
box on;
ylim([0 H_strait])

addpath('..\resources\helper_functions')
figure_settings(20,'Bold',1)
toc
