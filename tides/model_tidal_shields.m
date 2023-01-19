clc
clear
close all
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Author: Una Schneck (schneck.una@gmail.com)
% 
% This code 
%--------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
% External function requirements:
%   (1) make_wave_COE.m -- wind exchanges energy with liquid and produces waves that grow to
%   breaking
%   (2) shoal_wave_COE.m -- wave shoals along a constant slope and computes where the waves are
%   predicted to break
%   (3) figure_settings.m -- makes figures look nice (not strictly neccessary to make calculations)
%--------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
% Code produces:
%   (1) Shields diagram of particle Reynolds number vs Shields number for unidirectional flow within
%   the nearshore and straits of Titan's lakes for previous modelled tidal flows 
%   (2) Shields diagram with error for composition  (pure ethane to pure methane) 
%---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
% %Aprox time to run: <1 minute
%---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
% Cite: The Shoreline Shaping Capability of Waves and Tides at Titan's Lakes (Schneck et al.)
% 
% --------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
tic
%% --------------------------------------------------- 1 . CHARACTER OF GRAINS AND FLOW --------------------------------------------------------- %% 

g = 1.352; 

% SEDIMENT PROPERTIES (ICE ICE-ORGANICS ORGANICS)
rho_s_pos = [940 800 1500]; % kg/m3 [Ice Ice-Organic Organic] (Burr 2006, Witek and Czechowski 2014)
d50 = [6.35e-5:1e-4:0.1]; % m [Fine-Sand: 10 cm gravel]

% PREVIOUS TIDAL MODEL PARAMETERS AND RESULTS
rho_strait = 550; % kg/m3 (Vincent 2018)
rho_lake = 662; % kg/m3 (Vincent 2016)

dyn_vis_lake = 1736e-6; % Pa.S (Lorenz 2010)
kin_vis_lake = dyn_vis_lake/rho_lake; % m2/s
kin_vis_strait = 3e-7; % m2/s (Vincent 2018)

fluid_strait = [rho_strait kin_vis_strait];
fluid_lake = [rho_lake kin_vis_lake];

max_depth_strait = 15; % m (V-shaped basin) (Vincent 2018)
depth_strait = max_depth_strait/2; % use the average? the reported velocity would be a depth-averaged one
depth_lake = 3; % ~3 m (Vincent 2016, emailed to ask for better precision but not response yet

man_coef_max_strait = 0.03; % Vincent 2018
man_coef_min_strait = 0.06; % Vincent 2018
man_coef_lake = 0.03; % Vincent 2016 

u_max_strait = 0.64; % m/s (Vincent 2018) table 3
u_min_strait = 0.12; % m/s (Vincent 2018) table 3
u_max_lake = 0.046; % m/s (Vincent 2016)  
u_min_lake = 0.02; % m/s (Vincent 2016)


%% --------------------------------------------------- 2 . CALCULATE RE AND SHIELDS --------------------------------------------------------------- %% 

[Re_strait_max,shields_strait_max] = tidal_shields(u_max_strait,depth_strait,rho_s_pos,fluid_strait,man_coef_max_strait,d50,g);
[Re_strait_min,shields_strait_min] = tidal_shields(u_min_strait,depth_strait,rho_s_pos,fluid_strait,man_coef_min_strait,d50,g);
[Re_lake_max,shields_lake_max] = tidal_shields(u_max_lake,depth_lake, rho_s_pos,fluid_lake,man_coef_lake,d50,g);
[Re_lake_min,shields_lake_min] = tidal_shields(u_min_lake,depth_lake, rho_s_pos,fluid_lake,man_coef_lake,d50,g);

%% --------------------------------------------------- 3 . SHIELDS DIAGRAM ------------------------------------------------------------------------ %% 

figure;
set(gcf, 'Units', 'Normalized', 'OuterPosition', [0 0 1 1]);
[h1,h2,h3] = make_shields_diagram();
hold on;

%% --------------------------------------------------- 4 .  SHIELDS PARAMETER FOR TITAN GRAINS TO SHIELD DIAGRAM ---------------------------------- %% 

% PLOT
h4 = plot(Re_strait_max(1,:),shields_strait_max(1,:),'-b','LineWidth',5); % MAX STRAIT FOR ICE
h5 = plot(Re_strait_max(2,:),shields_strait_max(2,:),'-r','LineWidth',5); % MAX STRAIT FOR ICE-ORGANIC
h6 = plot(Re_strait_max(3,:),shields_strait_max(3,:),'-g','LineWidth',5); % MAX STRAIT FOR ORGANIC

h7 = plot(Re_lake_max(1,:),shields_lake_max(1,:),'--b','LineWidth',5); % MAX NEARSHORE FOR ICE
h8 = plot(Re_lake_max(2,:),shields_lake_max(2,:),'--r','LineWidth',5); % MAX NEARSHORE FOR ICE-ORGANIC
h9 = plot(Re_lake_max(3,:),shields_lake_max(3,:),'--g','LineWidth',5); % MAX NEARSHORE FOR ORGANIC

% FILL IN ERROR BAR FOR STRAIT
curve1 = shields_strait_max(1,:); % MAX STRAIT ICE 
curve2 = shields_strait_min(1,:); % MIN STRAIT ICE
x2 = [Re_strait_max(1,:),fliplr(Re_strait_max(1,:))];
inbetween = [curve1,fliplr(curve2)];
f1 = fill(x2,inbetween,'b','FaceAlpha',0.2,'LineStyle','none');

curve1 = shields_strait_max(2,:); % MAX STRAIT ICE-ORGANIC
curve2 = shields_strait_min(2,:); % MIN STRAIT ICE-ORGANIC
x2 = [Re_strait_max(2,:),fliplr(Re_strait_max(2,:))];
inbetween = [curve1,fliplr(curve2)];
f2 = fill(x2,inbetween,'r','FaceAlpha',0.2,'LineStyle','none');

curve1 = shields_strait_max(3,:); % MAX STRAIT ORGANIC
curve2 = shields_strait_min(3,:); % MIN STRAIT ORGANIC
x2 = [Re_strait_max(3,:),fliplr(Re_strait_max(3,:))];
inbetween = [curve1,fliplr(curve2)];
f3 = fill(x2,inbetween,'g','FaceAlpha',0.2,'LineStyle','none');

% FILL IN ERROR BAR FOR LAKE
curve1 = shields_lake_max(1,:); % MAX LAKE ICE 
curve2 = shields_lake_min(1,:); % MIN LAKE ICE
x2 = [Re_lake_max(1,:),fliplr(Re_lake_max(1,:))];
inbetween = [curve1,fliplr(curve2)];
ff1 = fill(x2,inbetween,'b','FaceAlpha',0.2,'LineStyle','none');

curve1 = shields_lake_max(2,:); % MAX LAKE ICE-ORGANIC
curve2 = shields_lake_min(2,:); % MIN LAKE ICE-ORGANIC
x2 = [Re_lake_max(2,:),fliplr(Re_lake_max(2,:))];
inbetween = [curve1,fliplr(curve2)];
ff2 = fill(x2,inbetween,'r','FaceAlpha',0.2,'LineStyle','none');

curve1 = shields_lake_max(3,:); % MAX LAKE ORGANIC
curve2 = shields_lake_min(3,:); % MIN LAKE ORGANIC
x2 = [Re_lake_max(3,:),fliplr(Re_lake_max(3,:))];
inbetween = [curve1,fliplr(curve2)];
ff3 = fill(x2,inbetween,'g','FaceAlpha',0.2,'LineStyle','none');
% 
% %% --------------------------------------------------- 5 .  INCLUDE EARTH RIVERS -------------------------------------------------------- %% 
%  (DATA FROM: S. BIRCH)
addpath('..\resources\outside_data')
load('sam_gravel_rivers_earth.mat') 
h = Depth; % depth of flow
S = Slope; % Rise/Run (good approximation for small angles)

% shear_stress = rho_w*g*h*slope
rho_w = 997; % kg/m3 for freshwater 
rho_s = 1552; %kg/m3 (silicate sand)
spec_weight = (rho_s/rho_w) - 1; %dimensionless
g = 9.81; %m/s2
kin_visc = 1.05e-6; %m2/s;

shear_stress = rho_w.*g.*h.*S;

%shields parameter
shields_earth_rivers = shear_stress./(g*(rho_s - rho_w).*D);

%Particle Reynolds Number
Re_particle = (D/kin_visc).*sqrt(spec_weight.*g.*D);

% clean out repeated values
shields_earth_rivers(77:144) = []; 
Re_particle(77:144) = [];

% add Earth rivers to plot
Re_e = Re_particle;   
shields_earth = shields_earth_rivers;

h10 = plot(Re_e,shields_earth,'d','MarkerFaceColor',[1 0.5 0],'MarkerEdgeColor',[1 0.5 0],'MarkerSize',10);


%% --------------------------------------------------- 6 .  PLOTTING PARAMETERS -------------------------------------------------------- %% 

% PLOT FORMATTING 
xlim([10^-2 10^7])
ylim([10^-3 10^3])
set(gca,'Children',[ff1 ff2 ff3 f1 f2 f3 h1 h2 h3 h4 h5 h6 h7 h8 h9 h10]) % order of data display to put Earth rivers behind the Titan curve
hold off;

%% --------------------------------------------------- 6 .  PURE METHANE/ETHANE VARIABILITY---------------------------------------------- %% 


% at 90K (Steckloff et al. 2020)
rho_methane = 540; % kg/m3
rho_ethane = 660; % kg/m3

kin_vis_methane = 3e-7; %m2/s
kin_vis_ethane = 1.8e-6; %m2/s

fluid_methane = [rho_methane kin_vis_methane];
fluid_ethane = [rho_ethane kin_vis_ethane];

[Re_strait_methane,shields_strait_methane] = tidal_shields(u_max_strait,depth_strait,rho_s_pos,fluid_methane,man_coef_max_strait,d50,g);
[Re_strait_ethane,shields_strait_ethane] = tidal_shields(u_max_strait,depth_strait,rho_s_pos,fluid_ethane,man_coef_max_strait,d50,g);

[Re_lake_methane,shields_lake_methane] = tidal_shields(u_max_lake,depth_lake, rho_s_pos,fluid_methane,man_coef_lake,d50,g);
[Re_lake_ethane,shields_lake_ethane] = tidal_shields(u_max_lake,depth_lake, rho_s_pos,fluid_ethane,man_coef_lake,d50,g);

figure;
set(gcf, 'Units', 'Normalized', 'OuterPosition', [0 0 1 1]);
[hh1,hh2,hh3] = make_shields_diagram();
hold on;

h4 = plot(Re_strait_max(1,:),shields_strait_max(1,:),'-b','LineWidth',5); % MAX STRAIT FOR ICE
h5 = plot(Re_strait_max(2,:),shields_strait_max(2,:),'-r','LineWidth',5); % MAX STRAIT FOR ICE-ORGANIC
h6 = plot(Re_strait_max(3,:),shields_strait_max(3,:),'-g','LineWidth',5); % MAX STRAIT FOR ORGANIC

h7 = plot(Re_lake_max(1,:),shields_lake_max(1,:),'--b','LineWidth',5); % MAX NEARSHORE FOR ICE
h8 = plot(Re_lake_max(2,:),shields_lake_max(2,:),'--r','LineWidth',5); % MAX NEARSHORE FOR ICE-ORGANIC
h9 = plot(Re_lake_max(3,:),shields_lake_max(3,:),'--g','LineWidth',5); % MAX NEARSHORE FOR ORGANIC

% FILL IN ERROR BAR FOR STRAIT
curve1 = shields_strait_methane(1,:); % METHANE ICE STRAIT
curve2 = shields_strait_ethane(1,:); % ETHANE ICE STRAIT
x2 = [Re_strait_methane(1,:),fliplr(Re_strait_methane(1,:))];
inbetween = [curve1,fliplr(curve2)];
f1 = fill(x2,inbetween,'b','FaceAlpha',0.2,'LineStyle','none');

curve1 = shields_strait_methane(2,:); % METHANE ICE-ORGANIC STRAIT
curve2 = shields_strait_ethane(2,:); % ETHANE ICE-ORGANIC STRAIT
x2 = [Re_strait_methane(2,:),fliplr(Re_strait_methane(2,:))];
inbetween = [curve1,fliplr(curve2)];
f2 = fill(x2,inbetween,'r','FaceAlpha',0.2,'LineStyle','none');

curve1 = shields_strait_methane(3,:); % METHANE STRAIT ORGANIC 
curve2 = shields_strait_ethane(3,:); % ETHANE STRAIT ORGANIC
x2 = [Re_strait_methane(3,:),fliplr(Re_strait_methane(3,:))];
inbetween = [curve1,fliplr(curve2)];
f3 = fill(x2,inbetween,'g','FaceAlpha',0.2,'LineStyle','none');

curve1 = shields_lake_methane(1,:); % METHANE LAKE ICE 
curve2 = shields_lake_ethane(1,:); % EHTNAE LAKE ICE
x2 = [Re_lake_methane(1,:),fliplr(Re_lake_methane(1,:))];
inbetween = [curve1,fliplr(curve2)];
ff1 = fill(x2,inbetween,'b','FaceAlpha',0.2,'LineStyle','none');

curve1 = shields_lake_methane(2,:); % METHANE LAKE ICE-ORGANIC
curve2 = shields_lake_ethane(2,:); % ETHANE LAKE ICE-ORGANIC
x2 = [Re_lake_methane(2,:),fliplr(Re_lake_methane(2,:))];
inbetween = [curve1,fliplr(curve2)];
ff2 = fill(x2,inbetween,'r','FaceAlpha',0.2,'LineStyle','none');

curve1 = shields_lake_methane(3,:); % MAX LAKE ORGANIC
curve2 = shields_lake_ethane(3,:); % MIN LAKE ORGANIC
x2 = [Re_lake_methane(3,:),fliplr(Re_lake_methane(3,:))];
inbetween = [curve1,fliplr(curve2)];
ff3 = fill(x2,inbetween,'g','FaceAlpha',0.2,'LineStyle','none');

% PLOT FORMATTING 
xlim([10^-2 10^7])
ylim([10^-3 10^3])

%set(gca,'Children',[ff1 ff2 ff3 f1 f2 f3 h1 h2 h3 h4 h5 h6 h7 h8 h9]) % order of data display to put Earth rivers behind the Titan curve
%legend([h1 h3 h4 h5 h6 h7 h8 h9 h10],'\theta_{cr,s}','\theta_{cr}','Strait,Ice','Strait,Ice-Organics','Strait,Organics','Nearshore,Ice','Nearshore,Ice-Organics','Nearshore,Organics','Terrestrial Rivers','Location','eastoutside','fontsize',14)
%legend([h1 h3 h4 h5 h6 h7 h8 h9],'\theta_{cr,s}','\theta_{cr}','Strait,Ice','Strait,Ice-Organics','Strait,Organics','Nearshore,Ice','Nearshore,Ice-Organics','Nearshore,Organics','Location','eastoutside','fontsize',14)
addpath('..\resources\helper_functions')
figure_settings(20,'Bold',1)
toc
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% WORK CITED
%   1. Burr 2006 -- "Sediment Transport by Liquid Surficial Flow: Application to Titan"
%   2. Witek and Czechowski 2014 -- "Dynamic Modelling of River Deltas on Titan and Earth"
%   3. Vincent 2018 -- "A Numerical Study of Tides in Titan's Northern Seas, Kraken and Ligeia Mare"
%   4. Vincent 2016 -- "Numerical Study of Tides in Ontario Lacus, a Hydrocarbon Lake on the Surface
%   of the Saturnian moon Titan"
%   5. Lorenz 2010 -- "Threshold of Wave Generationon Titan's Lakes and Seas: Effect of Viscosity
%   and Implication for Cassini Observations"
%   6. Steckloff et al. 2020 -- "Stratification Dynamics of Titan's Lakes via Methane Evaporation"
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
