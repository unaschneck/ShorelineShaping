clc
clear
close all
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Author: Una Schneck (schneck.una@gmail.com)
% 
% This code considers the threshold for grain entrainment on Titan by the oscillatory action of waves with
% the inclusion of a wave friction factor (fw) scaling term for the Shields parameter:
%                tau = 1/2*rho*g*fw*u^2 
% The exact form for the wave friction factor is not settled so in this code I am considering a
% range of fw that encompasses both hydraulically rough and hydraulically smooth:
%                (1) Rough Bed: (Soulsby & Whitehouse 1997)
%                       fw = 0.00251.*exp(5.21.*(((um*T)./(4*pi.*(d50))).^-0.19));
%                (2) Rough Bed: (Jonsson 1963)
%                       fw = 0.0604./(log(30*del./d50))     
%                           where del is the thickness of the laminar sublayer (goes to infinity as the grain size becomes comprable to the sublayer thickness, i.e. theory breakdowns)
%                (3) Smooth Bed: (Jonsson 1953)
%                       fw = 2*((dyn_visc/(0.5*rho*um*d0))^0.5);
%                (4) Smooth bed: (Schlicting 1968)
%                       fw = (2*(((w*dyn_visc)/rho_s)^0.5))/um
% Here the waves are being produced by surface winds (0.4 m/s to 3.3 m/s) via energy conservation
% besides viscous dissipation [make_wave_COE] and then allowed to shoal [shoal_wave_COE] along a
% constant slope consistent with measurements made at Ontario Lacus (Hayes et al. 2010)
%--------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
% External function requirements:
%   (1) make_wave_COE.m -- wind exchanges energy with liquid and produces waves that grow to
%   breaking
%   (2) shoal_wave_COE.m -- wave shoals along a constant slope and computes where the waves are
%   predicted to break
%   (3) figure_settings.m -- makes figures look nice (not strictly neccessary to make calculations)
%--------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
% Code produces:
%   (1) Shields diagram (dimensionless grain diameter vs Shields parameter) for surface winds of 0.4
%   m/s
%   (2) Shields diagram (dimensionless grain diameter vs Shields parameter) for surface winds of 3.3
%   m/s
%   (3) oscillatory_fw_wave_log.txt -- .txt log of commands run and some relevant final calculations
%---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
% %Aprox time to run: <5 minutes
%---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
% Cite: The Shoreline Shaping Capability of Waves and Tides at Titan's Lakes (Schneck et al.)
% 
% --------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% --------------------------------------------------- 1 . PREPPING THE RUN ------------------------------------------------ %% 
tic
% prepare log file for commands
diary log_fw_oscillatory_shields.txt;
diary on;

% ------------------------------------ EXPERIMENTAL VARIABLES AND CONSTANTS ------------------------------------------------ %  
max_depth_ontario = 50; % Mastrogiuseppe 2017
g = 1.352; % titan gravity
d = max_depth_ontario:-10^-3:10^-4; % waterdepth approaching shore; % water depth [m], deep water needs to be (d/L > 1/2), slope from Hayes 2010
near_shore = length(d);
max_fetch = 1e8; 

% ------------------------- EXPERIMENTAL CONSTANTS (from Earth experiments) ------------------------------------------------ %
cd = 0.002; %drag coef (dimensionless) (Lorenz and Hayes 2012)
a = 0.6; % pore space correction for medium sand
K = 0.77; % <<------ From earth experimental estimates

% ------------------------------------------------------ PROPERTIES OF FLOW ------------------------------------------------ %
rho_a = 1.2; %kg/m3 (atmospheric density at Titan)
rho_s = 940;% kg/m3 (ice grain density, water-ice)
rho_w = 590;% kg/m3 (liquid hydrocarbon density) Hayes 2016
d50 = [6.35e-5:1e-4:0.1];% diameters for [finegrain to 10 cm] [m]
kin_visc = 7.5e-7; %m2/s % Hayes 2016 (with Steckloff 2020)
dyn_visc = kin_visc*rho_w;
s = rho_s/rho_w;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% -----------------------------------------------------------2. MAKE WAVES ------------------------------------------------ %%

u_pos = [0.4 3.3]; % <--- surface wind speeds producing the waves [m/s] (from Titan GCM models Lora 2015)
alpha_0 = [0 0];

fprintf('wind speeds: %.2f and %.2f\n',u_pos(1),u_pos(2)); 

addpath('.\user_defined_functions\')
% MAKE WAVES
[T,H0,L0,C0] = make_wave_COE(u_pos,max_fetch,d,rho_a,rho_w,kin_visc,'Titan',0)
% SHOAL WAVES at breaking wave depth
[d0_breaking,um_breaking,~,~] = shoal_wave_COE(T,H0,L0,C0,alpha_0,rho_w,d,u_pos,'Titan',1,max_fetch,0)


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% --------------------------------- 3 . MAKE UNIDIRECTIONAL SHIELD DIAGRAM ------------------------------------------------ %%
D_star = (((g*(s-1))/(kin_visc^2))^(1/3)).*d50; % dimensionless grain diameter
D_star_wavelength = (((g*(s-1))/(kin_visc^2))^(1/3)).*(2.17/100); % dimensionless grain diameter at size of Ku-band wavelength
shield_orig = (0.24./D_star) + (0.055*(1 - exp(-0.02.*D_star))); % critical Shields number in terms of dimensionless grain diameter
shield_orig2 = (0.3./(1+1.2.*D_star)) + (0.055*(1 - exp(-0.02.*D_star))); % correction for critical Shields number (for oscillatory flow in unidirectional shields diagram) from Soulsby & Whitehouse 1997

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% -------------------------------------------4. Shields parameter with fw  ------------------------------------------------ %%
%      Soulsby & Whitehouse 1997 (Rough Bed)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% (1) SOULSBY AND WHITEHOUSE 1997 (Hydraulically-Rough Flows) (max fw = 0.3)

for i = 1:length(u_pos)
    fw_sw(i,:) = 0.00251.*exp(5.21.*(((um_breaking(i)*T(i))./(4*pi.*(d50))).^-0.19)); % wave friction factor
    shield_sw(i,:) = ((0.5*rho_w.*fw_sw(i,:).*(um_breaking(i)^2)))./(g*(rho_s - rho_w).*d50); % shields parameter with wave friction factor
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% --------------------------------------------------5. More fw equations ------------------------------------------------ %%
%       Schlichting 1968 (Smooth Bed), Jonnson 1967 (Smooth Bed), Jonsson 1963 (Rough Bed)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Schlichting 1968 and Jonnson 1967 Smooth are the same

for i = 1:length(u_pos)

    % wave friction factor:
    %fw_jon_smooth(i,:) = 2*((dyn_visc/(0.5*rho_w*um_breaking(i)*d0_breaking(i)))^0.5); % smooth fw (Jonnson 1967)
    del = 6.5*((dyn_visc*T(i)/(2*rho_w*pi))^0.5); % viscous sublayer thickness [m], using Li 1954 (as done in Le Roux 2003)
    fw_jon_rough(i,:) = 0.0604./((log10(30*del./d50)).^2); % there is a slight variation to this term in Jonsson 1980 as wel Eq. 7.10-11
    w = (2*pi)/T(i); % angular frequency
    fw_schlict(i,:) = (2*(sqrt((w*dyn_visc)/rho_w)))/um_breaking(i); % smooth fw (Schlichting 1968)

    % shields parameters with wave friction factor at the breaking wave depth
    %       tau = (1/2*rho*g*u^2)/((rho_s - rho)g*D_50
    %shield_jon_smooth(i,:) = (0.5*rho_w.*fw_jon_smooth(i,:).*(um_breaking(i)^2))./(g*(rho_s - rho_w).*d50); 
    shield_jon_rough(i,:) = (0.5*rho_w.*fw_jon_rough(i,:).*(um_breaking(i)^2))./(g*(rho_s - rho_w).*d50);
    shield_schlict(i,:) = (0.5*rho_w.*fw_schlict(i,:).*(um_breaking(i)^2))./(g*(rho_s - rho_w).*d50);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% -----------------------------------------------------6. Making the plots ------------------------------------------------ %% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% PLOT: Shields diagram (dimensionless grain diameter vs Shields parameter) 
figure;
h1 = plot(D_star,shield_orig,'-k','LineWidth',5);
hold on
%h2 = plot(D_star,shield_orig2,'--k','LineWidth',5);
h2 = plot(D_star,shield_sw(2,:),'--k','LineWidth',2);
set(gca, 'YScale', 'log')
set(gca, 'XScale', 'log')
xlabel('D_{**}')
ylabel('\theta')
X=[D_star,fliplr(D_star)];                
Y=[shield_schlict(2,:),fliplr(shield_sw(2,:))];  
ff1 = fill(X,Y,[0.1215969 0.470588 0.705882],'FaceAlpha',0.5);
h3 = plot(D_star,shield_schlict(2,:),':k','LineWidth',2);
h4 = plot(D_star,shield_schlict(1,:),':k','LineWidth',2);
grid on;
box on;
ylim([10^-3 10^3])
h5 = plot(D_star,shield_sw(1,:),'--k','LineWidth',2);
h6 = plot(D_star_wavelength.*ones(size([10^-3 1 10^3])),[10^-3 1 10^3],'--k');
X=[D_star,fliplr(D_star)];                
Y=[shield_schlict(1,:),fliplr(shield_sw(1,:)),];  
ff2 = fill(X,Y,[0.858777 .0985842 0.46287],'FaceAlpha',0.5);
lg = legend([h1 ff1 ff2],'$\theta_{cr}$','High Wind','Low Wind');
set(lg,'Interpreter','latex')


diary off; % close logging function
addpath('..\resources\helper_functions')
figure_settings(20,'Bold',1)
toc
disp('fw_oscillatory_shields.m completed')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% References
% (1) Madsen & Grant 1975 --- "The Threshold of Sediment Movement Under Oscillatory Waves; A Discussion"
% (2) Soulsby & Whitehouse 1997 --- "Threshold of Sediment Motion in Coastal Enviroments"
% (3) Schlichting & Kestin 1968 --- "Boundary Layer Theory"
% (4) Jonsson 1963 --- "Measurements in Turbulent Wave Boundary Layers"
% (5) Jonsson 1967 --- "Wave Boundary Layers and Friction Factors"
% (6) Li 1954 --- "Stability of Oscillatory Laminar Flow along a Wall"
% (7) Le Roux 2003 --- "Wave Friction Factor as Related to the Shields Parameter for Steady Currents"
% (8) Hayes et al. 2010 --- "Bathymetry and Absorbtivity of Titan's Ontario Lacus"
% (9) Mastrogiuseppe et al. 2017 --- "Bathymetry and Composition of Titan's Ontario Lacus Derived from
% Monte Carlo-Based Waveform inversion of Cassini RADAR Altimetry Data"
% (10) Lorenz and Hayes 2012 --- "The Growth of Wind-Waves in Titan's Hydrocarbon Seas"
% (11) Hayes 2016 --- "The Lakes and Seas of Titan"
% (12) Steckloff et al. 2020 --- "Stratification Dynamics of Titan's Lakes Via Methane Evaporation"
% (13) Lora et al. 2015 --- "GCM Simulations of Titan's Middle and Lower Atmosphere and Comparisons to
% Observations"
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


