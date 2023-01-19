clc
clear
close all
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Author: Una Schneck (schneck.una@gmail.com)
%
% Code finds the the modified-oscillatory Shields parameter for surface waves
% produced by a constant wind across the entire fetch of generic lake (Ontario Lacus)
% at Titan and Earth. Produces the modified Shields diagram (Komar & Miller 1973)
% to assess the flow competency at the wave breaking depth.
%
% Wave generation is a simple energy model with viscous dissipation (Lorenz & Hayes 2012)
%
%--------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
% External function requirements:
%   (1) make_wave_COE.m -- wind exchanges energy with liquid and produces waves that grow to
%   breaking
%   (2) shoal_wave_COE.m -- wave shoals along a constant slope and computes where the waves are
%   predicted to break
%   (3) figure_settings.m -- makes figures look nice (not strictly neccessary to make calculations)
%--------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
% Code produces:
%   (1) Shields diagram for Earth and Titan grains for variations of thresholds originally defined by Komar & Miller 1973
%---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
% %Aprox time to run: <5 minutes
%---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
% Cite: The Shoreline Shaping Capability of Waves and Tides at Titan's Lakes (Schneck et al.)
% 
% --------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% --------------------------------------------- WAVES IN TITAN CONDITIONS ------------------------------------------------ %% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 1. CONSTANTS AND VARIABLES 
%% --------------------------------------------------- 1 . PREPPING THE RUN ------------------------------------------------ %% 
tic
% prepare log file for commands
diary log_km_oscillatory_shields.txt;
diary on;

% --------------- EXPERIMENTAL VARIABLES AND CONSTANTS ----------------------- %%
max_depth_ontario = 50; % Mastrogiuseppe 2018
g = 1.352; % titan gravity
d = max_depth_ontario:-10^-3:10^-8; % waterdepth approaching shore;  water depth [m], deep water needs to be (d/L > 1/2), slope from Hayes 2010
near_shore = length(d);
max_fetch = 1e8; 
% ---------------- EXPERIMENTAL CONSTANTS (from Earth experiments) ------------ %%
cd = 0.002; %drag coef (dimensionless) (Lorenz and Hayes 2012)
a = 0.6; % pore space correction for medium sand
K = 0.77; % <<------ From earth experimental estimates

% ---------------- PROPERTIES OF FLOW ----------------------------------------- %
rho_a = 1.2; %kg/m3 (atmospheric density)
rho_s = 940;% kg/m3 (ice grain density)
rho_w = 590;% kg/m3 (liquid hydrocarbon density) Hayes 2016
d50 = [6.35e-5:1e-4:0.1];% diameters for [finegrain to 10 cm] [m]
kin_visc = 7.5e-7; %m2/s % Hayes 2016
dyn_visc = kin_visc*rho_w;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% ------------------------------- 2. CREATE WAVE USING ENERGY CONSERVATION ------------------------------------------------ %%  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
u_pos = [0.4 3.3]; % <--- surface wind speeds producing the waves [m/s] (range from Titan GCM models Lora 2015)
alpha_0 = [0 0];
fprintf('Titan: wind speeds: %.2f and %.2f\n',u_pos(1),u_pos(2)); 

addpath('.\user_defined_functions\')
% MAKE WAVES
[T,H0,L0,C0] = make_wave_COE(u_pos,max_fetch,d,rho_a,rho_w,kin_visc,'Titan',0)
% SHOAL WAVES at breaking wave depth
[d0_breaking,um_breaking,H_breaking_T,L_breaking_T] = shoal_wave_COE(T,H0,L0,C0,alpha_0,rho_w,d,u_pos,'Titan',1,max_fetch,0)
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% --------------------------------------- 3. OSCILLATORY SHIELDS PARAMETER ------------------------------------------------ %%  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% ------------------------------------------- KOMAR & MILLER THRESHOLD FORM ------------------------------------------------ %%   
% SHIELDS PARAMETER FOR WAVE (only entraining if this number is above the critical oscillatory
% shields number) [Komar & Miller 1973]
titan_shield_komar_winter = rho_w*(um_breaking(1).^2)./((rho_s - rho_w)*g.*d50); 
titan_shield_komar_summer = rho_w*(um_breaking(2).^2)./((rho_s - rho_w)*g.*d50); 

% The lines are not actually overlapping, just very close in non-dimensional log-space. Can see
% here if uncomment these two figures:
% figure;
% plot(sqrt(1./d50),titan_shield_komar_winter,'b')
% hold on
% plot(sqrt(1./d50),titan_shield_komar_summer,'r')
% set(gca, 'YScale', 'log')
% set(gca, 'XScale', 'log')
% xlabel('sqrt(1/d50)')
% ylabel('\theta')
% legend('Titan winter','Titan summer')
% figure;
% ax1 = plot(sqrt(d0_breaking(1)./d50),titan_shield_komar_winter,'-b');
% hold on
% ax2  = plot(sqrt(d0_breaking(2)./d50),titan_shield_komar_summer,'-r');
% set(gca, 'YScale', 'log')
% set(gca, 'XScale', 'log')
% xlabel('sqrt(d0/d50)')
% ylabel('\theta')

% put in cgs
rho_s = rho_s/1000; % g/cm3
rho_w = rho_w/1000; % g/cm3
g = g*100; %cm/s2
d50 = d50.*100; % cm
dyn_visc = dyn_visc*10*100; %poise
d0_breaking = d0_breaking.*100; % cm


% Critical Oscillatory Shields Parameter for Dingler 1979
shield_ding_winter = 0.04.*((d0_breaking(1)./d50).^(2/3)).*(((rho_w*(rho_s - rho_w)*g.*(d50.^3))./(dyn_visc^2)).^(-1/9));
shield_ding_summer = 0.04.*((d0_breaking(2)./d50).^(2/3)).*(((rho_w*(rho_s - rho_w)*g.*(d50.^3))./(dyn_visc^2)).^(-1/9));

% Variation on Komar & Miller's Critical Oscillatory Shields Parameter
a_prime = [0.3 0.1 0.21 0.21 0.21 0.392 0.39 0.14 0.24]; % Komar & Miller 1973, Komar & Miller 1975, Bagnold 1945, Manohar 1955, Sternberg & Larsen 1975, Larsen et al. 1981
x1 = [10^-1:10^-1:10^3];
curve_min = min(a_prime).*x1; 
curve_max = max(a_prime).*x1;
x2 = [x1,fliplr(x1)];
inBetween = [curve_min, fliplr(curve_max)];


figure('units','normalized','outerposition',[0 0 1 1])
h13 = fill(x2,inBetween,[.7 .7 .7]); % fill area around Komar & Miller 1973's estimation for critical oscillatory shields parameter
hold on;
h1 = plot(x1,a_prime(1).*x1,'-k','linewidth',4); % Komar and Miller 1973 critical shields parameter
h9 = plot(sqrt(d0_breaking(1)./d50),shield_ding_winter,':','Color','k','linewidth',4);
h10 = plot(sqrt(d0_breaking(2)./d50),shield_ding_summer,'-.','Color','k','linewidth',4);
h11 = plot(sqrt(d0_breaking(1)./d50),titan_shield_komar_winter,'--','Color',[0.9290 0.6940 0.1250],'linewidth',4);
h12 = plot(sqrt(d0_breaking(2)./d50),titan_shield_komar_summer,'-','Color',[0.9290 0.6940 0.1250],'linewidth',4);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% REPEAT PROCESS BUT FOR EARTH NOW
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clearvars -except h13 h1 h11 h12 h9 h10 h1 h13 x1 h15 h16 L_breaking_T H_breaking_T
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 1.b. ----------------------------------- CONSTANTS AND VARIABLES (EARTH) ------------------------------------------------ %% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% ----------------------------------- EXPERIMENTAL VARIABLES AND CONSTANTS ----------------------- %%
max_depth_ontario = 50; % Mastrogiuseppe 2016
d = max_depth_ontario:-10^-3:10^-8; % waterdepth approaching shore;  water depth [m], deep water needs to be (d/L > 1/2), slope from Hayes 2010
near_shore = length(d);
max_fetch = 1e8; 

% ------------------------- EXPERIMENTAL CONSTANTS (from Earth experiments) ------------ %%
cd = 0.002; %drag coef (dimensionless) (lorenz and Hayes 2012)
a = 0.6; % pore space correction for medium sand
K = 0.77; % <<------ From earth experimental estimates
g = 9.81;


% ------------------------------------------------------ PROPERTIES OF FLOW ----------------------------------------- %
rho_a = 1; %kg/m3 (atmospheric density)
rho_s = 2650;% kg/m3 (quartz grain density)
rho_w = 1000;% kg/m3 (liquid water density) 
d50 = [6.35e-5:1e-4:0.1];% diameters for [finegrain to 10 cm] [m]
kin_visc = 1.1e-5; %m2/s 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 2.b ----------------------CREATE WAVE USING ENERGY CONSERVATION  (EARTH) ------------------------------------------------ %% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

u_pos = [0.4 3.3]; % <--- surface wind speeds producing the waves [m/s] (from Titan GCM models Lora 2015)
alpha_0 = [0 0];
fprintf('Earth: wind speeds: %.2f and %.2f\n',u_pos(1),u_pos(2)); 

% MAKE WAVES
[T,H0,L0,C0] = make_wave_COE(u_pos,max_fetch,d,rho_a,rho_w,kin_visc,'Earth',0)
% SHOAL WAVES
[d0_breaking,um_breaking,H_breaking_E,L_breaking_E] = shoal_wave_COE(T,H0,L0,C0,alpha_0,rho_w,d,u_pos,'Earth',1,max_fetch,0)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% ------------------------------ 3.b OSCILLATORY SHIELDS PARAMETER (EARTH) ------------------------------------------------ %% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    

shields_earth_winter = rho_w*(um_breaking(1)^2)./((rho_s - rho_w)*g.*d50);
shields_earth_summer = rho_w*(um_breaking(2)^2)./((rho_s - rho_w)*g.*d50);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% --------------------------------------------------------------- 4. PLOTS ------------------------------------------------ %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

h15 = plot(sqrt(d0_breaking(1)./d50),shields_earth_winter,'--','Color',[0.36863 0.23529 0.60000],'LineWidth',4);
h16 = plot(sqrt(d0_breaking(2)./d50),shields_earth_summer,'-','Color',[0.36863 0.23529 0.60000],'LineWidth',4);
set(gca, 'YScale', 'log')
set(gca, 'XScale', 'log')
grid on;
xlabel('$\sqrt{\frac{d_0}{D}}$','interpreter','latex', 'FontWeight','bold','FontSize',40)
ylabel('$\theta$','interpreter','latex', 'FontWeight','bold','FontSize',40)
xlim([1 100])
set(gca, 'XDir','reverse')


% figure;
% plot(sqrt(1./d50),shields_earth_winter,'b')
% hold on
% plot(sqrt(1./d50),shields_earth_summer,'r')
% set(gca, 'YScale', 'log')
% set(gca, 'XScale', 'log')
% xlabel('sqrt(1/d50)')
% ylabel('\theta')
% legend('Earth winter','Earth summer')
% figure;
% ax1 = plot(sqrt(d0_breaking(1)./d50),shields_earth_winter,'-b');
% hold on
% ax2  = plot(sqrt(d0_breaking(2)./d50),shields_earth_summer,'-r');
% set(gca, 'YScale', 'log')
% set(gca, 'XScale', 'log')
% xlabel('sqrt(d0/d50)')
% ylabel('\theta')

diary off; % close logging function
addpath('..\resources\helper_functions')
figure_settings(20,'Bold',1)
%legend([h11 h15 h9 h10 h1 h13],'Titan Grains','Earth Grains','D1979 |u| = 0.4 m/s','D1979 |u| = 3.3 m/s','K&M1973','a'' variations','Location','northeast','fontsize',25)
lgd1 = legend([h11 h12 h15 h16],'Titan Low-Wind','Titan High-Wind','Earth Low-Wind','Earth High-Wind','Location','NorthEast');
lgd1.FontSize = 14;
a=axes('position',get(gca,'position'),'visible','off');
lg2 = legend(a,[h1 h13 h9 h10],'$\theta_{cr}$','$\pm \theta_{cr}$','$\theta_{cr}$ Low-Wind','$\theta_{cr}$ High-Wind','Location','SouthWest');
lg2.FontSize = 14;
set(lg2,'Interpreter','latex')

disp('km_oscillatory_shields.m completed')
toc
 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% WORK CITED
% 1.  Mastrogiuseppe 2018 -- "Bathymetry and Composition of Titan's Ontario Lacus Derived from Monte Carlo-Based Waveform Inversion of Cassini RADAR Altimetry Data"
% 2.  Hayes 2010 -- "Bathymetry and Absorbtivity of Titan's Ontario Lacus"
% 3.  Lorenz and Hayes 2012 -- "The Growth of Wind-Waves in Titan's Hydrocarbon Lakes"
% 4.  Hayes 2016 -- "The Lakes and Seas of Titan"
% 5.  Davidson-Arnott -- "Introduction to Coastal Processes and Geomorphology"
% 6.  CEM -- "Coastal Engineering Manuals  (https://www.publications.usace.army.mil/USACE-Publications/Engineer-Manuals/u43544q/636F617374616C20656E67696E656572696E67206D616E75616C/)
% 7.  Komar & Miller 1973 -- "The Threshold of Sediment Movement under Oscillatory Water Waves"
% 8.  Dingler 1979 -- "The Threshold of Grain Motion under Oscillatory Flow in a Labratory Wave Channel"
% 9. Komar & Miller 1975 -- "On the Comparison between the Threshold of Sediment Motion Under Waves and Unidirectional Currents with a Discussion of the Practical Evaluation of the Threshold: REPLY"
% 10. Bagnold 1945 -- "Motion of Waves in Shallow Water. Interaction between Waves and Sand Bottoms"
% 11. Manohar 1955 -- "Mechanics of Bottom Sediment Movement due to Wave Action"
% 12. Sternberg & Larsen 1975 -- "Threshold of Sediment Movement by Open Ocean Waves: Observations"
% 13. Larsen et al. 1981 -- "Field Investigations of the Threshold of Grain Motion by Ocean Waves and Currents"
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%