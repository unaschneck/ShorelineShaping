% 
clc
clear
close all
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Author: Una Schneck (schneck.una@gmail.com)
% 
% This code tests the functionality of the helper functions make_wave_coe.m and shoal_wave_coe.m
% Important tests:
%       (1) The deepwater wave height, wavelength, phase speed matches the shallow water wave character in deep water?
%                (this figure is produced in the shoal wave code)
%       (2) if max_fetch < near_shore distance then the point will be skipped
%--------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
% External function requirements:
%   (1) make_wave_COE.m -- wind exchanges energy with liquid and produces waves that grow to
%   breaking
%   (2) shoal_wave_COE.m -- wave shoals along a constant slope and computes where the waves are
%   predicted to break
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


% ------------------------------------ EXPERIMENTAL VARIABLES AND CONSTANTS ------------------------------------------------ %  
max_depth_ontario = 50; % Mastrogiuseppe 2017
g = 1.352; % titan gravity
d = max_depth_ontario:-10^-3:10^-4; % waterdepth approaching shore; % water depth [m], deep water needs to be (d/L > 1/2), slope from Hayes 2010
near_shore = length(d);
max_fetch = 1e6; 

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

% MAKE WAVES
[T,H0,L0,C0] = make_wave_COE(u_pos,max_fetch,d,rho_a,rho_w,kin_visc,'Titan',0);
% SHOAL WAVES at breaking wave depth
[d0_breaking,um_breaking] = shoal_wave_COE(T,H0,L0,C0,alpha_0,rho_w,d,u_pos,'Titan',1,max_fetch,0);

disp('test_wave_model.m completed')
