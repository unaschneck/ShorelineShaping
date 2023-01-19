clc
clear
close all
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Author: Una Schneck (schneck.una@gmail.com)
% 
% 
%--------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
% External function requirements:
%   (1) Need to smooth bathymetry using TADPOLE
%--------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
% Code produces:
%   (1) Ligeia Mare Entrainment Contours
%   (2) Ontario Lacus Entrainment Contours
%---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
% %Aprox time to run: <1 minute
%---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
% Cite: The Shoreline Shaping Capability of Waves and Tides at Titan's Lakes (Schneck et al.)
% 
% --------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% --------------------------------------------------- 1 . PREPPING THE RUN ------------------------------------------------ %% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

tic

addpath('..\resources\lake_bathymetries\SAR_bathy_smoothed')

load('lm_smoothed.mat')

[rows, columns, numberOfColorChannels] = size(lm_smooth);

x_extent = 1:columns; % 401
y_extent = 1:rows; % 401
[X,Y] = meshgrid(x_extent,y_extent);

figure;
contour3(X,Y,lm_smooth,[0:50:300],'-k','ShowText','on','LabelSpacing',1000)
hold on;
[Mcon,Ccon] = contour3(X, Y, lm_smooth,[6.7 6.7],'Color',[62/255 150/255 81/255],'LineWidth',2,'ShowText',0); % entrainment depth is found in entrainment_depth_bathtub.m
[Mcon1,Ccon1] = contour3(X, Y, lm_smooth,[15.8 15.8],'Color',[204/255 37/255 41/255],'LineWidth',2,'ShowText',0); % entrainment depth is found in entrainment_depth_bathtub.m
view(0,90)
%colormap winter
set(gca, 'XDir','reverse')
set(gca, 'YDir','reverse')
set(gca,'YTickLabel',[]);
set(gca,'XTickLabel',[]);
grid off;


load('ol_smoothed.mat')

[rows, columns, numberOfColorChannels] = size(ol_smooth);

x_extent = 1:columns; 
y_extent = 1:rows; 
[X,Y] = meshgrid(x_extent,y_extent);

figure;
contour3(X,Y,ol_smooth,[0:10:80],'-k','ShowText','on','LabelSpacing',1000)
hold on;
[Mcon,Ccon] = contour3(X, Y, ol_smooth,[6.7 6.7],'Color',[62/255 150/255 81/255],'LineWidth',2,'ShowText',0);
[Mcon1,Ccon1] = contour3(X, Y, ol_smooth,[15.8 15.8],'Color',[204/255 37/255 41/255],'LineWidth',2,'ShowText',0);
view(0,90)
set(gca, 'XDir','reverse')
set(gca, 'YDir','reverse')
set(gca,'YTickLabel',[]);
set(gca,'XTickLabel',[]);
grid off;
