function [processed_lake] = clean_lake_bathy(lake_data)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Objective: 
%   Function will produce a single lake bathymetry from raw interpreted bathymetry data with a
%   shoreline defined from edge-detection for the largest closed basin. All points outside this main
%   lake basin will be zeroed.
%   processed_lake data can be used in the tadpole diffusion model for a smoother lake bottom
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% inputs:
%   lake_data = array of depth measurements 
% output:
%   processed_lake = array of depth measurements for main basin to be used in tadpole diffusion
%   model for smoothing
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% --------------------------------------------------- 1 . LOAD IN INTERPRETED BATHY DATA ------------------------------------------------ %% 


%lake_data = smoothdata(ligeia_bathy,1,'gaussian','omitnan'); % smooth data

figure; % show original data
image(lake_data);
set(gca, 'XDir','reverse')
set(gca, 'YDir','reverse')
title('original data')
colorbar

%% --------------------------------------------------- 2 . FIND SHORELINE FROM EDGE-DETECTION---------------------------------------------- %% 

BW = im2bw(lake_data); % mask to ID lake boundary
B = bwboundaries(BW); % all boundaries of all closed polygons

figure; % show all the closed basins with the data
for k=1:length(B)
   b = B{k};
   plot(b(:,2),b(:,1),'LineWidth',3);
   axis([0 450 0 450])
   set(gca, 'XDir','reverse')
   set(gca, 'YDir','reverse')
   hold on
   title('all closed polygons')
%    drawnow
%    pause(0.1)
end

[max_size, max_index] = max(cellfun('size', B, 1)); % look for largest closed basin to assign shoreline values
X_perimeter = B{max_index}(:,2); % X-coordinates of shoreline perimeter
Y_perimeter = B{max_index}(:,1); % Y-coordinates of shoreline perimeter

 
% figure; % show shoreline (aka perimeter of largest basin)
% plot(x,y)
% title('coordinates of shoreline')

%% --------------------------------------------------- 3 . MESH FOR SURFACE PLOT ---------------------------------------------------------- %% 

% make grid of points for mesh
[rows, columns, numberOfColorChannels] = size(lake_data);
x_extent = 1:columns; 
y_extent = 1:rows; 
[X,Y] = meshgrid(x_extent,y_extent);
% need points within the lake
[in,on] = inpolygon(X,Y,X_perimeter,Y_perimeter);

% figure; % show grid of points within and outside the lake perimeter
% plot(X(~in),Y(~in),'ro')
% hold on
% plot(X(in|on),Y(in|on),'go')
% title('points within the lake')
% set(gca, 'XDir','reverse')
% set(gca, 'YDir','reverse')

% set all values outside the main shoreline to zero
lake_data(~in) = 0;   

processed_lake = lake_data;

figure; % show original data
image(processed_lake);
set(gca, 'XDir','reverse')
set(gca, 'YDir','reverse')
title('processed data for tadpole diffusion')
colorbar

end




