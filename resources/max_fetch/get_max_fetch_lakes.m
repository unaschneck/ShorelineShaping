clc
clear
close all

for lakes = 1:2 % [OL LM]

    shorelines{1} = [];
    if lakes == 1
        load('OL_COAST_SCHNECK.mat')
        x = X_cor;
        x = x(~isnan(x));
        x = vertcat(x,x(1));
        y = Y_cor;
        y = y(~isnan(y));
        y = vertcat(y,y(1));
        %ROTATE MAP SO NORTH POLE IS IN +Y UNIT DIRECTION 
        polyin = polyshape({x},{y});
        [xc,yc] = centroid(polyin);
        theta_rot = -atan2(yc,xc) + pi/2; % in Southern Hemisphere, (0,0) = South
    elseif lakes == 2
        load('LigeiaMare_data.mat')
        x = x_lm;
        y = y_lm;

        % ROTATE MAP SO NORTH POLE IS IN +Y UNIT DIRECTION 
        polyin = polyshape({x},{y});
        [xc,yc] = centroid(polyin);
        theta_rot = -atan2(yc,xc) - pi/2; % in Northern Hemisphere (0,0) = North Pole

    else
        disp('no lake data')
    end
    
    v = [x y]';
    
    center = repmat([xc; yc], 1, length(x));
    R = [cos(theta_rot) -sin(theta_rot); sin(theta_rot) cos(theta_rot)];
    s = v - center;     % shift points in the plane so that the center of rotation is at the origin
    so = R*s;           % apply the rotation about the origin
    vo = so + center;   % shift again so the origin goes back to the desired center of rotation
    x_rotated = vo(1,:);
    y_rotated = vo(2,:);
    
    
    dp0 = [0 0] - [xc yc];
    
    s1 = [0 0]' - center;     % shift points in the plane so that the center of rotation is at the origin
    so1 = R*s1;           % apply the rotation about the origin
    vo1 = so1 + center;   % shift again so the origin goes back to the desired center of rotation
    x0_rotated = vo1(1,end);
    y0_rotated = vo1(2,end);
    
    dp = [x0_rotated y0_rotated] - [xc yc];
    
    x = x_rotated';
    y = y_rotated';
    x = x+1e6;
    y = y+1e6;
                
    shorelines{1}(:,1) = x;
    shorelines{1}(:,2) = y;

    %include visilibity program to path
    addpath 'C:\Users\schne\OneDrive\Documents\MATLAB\MIT\Tidal Flows on Titan\Visilibity_ROSE\Fetch_VisiLibity-main'
    
    %parameters for visibility functions
    eps = 1e-4;
    epsilon = 0.1;
    snap_distance = 0.05; % must make sure snap_distance >= eps.
    
    [WaveArea, Fetch_dist, cosang,v_points] = fetch_VisiLibity(shorelines,eps,epsilon,snap_distance);
    my_fetch = {};
    for POI = 1:length(x)
    
        % fetch distance at position POI to each visible point
        my_fetch{POI}(:) = sqrt((x(POI)-v_points{1}{POI}(:,1)).^2 + (y(POI)-v_points{1}{POI}(:,2)).^2); % distance between observer and observable point on another coastline
        
    
        for vis = 1:length(v_points{1}{POI}(:,1))
            
            if my_fetch{POI}(vis)<1
                my_fetch{POI}(vis) = NaN;
                v_points{1}{POI}(vis,1) = NaN;
                v_points{1}{POI}(vis,2) = NaN;
            end
    
                x_fetch = x(POI) - v_points{1}{POI}(vis,1);
                y_fetch = y(POI) - v_points{1}{POI}(vis,2);
    
                ang_fetch{POI}(vis) = atan2(y_fetch,x_fetch);
    
        end
        
        maximum_fetch(lakes,POI) = nanmax(my_fetch{POI}(:));

    end


end
