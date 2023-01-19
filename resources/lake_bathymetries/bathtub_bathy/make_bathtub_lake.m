function [Xmesh,Ymesh,zDep] = make_bathtub_lake(bath_slope,max_fetch_pts)


x = max_fetch_pts(:,1);
y = max_fetch_pts(:,2);

figure;
plot(x,y)
hold on;

% 101 x 101 grid within lake basin
start_x = min(x)-1e3; % add 1e3 buffer around points 
end_x = max(x) + 1e3; % add 1e3 buffer around points 
x_frac_step = (end_x - start_x)/100; 

start_y = min(y)-1e3; % add 1e3 buffer around points 
end_y = max(y)+1e3; % add 1e3 buffer around points 
y_frac_step = (end_y - start_y)/100;

x_extent = start_x:x_frac_step:end_x; 
y_extent = start_y:y_frac_step:end_y;

[X,Y] = meshgrid(x_extent,y_extent);

polyin = polyshape([x y]); % make buffer around lake coordinates so fully captured by mesh?
polyout = polybuffer(polyin,2e3);

[in,on] = inpolygon(X,Y,polyout.Vertices(:,1),polyout.Vertices(:,2));

plot(X(~in),Y(~in),'ro')
plot(X(in|on),Y(in|on),'go')

X_lake = X(in|on); % includes points within the boundary and on the boundary
Y_lake = Y(in|on); % includes points within the boundary and on the boundary


% finding the minimum fetch distance of each point within the basin to assign it to a depth
for j = 1:length(X_lake)

    for i = 1:length(x)
        xv = x(i) - X_lake(j);
        yv = y(i) - Y_lake(j);
        dv(i) = sqrt((xv).^2 + (yv).^2);

    end
    [min_dis(j),imin(j)] = nanmin(dv);
end

[max_distance,imax] = max(min_dis); % find what point is furthest from shore

figure;
scatter(X_lake,Y_lake,[],min_dis')
colorbar;
hold on
%scatter(X_lake(imax),Y_lake(imax),[],'k','filled')
title('distance from shore')



min_dep = (bath_slope).*min_dis; % assume constant slope on all sides

xLon = linspace(min(X_lake), max(X_lake), 1E+3); % for mesh
yLat = linspace(min(Y_lake), max(Y_lake), 1E+3); % for mesh
[Xmesh,Ymesh] = meshgrid(xLon, yLat); % make mesh for interpolation
zDep = griddata(X_lake, Y_lake, min_dep, Xmesh, Ymesh,'linear'); % interpolates between points

[in,on] =inpolygon(Xmesh,Ymesh,x,y); % finds all points within or along the boundary of the basin
zDep(~in&~on) = NaN; % assigns all depth outside the basin to NaN

figure;
mesh(Xmesh, Ymesh, zDep);
hold on;
[Mcon,Ccon] = contour3(Xmesh, Ymesh, zDep,[0 0],'k', 'LineWidth',2,'ShowText',1,'LabelSpacing',2000);
h = plot(x,y,'--k','LineWidth',2);
z = get(h,'ZData');
set(h,'ZData',z+10) 
view(0,90)  % XYplane view
title('Bathtub Model for Lake Depth')
colorbar;
    
         


end