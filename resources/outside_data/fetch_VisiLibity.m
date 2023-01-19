 function [WaveArea, FetchArea,cosang,v_points] = fetch_VisiLibity(shorelines,eps,epsilon,snap_distance)

% initialize
WaveArea = cell(length(shorelines),1);
% Fetch_dist = cell(length(shorelines),1);
% cosang = cell(length(shorelines),1);

% find all fetch for each point
for k = 1:length(shorelines)
    for l = 1:length(shorelines{k})
        
        % find the shoreline point, l
        x_l = shorelines{k}(l,1);
        y_l = shorelines{k}(l,2);

        
        % Get the next and previous points.
        if l > 1
            l_ind_prev = l-1;
        else
            l_ind_prev = length(shorelines{k});
        end
        
        if l < length(shorelines{k})
            l_ind_next = l+1;
        else
            l_ind_next = 1;
        end
        
        xp = shorelines{k}(l_ind_prev, 1);
        yp = shorelines{k}(l_ind_prev, 2);
        xn = shorelines{k}(l_ind_next, 1);
        yn = shorelines{k}(l_ind_next, 2);
        
        % Get the observer position (a point eps in from the shoreline because VisiLibity doesn't let you use the boundary).
        if k == 1 % this is the lake
            [Pobs,nbi] = Pint([xn,yn],[x_l,y_l],[xp,yp],eps);
        else % this is an island
            [Pobs,nbi] = Pint([xp,yp],[x_l,y_l],[xn,yn],eps);
        end
        
        V = visibility_polygon(Pobs, shorelines, epsilon, snap_distance);
        v_points{k,1}{l,1} = V;

        % fetch & wave area & distance & cos(theta-phi)
        FetchArea{k,1}(l,1) = polyarea(V(:,1),V(:,2));
        Fetch_dist{k,1}{l,1} = sqrt(sum(([shorelines{k}(l,1),shorelines{k}(l,2)] - V).^2,2));
        
        % do not let POS = OBS
        i_zero = find(Fetch_dist{1}{k}<1);
        Fetch_dist{1}{k}(i_zero) = NaN;
        v_points{1}{k}(i_zero,1) = NaN;
        v_points{1}{k}(i_zero,2) = NaN;
        
        minFetch_dist = min(Fetch_dist{k,1}{l,1},200); % this is how we can limit fetch eventually
        % Wave weighting = (F)*cosang
        weighted_fetch_dist{k,1}{l,1} = ([shorelines{k}(l,1),shorelines{k}(l,2)] - V)*nbi'; % magnitude of fetch * magnitude of normal vector * cosang
        % cos(theta - phi) = dot product of slvec and losvec
        cosang{k,1}{l,1} = -weighted_fetch_dist{k,1}{l,1}./Fetch_dist{k,1}{l,1}; %[mag_fetch*mag_norm(which is 1)*cosang]/mag_fetch = cosang
        cosang{k,1}{l,1}(isnan(cosang{k,1}{l,1})) = 0;
%         cosang{k,1}{l,1}(cosang{k,1}{l,1}<-1e-7) = 0; % any less than 0 are artifacts of discritezation and not physical
        
        % remove fetch distances that are artifacts and recalculate fetch
        % area and fetch dist and cosang
        V(cosang{k,1}{l,1}<0,:) = [];
        FetchArea{k,1}(l,1) = polyarea(V(:,1),V(:,2));
        Fetch_dist{k,1}{l,1} = sqrt(sum(([shorelines{k}(l,1),shorelines{k}(l,2)] - V).^2,2));
        % Wave weighting = (F)*cosang
        weighted_fetch_dist{k,1}{l,1} = ([shorelines{k}(l,1),shorelines{k}(l,2)] - V)*nbi'; % magnitude of fetch * magnitude of normal vector * cosang
        % cos(theta - phi) = dot product of slvec and losvec
        cosang{k,1}{l,1} = -weighted_fetch_dist{k,1}{l,1}./Fetch_dist{k,1}{l,1}; %[mag_fetch*mag_norm(which is 1)*cosang]/mag_fetch = cosang
        cosang{k,1}{l,1}(isnan(cosang{k,1}{l,1})) = 0;
        
        Wavepts = [shorelines{k}(l,1),shorelines{k}(l,2)]+(V-[shorelines{k}(l,1),shorelines{k}(l,2)]).*cosang{k,1}{l,1};
        WaveArea{k,1}(l,1) = polyarea(Wavepts(:,1),Wavepts(:,2));
        

    end
end
end
