clc
clear
close all
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Author: Una Schneck (schneck.una@gmail.com)
% 
% This code will produce plots for the maximum entrainment depth (depth where the grains' Shield parameter crosses the critical Shields parameter as  defined 
% by Komar & Miller 1973. Code will produce maps for constant  V-shaped basin bathymetries of Ontario Lacus (Lake 1) and Ligeia Mare  (Lake 2). 
%
% The maximum fetches for the lakes is computed in a seperate doc because  takes a while to run since using Visilibity (get_max_fetches.m). The  coordinates 
% were pulled from polar composites of RADAR swathes and  linearly interpolated so the points are 1 pixel (351 m) apart.
%--------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
% External function requirements:
%   (1) figure_settings.m -- makes plots pretty
%--------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
% Code produces:
%   (1) Wave growth plots
%   (2) Wave shoaling plots
%   (3) Depth entrainment vs. D50 plots
%   (4) Contour plots for depth entrainment
%---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
% %Aprox time to run: The larger lakes (Ligeia Mare) can take ~3 hours for all the points
%---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
% Cite: The Shoreline Shaping Capability of Waves and Tides at Titan's Lakes (Schneck et al.)
% 
% --------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% --------------------------------------------------- 1 . PREPPING THE RUN ------------------------------------------------ %% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

tic

% CONSTANTS
g = 1.352; % titan gravity
%EXPERIMENTAL CONSTANTS (from Earth experiments)
cd = 0.002; %drag coef (dimensionless) (lorenz and Hayes 2012)
a = 0.6; % pore space correction for medium sand
%K = 0.77; % <<------ From earth experimental estimates
%PROPERTIES OF FLOW
rho_a = 1.2; %kg/m3 (atmospheric density)
rho_s = 940;% kg/m3 (ice grain density)     
rho_w = 590;% kg/m3 (liquid hydrocarbon density) Hayes 2016
d50 = [6.35e-5:1e-4:0.1];% diameters for [finegrain to 10 cm] [m]
kin_visc = 7.5e-7; %m2/s % Hayes 2016

% grains of interest are the smallest (fine-sand) and the biggest (gravel)
grainsizeindex = [1 length(d50)];

% Individual Lake Properties (from RADAR altimetry observations)
pos_bath_slope = [0.5e-3 2e-3]; % total range measured between OL and LM  
pos_max_depth = [50 160]; % [max(OL) max(LM)]

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% --------------------------------------------------- 2 . LAKE AND BATHY TO RUN ------------------------------------------------ %% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for bath_type = 2 % [shallow steep]


    for lakes = 2 % [OL LM]
        %% EXPERIMENTAL VARIABLES AND CONSTANTS %%
        bath_slope = pos_bath_slope(bath_type); % bathymetric slope of lake
        max_depth_lake = pos_max_depth(lakes); % maximum depth of lake
        d = max_depth_lake:-bath_slope:10^-4; % waterdepth approaching shore; % water depth [m], deep water needs to be (d/L > 1/2), slope from Hayes 2010
        near_shore = length(d); % size of the sloping bed's array
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% --------------------------------------------------- 3 . MAKE WAVES ------------------------------------------------ %% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


        %% CREATE WAVE USING ENERGY CONSERVATION %%
        
        % get maximum fetch at each point in a seperate doc (get_max_fetch_lakes.m) to save runtime
        addpath('..\resources\max_fetch')
        if lakes == 1
            load('ol_max_fetch.mat','maximum_fetch')
        elseif lakes == 2
            load('lm_max_fetch.mat','maximum_fetch')
        else
            disp('no lake data')
        end


        for ff = 1:1000:length(maximum_fetch) % different max fetch for each point on shoreline
            max_fetch = maximum_fetch(ff);
            
            % [summer winter] wind
            u_pos = [0.4 3.3];
            
            b = 1;
            
            % Waves growing iteratively in deepwater
            for season = 1:2
                
            x_meter = 2;
            
            dE_dx = 0;
            E_deep = 0;
            h_deep = 0;
            c_deep = 0;
            L_deep = 0;
            T_deep = 0;
        
        
                while x_meter <= max_fetch
            
                    c_deep(1,:) = 0;
                    E_deep(1,:) = 0;
                    dt(1,:) = 0;
                    visc_tau(1,:) = 1;
            
                    dE_dx(x_meter,:) = (0.5*cd*rho_a.*((u_pos(season) - c_deep(x_meter-1,:)).^2))*exp(-dt(x_meter-1,:)/visc_tau(x_meter-1,:));
            
                    E_deep(x_meter,:) = E_deep(x_meter-1,:) + dE_dx(x_meter,:);
            
                    h_deep(x_meter,:) = sqrt((8*E_deep(x_meter,:))/(rho_w*g));
            
                    L_deep(x_meter,:) = 7*h_deep(x_meter,:);
            
                    T_deep(x_meter,:) = sqrt((2*pi*L_deep(x_meter,:))/g);
                    
                    c_deep(x_meter,:) =L_deep(x_meter,:)/T_deep(x_meter,:);
            
                     dt(x_meter,:) = (2*x_meter)/c_deep(x_meter,:);
                     visc_tau(x_meter,:) = ((L_deep(x_meter,:)^2)*rho_w)/(8*pi^2*kin_visc);
            
                   
                    x_meter = x_meter+1;

                end
            

            H0_each(b,:) = h_deep;
            T_each(b,:) = T_deep;
            L_each(b,:) = L_deep;
            c_each(b,:) = c_deep;

            b = b + 1;

                if ff == 1
                    figure;    
                    plot(2:x_meter,h_deep,'--r','LineWidth',5)
                    hold on;
                    plot(2:x_meter,L_deep,'--g','LineWidth',5)
                    plot(2:x_meter,T_deep,'--b','LineWidth',5)
                    plot(2:x_meter,c_deep,'--m','LineWidth',5)
                    legend('h','L','T','c')
                    title('deepwater')
                    grid on;
                    xlabel('fetch [m]')
                    cap_d_deep = sprintf('Slope: %.4f, Season: %i, Lake %i',bath_slope,season,lakes);
                    title(cap_d_deep)
                end
            end

            
            if max_fetch>near_shore % if wave matures before shoaling
                H0_winter = H0_each(1,end-near_shore); %near_shore = length(d) in loop approaching shore
                H0_summer = H0_each(2,end-near_shore);
                
                T_winter = T_each(1,end-near_shore);
                T_summer = T_each(2,end-near_shore);
                
                L0_winter = L_each(1,end-near_shore);
                L0_summer = L_each(2,end-near_shore);
                
                c0_winter = c_each(1,end-near_shore);
                c0_summer = c_each(2,end-near_shore);
            else % if wave doesn't mature before shoaling
                H0_winter = H0_each(1,:); 
                H0_summer = H0_each(2,:);
                
                T_winter = T_each(1,:);
                T_summer = T_each(2,:);
                
                L0_winter = L_each(1,:);
                L0_summer = L_each(2,:);
                
                c0_winter = c_each(1,:);
                c0_summer = c_each(2,:);
            end
    
            H0_each(:) = [];
            T_each(:) = [];
            L_each(:) = [];
            c_each(:) = [];
            L = [];
            H = [];
            d0 = [];
            um = [];
        
            T = [T_winter T_summer];
            H0 = [H0_winter H0_summer];
            L0 = [L0_winter L0_summer];
            C0 = [c0_winter c0_summer];
            
        end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% --------------------------------------------------- 4. SHOAL WAVES ------------------------------------------------ %% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        for season = 1:2
        
            Cg0 = 0.5*C0(season); %CEM
            
            % SHOAL
            
            alpha_0 = deg2rad(0);
            
            L = [];
            alpha = [];
            C = [];
            n = [];
            Cg = [];
            KR = [];
            KS = [];
            H = [];

        
            %         shoaling wave properties (T = constant)
                for j = 1:length(d)
                    L(j) = L0(season).*sqrt(tanh((4*(pi^2).*d(j))./((T(season)^2).*g))); %shoaling wavelength (Eckart Aproximation, CEM)
                    alpha(j) = asin(tanh((2*pi.*d(j))./L(j)).*sin(alpha_0)); %shoaling wave angle (always zero because starts at zero) Davidson-Arnott 5.31)
                    C(j) = tanh((2*pi.*d(j))./L(j)).*C0(season); %shoaling wave celerity (Davidson Arnott 5.31)
                    n(j) = 0.5.*(1 + (((4*pi.*d(j))./L(j))./(sinh((4*pi.*d(j))./L(j))))); % shoaling coefficient (TABLE CEM)
                    Cg(j) = n(j).*C(j); % shoaling wave group velocity (TABLE CEM)
                    KR(j) = sqrt(cos(alpha_0)/cos(alpha(j))); %refraction coef  http://www.coastalwiki.org/wiki/Shallow-water_wave_theory
                    KS(j) = sqrt(Cg0/Cg(j)); % shoaling coef http://www.coastalwiki.org/wiki/Shallow-water_wave_theory
                    if cos(alpha_0) < 0 
                        H(j) = NaN;
                    else
                        H(j) = KR(j)*KS(j)*H0(season); % http://www.coastalwiki.org/wiki/Shallow-water_wave_theory
                    end
                    
        
                end

                figure;
                plot(d,L,'--r','LineWidth',5)
                hold on
                plot(d,C,'--g','LineWidth',5)
                plot(d,n,'--b','LineWidth',5)
                plot(d,Cg,'--m','LineWidth',5)
                plot(d,H,'--c','LineWidth',5)
                legend('L','C','n','Cg','H')
                title('shoaling')
                set(gca, 'XScale', 'log')
                set(gca,'Xdir','reverse')
                grid on;
                xlabel('depth [m]')
                cap_d_deep = sprintf('Slope: %.4f, Season: %i, Lake %i',bath_slope,season,lakes);
                title(cap_d_deep)

                d0(season,:) = H./sinh((2*pi.*d)./L); % wave orbital diameter
                um(season,:) = (pi.*d0(season,:))/T(season); % wave orbital velocity

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% --------------------------------------------------- 5. SHIELDS NUMBER ------------------------------------------------ %% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

                for kk = 1:length(d50)
                    
                    
                    shields{season,kk} = rho_w.*(um(season,:).^2)./((rho_s - rho_w)*g.*d50(kk));
                    % check what depth j the shield value exceeds k&m shields value for each
                    % d50 for each season i
                    KM_crash{season,kk} = 0.3.*sqrt(d0(season,:)./d50(kk));
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% --------------------------------------------------- 6 . ENTRAINMENT DEPTH ------------------------------------------------ %% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

                    if isempty(d(find(shields{season,kk}>KM_crash{season,kk},1,'first'))) % does not reach entrainment
                        d_crash(season,kk) = NaN;
                    else % does reach entrainment
                        d_crash(season,kk) = d(find(shields{season,kk}>KM_crash{season,kk},1,'first'));
                    end
                end
        end

        figure;
        plot(d50,d_crash(1,:),'--','LineWidth',5)
        hold on
        plot(d50,d_crash(2,:),'--','LineWidth',5)
        set(gca, 'YScale', 'log')
        set(gca, 'XScale', 'log')
        xlabel('D_{50}')
        ylabel('entrainment depth [m]')
        cap_d_d50 = sprintf('Slope: %.4f, Season: %i, Lake %i',bath_slope,season,lakes);
        title(cap_d_d50)
        grid on;
        box on;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% --------------------------------------------------- 7 . ENTRAINMENT DEPTH PLOTS ------------------------------------------------ %% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

            figure;
            for my_grain_size = 1:2
                grainsize = d50(grainsizeindex(my_grain_size));
        
                 
                if season == 2 && bath_type == 1 && lakes == 1 && ff == 1 || season == 2 && bath_type == 1 && lakes == 2 && ff == 1
                    figure();
                    plot(d50,d_crash(1,:),'linewidth',5)
                    hold on
                    plot(d50,d_crash(2,:),'linewidth',5)
                    legend('0.4 m/s','3.3 m/s')
                    xlabel('d50 m')
                    ylabel('depth of entrainment m')
                    set(gca, 'YScale', 'log')
                    set(gca, 'XScale', 'log')
                    grid on;
                    cap = sprintf('Slope: %.4f, Lake: %i',bath_slope,lakes);
                    title(cap)
                end
        
                % create basin maps for entrainment distance
                addpath('..\resourceslake_bathymetries\bathtub_bathy')
                if lakes == 1
                    load('ol_bathtub_0.000500_slope.mat','x','y')
                elseif lakes == 2
                    load('lm_bathtub_0.000500_slope.mat','x','y') 
                else
                    disp('no lake data')
                end
              
                
                for i = 1:length(x)
                
                    if i == 1
                        x4(i) = (((y(i+1) - y(end))));
                        y4(i) = (((x(end) - x(i+1))));
                    elseif i == length(x)
                        x4(i) = (((y(1) - y(i-1))));
                        y4(i) = (((x(i-1) - x(1))));
                    else
                        x4(i) = (((y(i+1) - y(i-1))));
                        y4(i) = (((x(i-1) - x(i+1))));
                    end
                    v = [x4(i) y4(i)];
                    v = v./norm(v);
                    u1(i) = (d_crash(season,grainsizeindex(my_grain_size))/bath_slope)*v(1); % entrainment distance (x-commponent)
                    v1(i) = (d_crash(season,grainsizeindex(my_grain_size))/bath_slope)*v(2); % entrainment distance (y-component)
                    xinner(i) = x(i) - (d_crash(season,grainsizeindex(1))/bath_slope)*v(1);
                    yinner(i) = y(i) - (d_crash(season,grainsizeindex(1))/bath_slope)*v(2);
                    xouter(i) = x(i) - (d_crash(season,grainsizeindex(end))/bath_slope)*v(1);
                    youter(i) = y(i) - (d_crash(season,grainsizeindex(end))/bath_slope)*v(2);
                
                end
                
                
                
                [in,on] = inpolygon(xinner,yinner,x,y);
                
                plot(x,y,'-k','linewidth',3)
                hold on;
                con_x = x - u1';
                con_y = y - v1';
                if my_grain_size == 1
                    plot(con_x(in), con_y(in),'Color',[62/255 150/255 81/255],'LineWidth',4)
                elseif my_grain_size == 2
                    plot(con_x(in),con_y(in),'Color',[204/255 37/255 41/255],'LineWidth',4)
                else             
                    
                end
                set(gca,'xticklabel',[])
                set(gca,'yticklabel',[])
                box on;

            end
    end
    if bath_slope == 0.5e-3 && lakes == 1
        load('ol_bathtub_0.000500_slope.mat','Xmesh','Ymesh','zDep','x','y') % from make_synthetic_bathy.m
    elseif bath_slope == 2e-3 && lakes == 1
        load('ol_bathtub_0.002000_slope.mat','Xmesh','Ymesh','zDep','x','y') % from make_synthetic_bathy.m
    elseif bath_slope == 0.5e-3 && lakes == 2
        load('lm_bathtub_0.000500_slope.mat','Xmesh','Ymesh','zDep','x','y') % from make_synthetic_bathy_lm.m
    elseif bath_slope == 2e-3 && lakes == 2
        load('lm_bathtub_0.002000_slope.mat','Xmesh','Ymesh','zDep','x','y') % from make_synthetic_bathy_lm.m
    else
        disp('no bathy map available')
    end

     plot(x,y,'k','LineWidth',4)
     hold on;
    [l1,lh1] = contour3(Xmesh, Ymesh, zDep, 0:20:max(max(abs(zDep))), 'LineStyle','--','LineWidth',1,'ShowText',1,'LabelSpacing',1200,'LineColor','k');

    view(0,90)
    set(gca,'YTickLabel',[]);
    set(gca,'XTickLabel',[]);
    box on;
    grid on;
    cap = sprintf('Slope: %.4f, Season %i',bath_slope,season);
    title(cap);
    progressbar(bath_type/length(pos_bath_slope))
end


figure
plot(x,y,'k','LineWidth',4)
hold on;
[l1,lh1] = contour3(Xmesh, Ymesh, zDep, 0:50:max(max(abs(zDep))), 'LineStyle','--','LineWidth',1,'ShowText',1,'LabelSpacing',1200,'LineColor','k');
hold on;
[M1,cc1] = contour3(Xmesh,Ymesh,zDep,[d_crash(2,1) d_crash(2,1)], 'LineStyle','-','LineWidth',4,'ShowText',0,'LabelSpacing',1200,'Color',[62/255 150/255 81/255]);
[M2,cc2] = contour3(Xmesh,Ymesh,zDep,[d_crash(2,end) d_crash(2,end)],'LineStyle','-','LineWidth',4,'ShowText',0,'LabelSpacing',1200,'Color',[204/255 37/255 41/255]);
view(0,90)
set(gca,'YTickLabel',[]);
set(gca,'XTickLabel',[]);

addpath('..\resources\helper_functions')
figure_settings(20,'Bold',1)
toc
