function [T,H0,L0,C0] = make_wave_COE(u_pos,max_fetch,d,rho_a,rho_w,kin_visc,planet,plotting_on)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Objective:
%   This function will grow linear Airy waves for the single period at the precipice of 
%   breaking. The only source of dissipation of energy is viscous dissipation. The 
%   results can then be used in shoal_wave_COE to find the shoaling wave properties.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% inputs:
%   u_pos = near surface wind speeds [m/s]
%   max_fetch = largest fetch of basin of interest (chosen to allow waves to fully mature) [m]
%   d = water depth [m]
%   rho_a = atmospheric density [kg/m3]
%   rho_w = liquid density [kg/m3]
%   kin_visc = kinematic viscocity of liquid [m2/s]
%   planet = 'Titan' or 'Earth'
%   plotting_on = plot out results of wave maturation (yes = 1)
% outputs:
%   T = period of wave at distance of max_fetch [s]
%   H0 = wave height at distance of max_fetch [m]
%   L0 = wavelength at distance of max_fetch [m]
%   C0 = wave phase speed at distance of max_fetch [m]
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

disp('make_wave_coe: initiated...')

near_shore = length(d); % region where wave will begin to shoal (no more energy transfer from wind)

if max_fetch <= near_shore
    disp('WARNING -- make_wave_coe: max_fetch too small. skipping this point.')
    T = [NaN NaN];
    H0 = [NaN NaN];
    L0 = [NaN NaN];
    C0 = [NaN NaN];
    return
end
cd = 0.002; %drag coef (dimensionless) (lorenz and Hayes 2012)

b = 1;
min_fetch = 1;
if planet == 'Titan'
    g = 1.352; % titan gravity
elseif planet == 'Earth'
    g = 9.81;
else
    disp('make_wave_coe: ERROR -- planet name is not correct. Either use Earth or Titan')
end

% Waves growing iteratively in deepwater
for season = 1:length(u_pos)
    
x_meter = 2;

dE_dx = 0;
E_deep = 0;
h_deep = 0;
L_deep = 0;
c_deep = 0;

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

fprintf('make_wave_COE: %.2f m/s complete...\n',u_pos(season))

end

% clean out any data from wave field smaller than the smallest fetch
H0_each(:,1:min_fetch-1) = [];
T_each(:,1:min_fetch-1) = [];
L_each(:,1:min_fetch-1) = [];
c_each(:,1:min_fetch-1) = [];

if plotting_on == 1
% Uncomment to see deepwater wave growth over fetch:
    figure;
    plot(1:length(T_each(1,:)),T_each(1,:),'LineWidth',4,'Color','r','LineStyle','--')
    hold on
    plot(1:length(H0_each(1,:)),H0_each(1,:),'LineWidth',4,'Color','b','LineStyle','--')
    plot(1:length(L_each(1,:)),L_each(1,:),'LineWidth',4,'Color','g','LineStyle','--')
    plot(1:length(c_each(1,:)),c_each(1,:),'LineWidth',4,'Color','m','LineStyle','--')
    hold off;
    legend('T','H0','L','C','location','westoutside')
    xlabel('fetch')
    title('wave growth low speed')
    
    figure;
    plot(1:length(T_each(2,:)),T_each(2,:),'LineWidth',4,'Color','r','LineStyle','--')
    hold on
    plot(1:length(H0_each(2,:)),H0_each(2,:),'LineWidth',4,'Color','b','LineStyle','--')
    plot(1:length(L_each(2,:)),L_each(2,:),'LineWidth',4,'Color','g','LineStyle','--')
    plot(1:length(c_each(2,:)),c_each(2,:),'LineWidth',4,'Color','m','LineStyle','--')
    hold off;
    xlabel('fetch')
    legend('T','H0','L','C','location','westoutside')
    title('wave growth high speed')
end

% clean out data where wave field will begin to shoal
H0_winter = H0_each(1,end-near_shore); 
H0_summer = H0_each(2,end-near_shore);

T_winter = T_each(1,end-near_shore);
T_summer = T_each(2,end-near_shore);

L0_winter = L_each(1,end-near_shore);
L0_summer = L_each(2,end-near_shore);

c0_winter = c_each(1,end-near_shore);
c0_summer = c_each(2,end-near_shore);

T = [T_winter T_summer];
H0 = [H0_winter H0_summer];
L0 = [L0_winter L0_summer];
C0 = [c0_winter c0_summer];

disp('make_wave_COE: completed...')

end