function [d0_breaking,um_breaking,H_breaking,L_breaking] = shoal_wave_COE(T,H0,L0,C0,alpha_0,rho_w,d,u_pos,planet,break_pt,max_fetch,plotting_on)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Objective: 
% This function will shoal waves to the breaking depth and calculate the orbital diameter and
% orbital velocity of water particles at the point of breaking.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% inputs:
%   T = wave period [s]
%   H0 = deepwater wave height [m]
%   L0 = deepwater wavelength [m]
%   C0 = deepwater phase speed [m/s]
%   alpha_0 = deepwater wave angle [rads]
%   rho_w = liquid density [kg/m3]
%   d = water depth [m]
%   u_pos = near-surface wind-speed [m/s]
%   planet = 'Earth' or 'Titan'
%   break_pt = depth where shoaling calculation stops (1 = breaking wave depth, 0 = shoreline)
%   plotting_on = plot out results of wave breaking (yes = 1)
% outputs:
%   d0_breaking = orbital diameter at breaking wave depth [m]
%   um_breaking = orbital velocity at breaking wave depth [m/s]
%   H_breaking = wave height at breaking wave depth [m]
%   L_breaking = wave length at breaking wave depth [m]
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp('shoal_wave_COE initiated...')

% choose gravity to use
if planet == 'Titan'
    g = 1.352; % titan gravity
elseif planet == 'Earth'
    g = 9.81;
else
    disp('make_wave_coe: ERROR -- planet name is not correct. Either use Earth or Titan')
end

% shoal wave 
for i = 1:length(u_pos)

    Cg0 = 0.5*C0(i); %CEM
    
    % SHOAL

    % breaking wave depth
    hb = (1.28*H0(i))/(3.3*(H0(i)/L0(i))^0.33); % min depth before breaking (Davidson-Arnott 5.34 & 5.36)
    % shoaling wave properties (T = constant)
        for j = 1:length(d)
            L(j) = L0(i).*sqrt(tanh((2*pi.*d(j))./(L0(i)))); %shoaling wavelength (Eckart Aproximation, CEM)
            alpha(j) = asin(tanh((2*pi.*d(j))./L(j)).*sin(alpha_0(i))); %shoaling wave angle (always zero because starts at zero) Davidson-Arnott 5.31)
            C(j) = tanh((2*pi.*d(j))./L(j)).*C0(i); %shoaling wave celerity (Davidson Arnott 5.31)
            n(j) = 0.5.*(1 + (((4*pi.*d(j))./L(j))./(sinh((4*pi.*d(j))./L(j))))); % shoaling coefficient (TABLE CEM)
            Cg(j) = n(j).*C(j); % shoaling wave group velocity (TABLE CEM)
            KR(j) = sqrt(cos(alpha_0(i))/cos(alpha(j))); %refraction coef  http://www.coastalwiki.org/wiki/Shallow-water_wave_theory
            KS(j) = sqrt(Cg0/Cg(j)); % shoaling coef http://www.coastalwiki.org/wiki/Shallow-water_wave_theory
            if cos(alpha_0(i)) < 0 
                H(j) = NaN;
            else
                H(j) = KR(j)*KS(j)*H0(i); % http://www.coastalwiki.org/wiki/Shallow-water_wave_theory
            end
            E(j) = (1/8)*rho_w*g.*((H(j)).^2); % wave energy
        end

    if plotting_on == 1
        figure;
        plot(d,L,'LineWidth',4)
        hold on;
        plot(d,alpha,'LineWidth',4)
        plot(d,C,'LineWidth',4)
        plot(d,H,'LineWidth',4)
        plot(d,Cg,'LineWidth',4)
        plot(d+d(1),L0(i).*ones(size(d)),'LineWidth',4)
        plot(d+d(1),alpha_0(i).*ones(size(d)),'LineWidth',4)
        plot(d+d(1),C0(i).*ones(size(d)),'LineWidth',4)
        plot(d+d(1),H0(i).*ones(size(d)),'LineWidth',4)
        plot(d+d(1),Cg0.*ones(size(d)),'LineWidth',4)
        hold off;
        legend('L','\alpha','C','H','Cg','L0','\alpha_{0}','C0','H0','Cg0','Location','eastoutside')
        cap = sprintf('wave shoaling u_{10} = %.2f',u_pos(i));
        title(cap)
        grid on;
        xlim([0 d(1)+100])
        set(gca, 'XScale', 'log')
        set(gca, 'XDir','reverse')
        xlabel('distance from shore [m]')
            
        figure; 
        plot(d,E,'LineWidth',4)
        cap = sprintf('wave energy shoaling u_{10} = %.2f',u_pos(i));
        title(cap)
        set(gca, 'XDir','reverse')
        hold off;
        set(gca, 'XScale', 'log')
        grid on;
        xlabel('distance from shore [m]')
        ylabel('Shoaling wave energy')
    end

    
    if break_pt == 1
        % find d0 and um at the empirical breaking wave depth
        if any(d<=hb) % wave breaks before shoreline
            H_breaking(i) = H(find(d<=hb,1,'first'));
            L_breaking(i) = L(find(d<=hb,1,'first'));
            d0_breaking(i) = H_breaking(i)/sinh((2*pi*hb)/L_breaking(i));
            um_breaking(i) = (pi*d0_breaking(i))/T(i);
            fprintf('Breaking wave depth for %.2f m/s is at %.1f m depth\n',u_pos(i),hb)

        else % for very energetic waves, the breaking depth can be right at the shore, setting depth to ~0 m
            disp('WARNING: BREAKING DEPTH MAY BE POORLY DEFINED. SETTING BREAKING DEPTH TO ZERO DEPTH')
            H_breaking(i) = H(end);
            L_breaking(i) = L(end);
            d0_breaking(i) = H_breaking(i)/sinh((2*pi*hb)/L_breaking(i));
            um_breaking(i) = (pi*d0_breaking(i))/T(i);
        end
        
    else % calculate wave breaking up to the shoreline
        d0_breaking(i,:) = H./sinh((2*pi.*d)./L); % wave orbital diameter
        um_breaking(i,:) = (pi.*d0_breaking(i,:))/T(i); % wave orbital velocity      
    end
    
    fprintf('shoal_wave_coe:  %.2f m/s complete...\n',u_pos(i))
end
disp('shoal_wave_COE completed...')
end