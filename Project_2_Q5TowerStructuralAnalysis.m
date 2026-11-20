%% Project 2 - Deliverable 5: Tower Geometry and Self-Weight
% Geometry: Tower Appendix (0-77.7 m); hub height is 80.4 m.
% Steel shell only. No rotor/nacelle mass or extra tower height included.
clear; clc; close all;

%% Tower dimensions from Appendix [mm]
z_mm = [0 2775 5550 8325 11100 13875 16650 19425 22200 24975 ...
    27750 30525 33300 36075 38850 41625 44400 47175 49950 52725 ...
    55500 58275 61050 63825 66600 69375 72150 74925 77700]';
D_mm = [4120 4120 4120 4120 4120 4120 4120 4120 4120 4120 ...
    4109 4050 3991 3931 3872 3813 3754 3695 3636 3577 ...
    3518 3459 3400 3341 3282 3223 3165 3106 3047]';
t_mm = [39 39 39 35 31 29 27 25 24 22 22 22 21 21 20 ...
    20 19 19 18 18 17 17 16 15 15 14 13 13 19]';

z = z_mm/1000; Do = D_mm/1000; t = t_mm/1000;
Di = Do - 2*t;
rho_steel = 7850;       % kg/m^3, ASTM A572 Grade 50
sigma_y = 345e6;        % Pa (for later analysis)
sigma_u = 450e6;        % Pa (for later analysis)
g = 9.81;               % m/s^2
H_hub = 80.4;           % m; not used to extend tower geometry
H_tower_data = z(end);

assert(numel(z)==29 && numel(Do)==29 && numel(t)==29 && ...
    all(diff(z)>0) && all(Do>0) && all(t>0) && all(Di>0) && all(t<Do/2), ...
    'Check tower geometry data.');

%% Cross-sectional properties
A = pi/4*(Do.^2 - Di.^2);           % m^2
I = pi/64*(Do.^4 - Di.^4);         % m^4, bending
J = pi/32*(Do.^4 - Di.^4);         % m^4, polar
m_per_length = rho_steel*A;       % kg/m
w_per_length = m_per_length*g;    % N/m

% Linearly interpolate diameter and thickness, then integrate shell mass.
z_fine = linspace(z(1),z(end),10001)';
Do_fine = interp1(z,Do,z_fine,'linear');
t_fine = interp1(z,t,z_fine,'linear');
Di_fine = Do_fine - 2*t_fine;
A_fine = pi/4*(Do_fine.^2 - Di_fine.^2);
I_fine = pi/64*(Do_fine.^4 - Di_fine.^4);
m_fine = rho_steel*A_fine;
m_tower = trapz(z_fine,m_fine);
W_tower = m_tower*g;
m_cumulative = cumtrapz(z_fine,m_fine);

%% Results and section properties
fprintf('Deliverable 5: Tower Analysis\n');
fprintf('Base cross-sectional area: %.6f m^2\n', A(1));
fprintf('Base moment of inertia: %.6f m^4\n', I(1));
fprintf('Tower shell mass: %.1f kg\n', m_tower);
fprintf('Tower self-weight: %.1f kN\n', W_tower/1000);
fprintf('Note: Tower mass excludes rotor, nacelle, and other components.\n');

stationTable = table(z,Do,t,Di,A,I,J,m_per_length,w_per_length, ...
    'VariableNames',{'Height_m','OuterDiameter_m','WallThickness_m', ...
    'InnerDiameter_m','Area_m2','BendingI_m4','PolarJ_m4', ...
    'MassPerLength_kg_m','WeightPerLength_N_m'});
    openvar('stationTable')

%% Figures
figTitles = {'Tower Outer Diameter vs. Height','Tower Wall Thickness vs. Height', ...
    'Tower Bending Moment of Inertia vs. Height','Tower Cumulative Steel Shell Mass'};
yLabels = {'Outer diameter, D_o (m)','Wall thickness, t (mm)', ...
    'Area moment of inertia, I (m^4)', ...
    'Cumulative shell mass from base (metric tons)'};
xData = {z,z,z_fine,z_fine};
yData = {Do,t*1000,I_fine,m_cumulative/1000};
for k = 1:4
    figure(k); clf;
    if k<=2
        plot(xData{k},yData{k},'o-','LineWidth',1.6,'MarkerSize',4);
    elseif k==3
        plot(z_fine,I_fine,'LineWidth',1.8); hold on;
        plot(z,I,'o','MarkerSize',4); hold off;
        legend('Interpolated sections','Appendix stations','Location','best');
    else
        plot(xData{k},yData{k},'LineWidth',1.8);
    end
    grid on; xlabel('Height above base, z (m)');
    ylabel(yLabels{k}); title(figTitles{k});
end
