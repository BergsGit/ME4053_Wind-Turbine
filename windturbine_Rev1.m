clear; clc; close all;
%ME4053 Wind Turbine Project
%Group 26
%
%Authors:
%...
%%
%_____Table of Contents____
% 1. Deliverables
%       1.
%       2.
%       3.
%       4.
%       5.
% 2. Functions
%       Wind Turbine Analysis
%           Import Data and Paramaters
%           Blade Analysis                  (... - ...)
%
%

%% 1. Deliverables

% Load blade geometry and airfoil lookup tables
bladeData     = readtable('BladeProfile.csv');
airfoils      = LoadAirfoils();

% DELIVERABLE 1: V = 10 m/s, Omega = 14 RPM, Pitch = 0 deg
    V_part1     = 10.0;                  % [m/s]
    omega_part1 = 14.0 * pi / 30;        % [rad/s]
    Pitch_Angle1 = 0.0;                  % [deg]

    [P_part1, T_part1, Cp_part1, Ct_part1] = solveBEMs(Pitch_Angle1,...
        V_part1, omega_part1, bladeData, airfoils);

    fprintf('=== Deliverable 1 Results ===\n');
    fprintf('Coefficient of Power (Cp):  %.4f\n', Cp_part1);
    fprintf('Coefficient of Thrust (Ct): %.4f\n\n', Ct_part1);

%  DELIVERABLE 2: Blade Pitch Control at V = 8.5 m/s, TipSpeedRatio = 7.10
    U_wind2  = 8.5;                 % [m/s] Assigned wind speed
    TSR      = 7.10;                % [] Tip Speed Ratio
    omega2   = TSR * U_wind2 / 48;  % [rad/s] Rotational speed from TSR

    % Run the pitch sweep function (-10 to 25 deg in 0.05 deg increments)
    theta_range = -10:0.05:20;
    [~, ~, Cp_sweep, ~] = sweepPitch(theta_range, ...
        U_wind2, omega2, bladeData, airfoils);

    %Find maximum Cp value and its corresponding pitch angle
    [max_Cp, index] = max(Cp_sweep);
    maxCp_angle = theta_range(index);

    fprintf('=== Deliverable 2 Results ===\n');
    fprintf('Max Coefficient of Power:  %.4f\n', max_Cp);
    fprintf('Pitch Angle at Max C_p:    %.4f\n\n', maxCp_angle);

    %Plot figure
    figure(1)
    plot(theta_range, Cp_sweep);
    grid on;
    title('C_p vs. Pitch Angle')
    xlabel('Pitch Angle (degrees)')
    ylabel('C_p')
    set(gca, 'FontName', 'Times New Roman', 'FontSize', 12);

%  DELIVERABLE 3: Blade Pitch Control at V = 8.3 m/s
    U_wind3  = 8.3;                  % [m/s] Assigned wind speed
    TSR      = 5.85:0.05:9.35;       % [] Tip Speed Ratio 
    omega3   = TSR .* U_wind2 / 48;  % [rad/s] Rotational speed from TSR

    % Run the pitch sweep function (-10 to 25 deg in 0.1 deg increments)
    theta_range = -10:0.1:20;
    
    %Array to store Cp data in
    Cp_data = zeros(numel(omega3), numel(theta_range));
    
    %sweep for each TSR value
    for i = 1:numel(omega3)
        %Sweep each pitch angle value
        [~, ~, Cp_sweep, ~] = sweepPitch(theta_range, ...
            U_wind3, omega3(i), bladeData, airfoils);
        %Store data
        Cp_data(i,:) = Cp_sweep;
    end

    %Find maximum Cp value and its corresponding pitch angle and TSR
    max_Cp = max(Cp_data(:));
    [row, col] = find(Cp_data == max_Cp);
    
    optimal_TSR        = TSR(row);
    optimal_PitchAngle = theta_range(col);

    fprintf('=== Deliverable 3 Results ===\n');
    fprintf('Max Coefficient of Power:     %.4f\n', max_Cp);
    fprintf('Pitch Angle at Max C_p:       %.4f\n', optimal_PitchAngle);
    fprintf('Tip Speed Ratio at Max C_p:   %.4f\n\n', optimal_TSR);
    
    %Plot Figure
    figure(3)
    contourf(theta_range, TSR, Cp_data);
    title('C_p vs. Pitch Angle & Tip Speed Ratio')
    xlabel('Pitch Angle (degrees)')
    ylabel('Tip Speed Ratio')
    c = colorbar;
    c.Label.String = 'C_p';
    set(gca, 'FontName', 'Times New Roman', 'FontSize', 12);

%  DELIVERABLE 4: Blade Pitch Control at V = 19 m/s, Omega = 15.5 RPM
    V_wind   = 19.0;               % [m/s] Assigned wind speed
    rpm      = 15.5;               % [RPM] Maximum rotational speed
    omega    = rpm * pi / 30;      % [rad/s]
    P_target = 2.5e6;              % [W] Rated power limit (2.5 MW)

    % Run the pitch sweep function (0 to 25 deg in 0.5 deg increments)
    theta_range = -7.5:0.5:25;
    [P_sweep, T_sweep, Cp_sweep, Ct_sweep] = sweepPitch(theta_range, ...
        V_wind, omega, bladeData, airfoils);

    % Pinpoint the exact pitch angle where Power = 2.5 MW using fzero
    objFun      = @(theta_p) solveBEMs(theta_p, V_wind, omega, ...
        bladeData, airfoils) - P_target;
    theta_exact = fzero(objFun, [0, 25]);

    % Evaluate final operating state and span wise loads
    [P_final, T_final, Cp_final, Ct_final, radialLoads] = ...
        solveBEMs(theta_exact, V_wind, omega, bladeData, airfoils);

    fprintf('=== Deliverable 4 Results ===\n');
    fprintf('Exact Pitch Angle (theta_p): %.3f deg\n', theta_exact);
    fprintf('Aerodynamic Power:           %.4f MW\n', P_final / 1e6);
    fprintf('Rotor Thrust:                %.2f kN\n', T_final / 1e3);
    fprintf('Power Coefficient (Cp):      %.4f\n', Cp_final);
    fprintf('Thrust Coefficient (Ct):     %.4f\n', Ct_final);

% DELIVERABLE 5: 
%% 2. Functions

%_____Import Functions_____
function airfoils = LoadAirfoils()
% Function: LoadAirfoils
% Input:    none
% Output:   airfoils (struct containing airfoil tables)
% Process:  For each unique airfoil, the corrseponding csv is read in and
%           stored in the data struct under its code.
    airfoils.DU97 = readtable('DU97-W-300.csv');
    airfoils.DU91 = readtable('DU91-W2-250.csv');
    airfoils.DU93 = readtable('DU93-W-210.csv');
    airfoils.DU96 = readtable('DU96-W-180.csv');
end

function afTable = getAirfoilTable(afName, airfoils)
% Function: getAirfoilTable
% Input:    afName (airfoil name), airfoils (struct containing airfoil
%           tables)
% Output:   afTable (table associated with the inputted airfoil name)
% Process:  Switch statement checks the airfoil name and the matching
%           airfoil is returned.
    switch afName
        case "DU 97-W-300"
            afTable = airfoils.DU97;
        case "DU 91-W2-250"
            afTable = airfoils.DU91;
        case "DU 93-W-210"
            afTable = airfoils.DU93;
        case "DU 96-W-180"
            afTable = airfoils.DU96;
        otherwise
            error('Unknown airfoil name: %s', afName);
    end
end
%

%_____Wind Turbine Blade Analysis_____
function [P_total, T_total, Cp, Ct, radialLoads] = solveBEMs(theta_p, V_wind, omega, bladeData, airfoils)
% Computes rotor power, thrust, Cp, Ct, and span wise loads
% Function: solveBEMs
% Input:    theta_p (pitch angle), V_wind (wind speed), omega (rotatianal 
%           velocity), bladeData (struct containing blade paramaters), 
%           airfoils (struct containing airfoil tables)
% Output:   P_total (output power value), T_total (array torque value),
%           Cp (coefficient of power), Ct (coefficient of torque),
%           radialLoads (array containing radial load at each r value)
% Process:  First, tower data, paramaters, and air properties are defined.
%           Next, the Tip speed ratio is found and used to define the local
%           tip speed ratio which in turn is used to find the axial'
%           induction factor for each radial position on the blade. The 
%           relative velocity is then found throuh the vector sum of the 
%           converted radial velocity and windspeed. This velocity can then
%           be used to find the reynolds number and AoA.
%
%           The arrays for the output paramaters are then initizialized to
%           the size of the 'r' array. And their values found through
%           analysis and linear interpolation of the imported tables.
%           trapz() is then used to find the sum of these arrays and the
%           final output paramaters are calculated.
%
    % Geometry and fluid properties from reference script
    r     = bladeData.DistanceFromCenterOfRotation / 1000;  % [m] (43x1)
    c     = bladeData.ChordLength / 1000;                   % [m] (43x1)
    beta  = bladeData.BladeTwist;                           % [deg] (43x1)
    R     = max(r);                                         % [m] 48 m
    B     = 3;                                   % Number of blades
    A_rot = (R^2) * pi;                          % Swept area [m^2]
    rho   = 1.1;                                 % Air density [kg/m^3]
    mu    = 1.8e-5;                              % Dynamic viscosity [Pa*s]
    a     = 1 / 3;                               % Axial induction factor

    % Local tip speed ratio and angular induction factor
    TSR   = R * omega / V_wind;
    TSR_r = TSR .* r / R;
    ap    = -0.5 + 0.5 .* sqrt(1 + (4 ./ (TSR_r.^2)) .* a .* (1 - a));

    % Relative velocity components and inflow angle
    V_x   = V_wind .* (1 - a);
    V_y   = r .* omega .* (1 + ap);
    V_rel = sqrt(V_y.^2 + V_x.^2);
    phi   = 90 - rad2deg(atan(V_y ./ V_x));                 % [deg]

    % Angle of attack and local Reynolds number
    alpha_deg = phi - (theta_p + beta);
    Re        = rho .* V_rel .* r ./ mu;

    nStations = height(bladeData);
    F_theta   = zeros(1, nStations);
    F_z       = zeros(1, nStations);
    CL_vec    = zeros(1, nStations);
    CD_vec    = zeros(1, nStations);
    FL_vec    = zeros(1, nStations);
    FD_vec    = zeros(1, nStations);

    for i = 1:nStations
        afName = strtrim(string(bladeData{i, 5}));

        if afName == "circle"
            % Uses provided cylinderCD.m function in project folder
            CD_vec(i)  = cylinderCD(Re(i));
            FD_vec(i)  = CD_vec(i) * 0.5 * rho * (V_rel(i)^2) * c(i);
            angleFnet  = 90 - phi(i);
            F_theta(i) = FD_vec(i) * sind(angleFnet);
            F_z(i)     = FD_vec(i) * cosd(angleFnet);
        else
            afTable = getAirfoilTable(afName, airfoils);
            CL_val  = interp1(afTable.AoA, afTable.CL, alpha_deg(i), 'spline');
            CD_val  = interp1(afTable.AoA, afTable.CD, alpha_deg(i), 'spline');
            [F_theta(i), F_z(i), CL_vec(i), CD_vec(i), FL_vec(i), FD_vec(i)] = ...
                calcElementForces(CL_val, CD_val, rho, V_rel(i), c(i), phi(i));
        end
    end

    % Integrate torque and thrust across span using trapz
    Q_total = B * trapz(r, F_theta' .* r);
    P_total = Q_total * omega;
    T_total = B * trapz(r, F_z);

    Ct = T_total / (0.5 * rho * (V_wind^2) * A_rot);
    Cp = P_total / (0.5 * rho * (V_wind^3) * A_rot);

    % Package radial loads for Deliverable 5
    radialLoads.r       = r;
    radialLoads.F_theta = F_theta';
    radialLoads.F_z     = F_z';
end


function [F_theta, F_z, CL, CD, FL, FD] = calcElementForces(CL, CD, rho, V_rel, c, phi)
% Function: calcElementForces
% Input:    CL (coeffient of lift), CD (coefficient of drag), rho (air
%           density), V_rel (relative air velocity), c (chord length), 
%           phi (AoA)
% Output:   F_theta (force in the theta direction), F_z (force in the z
%           direction), CL (lift coefficient), CD (drag coefficient),
%           L (lift force), FD (drag force)
% Process:  First, the lift and drag forces are found through the
%           definition of their coefficient equations. The drag and lift
%           forces are then split into their theta and z contributions, and
%           the Ftheta and Fz are computed through the vector sum.
    FL = CL * 0.5 * rho * (V_rel^2) * c;
    FD = CD * 0.5 * rho * (V_rel^2) * c;

    F_z     = (FL * cosd(phi)) + (FD * sind(phi));
    F_theta = (FL * sind(phi)) - (FD * cosd(phi));
end

function [P_vec, T_vec, Cp_vec, Ct_vec] = sweepPitch(theta_vec, V_wind, ...
    omega, bladeData, airfoils)
% Function: sweepPitch
% Input:    theta_vec (array containing pitch angles), V_wind (wind speed),
%           omega (rotatianal velocity), bladeData (struct containing blade
%           paramaters), airfoils (struct containing airfoil tables)
% Output:   P_vec (array containing output power values), T_vec (array
%           containing torque values), Cp_vec (array containing Cp values),
%           Ct_vec (array containing Ct values)
% Process:  First arrays for each of the output variables are initialized
%           each with same size as the inputted pitch angles. These arrays
%           are then populated with their corresponding values from the
%           solveBEMS function() and returned.
    nSteps = length(theta_vec);
    P_vec  = zeros(1, nSteps);
    T_vec  = zeros(1, nSteps);
    Cp_vec = zeros(1, nSteps);
    Ct_vec = zeros(1, nSteps);

    for k = 1:nSteps
        [P_vec(k), T_vec(k), Cp_vec(k), Ct_vec(k)] = ...
            solveBEMs(theta_vec(k), V_wind, omega, bladeData, airfoils);
    end
end
%