%% Wind Turbine BEM Analysis (Deliverable 1 & 4)
clear; clc; close all;

% Load blade geometry and airfoil lookup tables
bladeData     = readtable('BladeProfile.csv');
airfoils.DU97 = readtable('DU97-W-300.csv');
airfoils.DU91 = readtable('DU91-W2-250.csv');
airfoils.DU93 = readtable('DU93-W-210.csv');
airfoils.DU96 = readtable('DU96-W-180.csv');

%  DELIVERABLE 1: V = 10 m/s, Omega = 14 RPM, Pitch = 0 deg
V_part1     = 10.0;                  % [m/s]
omega_part1 = 14.0 * pi / 30;        % [rad/s]
theta_part1 = 0.0;                   % [deg]

[P_part1, T_part1, Cp_part1, Ct_part1] = solveBEMs(theta_part1, V_part1, omega_part1, bladeData, airfoils);

fprintf('=== Deliverable 1 Results (V = 10 m/s, RPM = 14, Pitch = 0 deg) ===\n');
fprintf('Coefficient of Power (Cp):  %.4f\n', Cp_part1);
fprintf('Coefficient of Thrust (Ct): %.4f\n\n', Ct_part1);

%  DELIVERABLE 4: Blade Pitch Control at V = 19 m/s, Omega = 15.5 RPM
V_wind   = 19.0;               % [m/s] Assigned wind speed
rpm      = 15.5;               % [RPM] Maximum rotational speed
omega    = rpm * pi / 30;      % [rad/s]
P_target = 2.5e6;              % [W] Rated power limit (2.5 MW)

% Run the pitch sweep function (0 to 25 deg in 0.5 deg increments)
theta_range = -7.5:0.5:25;
[P_sweep, T_sweep, Cp_sweep, Ct_sweep] = sweepPitch(theta_range, V_wind, omega, bladeData, airfoils);

% Pinpoint the exact pitch angle where Power = 2.5 MW using fzero
objFun      = @(theta_p) solveBEMs(theta_p, V_wind, omega, bladeData, airfoils) - P_target;
theta_exact = fzero(objFun, [0, 25]);

% Evaluate final operating state and span wise loads
[P_final, T_final, Cp_final, Ct_final, radialLoads] = solveBEMs(theta_exact, V_wind, omega, bladeData, airfoils);

fprintf('=== Deliverable 4 Results (V = %.1f m/s, RPM = %.1f) ===\n', V_wind, rpm);
fprintf('Exact Pitch Angle (theta_p): %.3f deg\n', theta_exact);
fprintf('Aerodynamic Power:           %.4f MW\n', P_final / 1e6);
fprintf('Rotor Thrust:                %.2f kN\n', T_final / 1e3);
fprintf('Power Coefficient (Cp):      %.4f\n', Cp_final);
fprintf('Thrust Coefficient (Ct):     %.4f\n', Ct_final);


%% ========================================================================
%  LOCAL FUNCTIONS
%  ========================================================================

function [P_vec, T_vec, Cp_vec, Ct_vec] = sweepPitch(theta_vec, V_wind, omega, bladeData, airfoils)
% SWEEPPITCH Sweeps across a vector of blade pitch angles and returns rotor performance
    nSteps = length(theta_vec);
    P_vec  = zeros(1, nSteps);
    T_vec  = zeros(1, nSteps);
    Cp_vec = zeros(1, nSteps);
    Ct_vec = zeros(1, nSteps);

    for k = 1:nSteps
        [P_vec(k), T_vec(k), Cp_vec(k), Ct_vec(k)] = solveBEMs(theta_vec(k), V_wind, omega, bladeData, airfoils);
    end
end


function [P_total, T_total, Cp, Ct, radialLoads] = solveBEMs(theta_p, V_wind, omega, bladeData, airfoils)
% Computes rotor power, thrust, Cp, Ct, and span wise loads
    % Geometry and fluid properties from reference script
    r     = bladeData.DistanceFromCenterOfRotation / 1000;  % [m] (43x1)
    c     = bladeData.ChordLength / 1000;                   % [m] (43x1)
    beta  = bladeData.BladeTwist;                           % [deg] (43x1)
    R     = max(r);                                         % [m] 48 m
    B     = 3;                                              % Number of blades
    A_rot = (R^2) * pi;                                     % Swept area [m^2]
    rho   = 1.1;                                            % Air density [kg/m^3]
    mu    = 1.8e-5;                                         % Dynamic viscosity [Pa*s]
    a     = 1 / 3;                                          % Axial induction factor

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
    Re        = rho .* V_wind .* r ./ mu;

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
            CL_val  = interp1(afTable.AoA, afTable.CL, alpha_deg(i), 'linear');
            CD_val  = interp1(afTable.AoA, afTable.CD, alpha_deg(i), 'linear');
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
% Computes normal (Fz) and tangential (Ftheta) forces per unit span
    FL = CL * 0.5 * rho * (V_rel^2) * c;
    FD = CD * 0.5 * rho * (V_rel^2) * c;

    F_z     = (FL * cosd(phi)) + (FD * sind(phi));
    F_theta = (FL * sind(phi)) - (FD * cosd(phi));
end


function afTable = getAirfoilTable(afName, airfoils)
% Maps BladeProfile airfoil name to its corresponding table
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