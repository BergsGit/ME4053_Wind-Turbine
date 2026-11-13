%ME4053 Wind Turbine Project

%% Clipper Liberty C96 Wind Turbine - BEM Aerodynamic Analysis
% Solves Deliverable 1a (Cp) and 1b (Ct) at U_inf = 10 m/s, 14 RPM, 0 deg pitch.
% Reads blade geometry from 'BladeProfile.csv' and airfoil polars from CSV file(s).

clear; clc; close all;

%% 1. Operating Conditions & Turbine Specifications
U_inf     = 10.0;              % Freestream wind speed [m/s]
rpm       = 14.0;              % Rotor speed [RPM]
pitch_deg = 0.0;               % Blade pitch angle [deg]
B         = 3;                 % Number of blades [-]
R         = 48.0;              % Rotor radius [m]
rho       = 1.225;             % Air density [kg/m^3]

omega     = rpm * 2 * pi / 60; % Angular velocity [rad/s]
TSR       = omega * R / U_inf; % Tip-speed ratio [-]

% BEM & Airfoil CSV options
use_prandtl_tip  = true;            % Prandtl tip-loss factor
use_prandtl_hub  = true;            % Prandtl hub-loss factor
use_glauert      = true;            % Glauert high-induction correction (a > 0.2)
circle_Cd        = 0.50;            % Drag coefficient for cylindrical root stations
default_csv      = 'DU96-W-180.csv';% Fallback / single-airfoil CSV file
force_single_csv = false;           % Set true to force DU96-W-180.csv for all airfoil stations

%% 2. Load Blade Geometry from BladeProfile.csv
opts = detectImportOptions('BladeProfile.csv');
opts = setvartype(opts, 'Airfoil', 'string');
bladeTbl = readtable('BladeProfile.csv', opts);

r         = bladeTbl.DistanceFromCenterOfRotation / 1000.0; % Convert mm -> m
c         = bladeTbl.ChordLength / 1000.0;                  % Convert mm -> m
twist_deg = bladeTbl.BladeTwist;                            % Local twist [deg]
airfoils  = erase(bladeTbl.Airfoil, " ");                   % Normalize names (e.g., "DU96-W-180")

R_hub     = r(1);                                           % Hub radius = 1.3 m
N         = numel(r);

%% 3. Load Airfoil Polar Data from CSV File(s)
polars = load_airfoil_csvs(airfoils, default_csv, force_single_csv);

%% 4. Station-by-Station BEM Solver
a         = 0.25 * ones(N, 1);  % Axial induction factor [-]
ap        = 0.01 * ones(N, 1);  % Tangential induction factor [-]
phi_deg   = zeros(N, 1);        % Inflow angle [deg]
alpha_deg = zeros(N, 1);        % Angle of attack [deg]
Cl_arr    = zeros(N, 1);        % Lift coefficient [-]
Cd_arr    = zeros(N, 1);        % Drag coefficient [-]
F_arr     = ones(N, 1);         % Combined Prandtl loss factor [-]
dT_dr     = zeros(N, 1);        % Thrust per unit span [N/m]
dQ_dr     = zeros(N, 1);        % Torque per unit span [N*m/m]

max_iter  = 500;
tol       = 1e-7;
relax     = 0.5;                % Under-relaxation factor for numerical stability

for i = 1:N
    % Enforce zero aerodynamic load at exact hub and tip boundaries where F = 0
    if (use_prandtl_hub && r(i) <= R_hub) || (use_prandtl_tip && r(i) >= R)
        a(i)     = 1.0;
        ap(i)    = 0.0;
        F_arr(i) = 0.0;
        dT_dr(i) = 0.0;
        dQ_dr(i) = 0.0;
        continue;
    end
    
    beta_deg = twist_deg(i) + pitch_deg;      % Local pitch + twist [deg]
    sigma    = B * c(i) / (2 * pi * r(i));    % Local solidity [-]
    af_name  = airfoils(i);
    
    for iter = 1:max_iter
        % Inflow angle phi [rad]
        phi = atan2(U_inf * (1 - a(i)), omega * r(i) * (1 + ap(i)));
        phi = max(phi, 1e-6);
        
        alpha = rad2deg(phi) - beta_deg;
        
        % Interpolate Cl and Cd from loaded CSV polar table
        if af_name == "circle"
            Cl = 0.0;
            Cd = circle_Cd;
        else
            tbl = polars.(matlab.lang.makeValidName(af_name));
            Cl  = interp1(tbl.AoA, tbl.CL, alpha, 'linear', 'extrap');
            Cd  = interp1(tbl.AoA, tbl.CD, alpha, 'linear', 'extrap');
        end
        
        sin_p = sin(phi);
        cos_p = cos(phi);
        
        % Normal and tangential force coefficients
        Cn = Cl * cos_p + Cd * sin_p;
        Ct = Cl * sin_p - Cd * cos_p;
        
        % Prandtl tip and hub loss factors
        F_tip = 1.0;
        if use_prandtl_tip
            f_tip = (B / 2) * (R - r(i)) / (r(i) * sin_p);
            F_tip = (2 / pi) * acos(min(1.0, exp(-min(50.0, f_tip))));
        end
        
        F_hub = 1.0;
        if use_prandtl_hub
            f_hub = (B / 2) * (r(i) - R_hub) / (r(i) * sin_p);
            F_hub = (2 / pi) * acos(min(1.0, exp(-min(50.0, f_hub))));
        end
        
        F = max(1e-6, F_tip * F_hub);
        
        % Axial induction update (with Hansen/Glauert correction for a > 0.2)
        K     = 4 * F * (sin_p^2) / (sigma * max(1e-6, Cn));
        a_new = 1 / (K + 1);
        
        if use_glauert && a_new > 0.2
            ac    = 0.2;
            term  = (K * (1 - 2*ac) + 2)^2 + 4 * (K * (ac^2) - 1);
            a_new = 0.5 * (2 + K * (1 - 2*ac) - sqrt(max(0.0, term)));
        end
        
        % Tangential induction update
        if abs(Ct) > 1e-8
            Kp     = 4 * F * sin_p * cos_p / (sigma * Ct);
            ap_new = 1 / (Kp - 1);
        else
            ap_new = 0.0;
        end
        
        % Check convergence
        if abs(a_new - a(i)) < tol && abs(ap_new - ap(i)) < tol
            a(i)  = a_new;
            ap(i) = ap_new;
            break;
        end
        
        a(i)  = (1 - relax) * a(i)  + relax * a_new;
        ap(i) = (1 - relax) * ap(i) + relax * ap_new;
    end
    
    % Store converged station parameters and elemental loads
    V_rel_sq     = (U_inf * (1 - a(i)))^2 + (omega * r(i) * (1 + ap(i)))^2;
    dT_dr(i)     = 0.5 * rho * B * c(i) * V_rel_sq * Cn;
    dQ_dr(i)     = 0.5 * rho * B * c(i) * V_rel_sq * Ct * r(i);
    phi_deg(i)   = rad2deg(phi);
    alpha_deg(i) = alpha;
    Cl_arr(i)    = Cl;
    Cd_arr(i)    = Cd;
    F_arr(i)     = F;
end

%% 5. Integrate Spanwise Loads & Compute Cp and Ct
T_total = trapz(r, dT_dr);              % Total rotor thrust [N]
Q_total = trapz(r, dQ_dr);              % Total rotor torque [N*m]
P_total = Q_total * omega;              % Mechanical rotor power [W]

A_swept = pi * R^2;                     % Swept rotor area [m^2]
Cp      = P_total / (0.5 * rho * A_swept * U_inf^3);
Ct      = T_total / (0.5 * rho * A_swept * U_inf^2);

fprintf('======================================================\n');
fprintf('  Clipper Liberty C96 BEM Results (U = %.1f m/s, %.1f RPM)\n', U_inf, rpm);
fprintf('======================================================\n');
fprintf('  Tip-Speed Ratio (TSR) : %.4f\n', TSR);
fprintf('  Rotor Power (P)       : %.3f MW\n', P_total / 1e6);
fprintf('  Rotor Thrust (T)      : %.2f kN\n', T_total / 1e3);
fprintf('  1a. Power Coeff (Cp)  : %.4f\n', Cp);
fprintf('  1b. Thrust Coeff (Ct) : %.4f\n', Ct);
fprintf('======================================================\n');

%% 6. Diagnostic Spanwise Plots
figure('Name', 'C96 BEM Spanwise Distributions', 'Color', 'w');

subplot(2, 2, 1);
plot(r / R, a, 'b-o', 'LineWidth', 1.4, 'MarkerSize', 4); hold on;
plot(r / R, ap, 'r-s', 'LineWidth', 1.4, 'MarkerSize', 4);
grid on; xlabel('r / R [-]'); ylabel('Induction Factor [-]');
legend('Axial (a)', 'Tangential (a'')', 'Location', 'best');
title('Induction Factors');

subplot(2, 2, 2);
plot(r(7:end) / R, alpha_deg(7:end), 'k-o', 'LineWidth', 1.4, 'MarkerSize', 4);
grid on; xlabel('r / R [-]'); ylabel('\alpha [deg]');
title('Angle of Attack (Airfoil Stations)');

subplot(2, 2, 3);
plot(r / R, dT_dr / 1e3, 'b-', 'LineWidth', 1.6);
grid on; xlabel('r / R [-]'); ylabel('dT/dr [kN/m]');
title('Spanwise Thrust Distribution');

subplot(2, 2, 4);
plot(r / R, dQ_dr / 1e3, 'r-', 'LineWidth', 1.6);
grid on; xlabel('r / R [-]'); ylabel('dQ/dr [kN\cdotm/m]');
title('Spanwise Torque Distribution');

%% Helper Function: Load Airfoil Polar CSV Files
function polars = load_airfoil_csvs(airfoils, default_csv, force_single_csv)
    unique_afs = unique(airfoils(airfoils ~= "circle"));
    default_tbl = readtable(default_csv);
    
    for k = 1:numel(unique_afs)
        af = unique_afs(k);
        field = matlab.lang.makeValidName(af);
        csv_file = af + ".csv";
        
        if ~force_single_csv && isfile(csv_file)
            polars.(field) = readtable(csv_file);
            fprintf('Loaded airfoil polar from: %s\n', csv_file);
        else
            polars.(field) = default_tbl;
            fprintf('Using %s for airfoil: %s\n', default_csv, af);
        end
    end
end