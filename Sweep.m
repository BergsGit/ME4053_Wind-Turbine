clear; clc; close all;
%ME4053 Wind Turbine Project
%Group 26
%
%Authors:
%...
%%
%_____Table of Contents____
% 1. Deliverables
%       1. Deliverable 1: Baseline Cp and Ct (V = 10 m/s, 14 RPM, 0 deg)
%       2. Deliverable 2: Pitch Angle Sweep at Fixed TSR (V = 8.5 m/s)
%       3. Deliverable 3: Pitch Angle & TSR Optimization (V = 8.3 m/s)
%       4. Deliverable 4: Rated Power Pitch Control (V = 19 m/s, 15.5 RPM)
%       5. Deliverable 5: Tower Structural Analysis (Deflection, Stress, Vibration, Buckling, Fatigue)
% 2. Functions
%       Import Data and Parameters (LoadAirfoils, getAirfoilTable)
%       Blade Analysis (solveBEMs, calcElementForces, sweepPitch)
%       Tower Analysis (analyzeTower)

%% 1. Deliverables

% Load blade geometry, tower geometry, and airfoil lookup tables
bladeData = readtable('BladeProfile.csv');
towerData = readtable('towerSpecs.csv');
airfoils  = LoadAirfoils();

% DELIVERABLE 1: V = 10 m/s, Omega = 14 RPM, Pitch = 0 deg
    V_part1      = 10.0;                  % [m/s]
    omega_part1  = 14.0 * pi / 30;        % [rad/s]
    Pitch_Angle1 = 0.0;                   % [deg]

    [P_part1, T_part1, Cp_part1, Ct_part1] = solveBEMs(Pitch_Angle1, ...
        V_part1, omega_part1, bladeData, airfoils);

    fprintf('=== Deliverable 1 Results ===\n');
    fprintf('Coefficient of Power (Cp):  %.4f\n', Cp_part1);
    fprintf('Coefficient of Thrust (Ct): %.4f\n\n', Ct_part1);

% DELIVERABLE 2: Blade Pitch Control at V = 8.5 m/s, TipSpeedRatio = 7.10
    U_wind2  = 8.5;                 % [m/s] Assigned wind speed
    TSR      = 7.10;                % [] Tip Speed Ratio
    omega2   = TSR * U_wind2 / 48;  % [rad/s] Rotational speed from TSR

    % Run the pitch sweep function (-10 to 20 deg in 0.05 deg increments)
    theta_range = -10:0.05:20;
    [~, ~, Cp_sweep, ~] = sweepPitch(theta_range, ...
        U_wind2, omega2, bladeData, airfoils);

    % Find maximum Cp value and its corresponding pitch angle
    [max_Cp, index] = max(Cp_sweep);
    maxCp_angle = theta_range(index);

    fprintf('=== Deliverable 2 Results ===\n');
    fprintf('Max Coefficient of Power:  %.4f\n', max_Cp);
    fprintf('Pitch Angle at Max C_p:    %.4f deg\n\n', maxCp_angle);

    % Plot Figure 1
    figure(1)
    plot(theta_range, Cp_sweep, 'LineWidth', 1.5);
    grid on;
    xlabel('Pitch Angle (degrees)')
    ylabel('C_p')
    set(gca, 'FontName', 'Times New Roman', 'FontSize', 12);

% DELIVERABLE 3: Blade Pitch Control at V = 8.3 m/s
    U_wind3  = 8.3;                  % [m/s] Assigned wind speed
    TSR      = 5.85:0.05:9.35;       % [] Tip Speed Ratio 
    omega3   = TSR .* U_wind3 / 48;  % [rad/s] Rotational speed from TSR

    % Run the pitch sweep function (-10 to 20 deg in 0.1 deg increments)
    theta_range = -10:0.1:20;
    
    % Array to store Cp data in
    Cp_data = zeros(numel(omega3), numel(theta_range));
    
    % Sweep for each TSR value
    for i = 1:numel(omega3)
        [~, ~, Cp_sweep, ~] = sweepPitch(theta_range, ...
            U_wind3, omega3(i), bladeData, airfoils);
        %Store data
        Cp_data(i,:) = Cp_sweep;
    end

    % Find maximum Cp value and its corresponding pitch angle and TSR
    max_Cp = max(Cp_data(:));
    [row, col] = find(Cp_data == max_Cp);
    
    optimal_TSR        = TSR(row);
    optimal_PitchAngle = theta_range(col);

    fprintf('=== Deliverable 3 Results ===\n');
    fprintf('Max Coefficient of Power:     %.4f\n', max_Cp);
    fprintf('Pitch Angle at Max C_p:       %.4f deg\n', optimal_PitchAngle);
    fprintf('Tip Speed Ratio at Max C_p:   %.4f\n\n', optimal_TSR);
    
    % Plot Figure 2
    figure(2)
    contourf(theta_range, TSR, Cp_data, 20, 'LineColor', 'none');
    xlabel('Pitch Angle (degrees)')
    ylabel('Tip Speed Ratio')
    c = colorbar;
    c.Label.String = 'C_p';
    set(gca, 'FontName', 'Times New Roman', 'FontSize', 12);

% DELIVERABLE 4: Blade Pitch Control at V = 19 m/s, Omega = 15.5 RPM
    V_wind   = 19.0;               % [m/s] Assigned wind speed
    rpm      = 15.5;               % [RPM] Maximum rotational speed
    omega    = rpm * pi / 30;      % [rad/s]
    P_target = 2.5e6;              % [W] Rated power limit (2.5 MW)

    % Run the pitch sweep function (-7.5 to 25 deg in 0.5 deg increments)
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
    Q_final = P_final / omega;     % [N*m] Rotor torque at 2.5 MW operating point

    fprintf('=== Deliverable 4 Results ===\n');
    fprintf('Exact Pitch Angle (theta_p): %.3f deg\n', theta_exact);
    fprintf('Aerodynamic Power:           %.4f MW\n', P_final / 1e6);
    fprintf('Rotor Thrust:                %.2f kN\n', T_final / 1e3);
    fprintf('Rotor Torque:                %.2f kN*m\n', Q_final / 1e3);
    fprintf('Power Coefficient (Cp):      %.4f\n', Cp_final);
    fprintf('Thrust Coefficient (Ct):     %.4f\n\n', Ct_final);

% DELIVERABLE 5: Tower Structural Analysis at Deliverable 4 Loading
    towerResults = analyzeTower(towerData, bladeData, V_wind, rpm, T_final, Q_final);


%% 2. Functions

%_____Import Functions_____
function airfoils = LoadAirfoils()
% Function: LoadAirfoils
% Input:    none
% Output:   airfoils (struct containing airfoil tables)
% Process:  For each unique airfoil, the corresponding csv is read in and
%           stored in the data struct under its code.
    airfoils.DU97 = readtable('DU97-W-300.csv');
    airfoils.DU91 = readtable('DU91-W2-250.csv');
    airfoils.DU93 = readtable('DU93-W-210.csv');
    airfoils.DU96 = readtable('DU96-W-180.csv');
end

function afTable = getAirfoilTable(afName, airfoils)
% Function: getAirfoilTable
% Input:    afName (airfoil name), airfoils (struct containing airfoil tables)
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


%_____Wind Turbine Blade Analysis_____
function [P_total, T_total, Cp, Ct, radialLoads] = solveBEMs(theta_p, V_wind, omega, bladeData, airfoils)
% Computes rotor power, thrust, Cp, Ct, and span wise loads
% Function: solveBEMs
% Input:    theta_p (pitch angle), V_wind (wind speed), omega (rotational 
%           velocity), bladeData (table containing blade parameters), 
%           airfoils (struct containing airfoil tables)
% Output:   P_total (output power value), T_total (total thrust value),
%           Cp (coefficient of power), Ct (coefficient of thrust),
%           radialLoads (struct containing radial load at each r value)
% Process:  First, tower data, parameters, and air properties are defined.
%           Next, the Tip speed ratio is found and used to define the local
%           tip speed ratio which in turn is used to find the axial'
%           induction factor for each radial position on the blade. The 
%           relative velocity is then found through the vector sum of the 
%           converted radial velocity and wind speed. This velocity can then
%           be used to find the reynolds number and AoA.
%
%           The arrays for the output parameters are then initialized to
%           the size of the 'r' array. And their values found through
%           analysis and linear interpolation of the imported tables.
%           trapz() is then used to find the sum of these arrays and the
%           final output parameters are calculated.
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

    % Package radial loads
    radialLoads.r       = r;
    radialLoads.F_theta = F_theta';
    radialLoads.F_z     = F_z';
end

function [F_theta, F_z, CL, CD, FL, FD] = calcElementForces(CL, CD, rho, V_rel, c, phi)
% Function: calcElementForces
% Input:    CL, CD, rho, V_rel, c, phi
% Output:   F_theta, F_z, CL, CD, FL, FD
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
% Input:    theta_vec, V_wind, omega, bladeData, airfoils
% Output:   P_vec, T_vec, Cp_vec, Ct_vec
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


function results = analyzeTower(towerData, bladeData, V_wind, rpm, T_final, Q_final)
% ANALYZETOWER Simplified structural analysis for Deliverable 5 (Deflection, Stress,
% Vibration, Buckling, and Fatigue) using direct beam integration & Rayleigh's method.

    % 1. Material & Turbine Properties (ASTM A572 Grade 50 Steel)
    rho_air   = 1.1;          % [kg/m^3]
    mu_air    = 1.8e-5;       % [N*s/m^2]
    g         = 9.81;         % [m/s^2]
    E_steel   = 200e9;        % [Pa]
    nu_steel  = 0.30;         % [-]
    rho_steel = 7850;         % [kg/m^3]
    Su        = 450e6;        % [Pa] Ultimate tensile strength
    Sy        = 345e6;        % [Pa] Yield strength

    m_rotor   = 55.7e3;       % [kg]
    m_blade   = 10.7e3;       % [kg]
    m_nacelle = 94.3e3;       % [kg]
    m_top     = m_rotor + 3*m_blade + m_nacelle; % [kg] Total RNA top mass (182,100 kg)
    H_hub     = 80.4;         % [m] Hub height
    omega     = rpm * pi / 30;% [rad/s]

    % 2. Tower Geometry (Convert mm to m)
    z_tower = towerData{:, 1} / 1000;        % [m]
    OD      = towerData{:, 2} / 1000;        % [m]
    t_wall  = towerData{:, 3} / 1000;        % [m]
    ID      = OD - 2 .* t_wall;              % [m]
    c_out   = OD / 2;                        % [m]
    A_sec   = (pi / 4)  .* (OD.^2 - ID.^2);  % [m^2]
    I_sec   = (pi / 64) .* (OD.^4 - ID.^4);  % [m^4]
    J_sec   = 2 .* I_sec;                    % [m^4]

    % Assumption 4: ABL 1/7 power-law wind profile & distributed drag on tower
    V_tower  = V_wind .* (z_tower ./ H_hub) .^ (1/7);
    Re_tower = rho_air .* V_tower .* OD ./ mu_air;
    Cd_tower = zeros(size(z_tower));
    Cd_tower(Re_tower > 0) = arrayfun(@cylinderCD, Re_tower(Re_tower > 0));
    f_wind   = 0.5 * rho_air * (V_tower.^2) .* Cd_tower .* OD; % [N/m]

    % 3. Internal Loads & Static Deflection (Direct Beam Integration)
    nNodes  = length(z_tower);
    V_shear = zeros(nNodes, 1);
    M_bend  = zeros(nNodes, 1);
    W_above = zeros(nNodes, 1);

    for i = 1:(nNodes - 1)
        V_shear(i) = T_final + trapz(z_tower(i:end), f_wind(i:end));
        M_bend(i)  = T_final * (H_hub - z_tower(i)) + ...
                     trapz(z_tower(i:end), f_wind(i:end) .* (z_tower(i:end) - z_tower(i)));
        W_above(i) = m_top * g + trapz(z_tower(i:end), rho_steel .* A_sec(i:end) .* g);
    end
    V_shear(end) = T_final;
    M_bend(end)  = T_final * (H_hub - z_tower(end));
    W_above(end) = m_top * g;

    % Double integration of curvature M / (E*I) with fixed base (slope(0)=0, defl(0)=0)
    slope      = cumtrapz(z_tower, M_bend ./ (E_steel .* I_sec));
    deflection = cumtrapz(z_tower, slope); % [m]

    % 4. Combined Stress Analysis
    sigma_axial = W_above ./ A_sec;
    sigma_bend  = (M_bend .* c_out) ./ I_sec;
    sigma_comp  = sigma_axial + sigma_bend;
    tau_torsion = (Q_final .* c_out) ./ J_sec;
    tau_shear   = 2 .* V_shear ./ A_sec;
    sigma_vm    = sqrt(sigma_comp.^2 + 3 .* (tau_torsion + tau_shear).^2);

    [max_bend, idx_bend] = max(sigma_bend);
    [max_vm, idx_vm]     = max(sigma_vm);
    SF_yield             = Sy / max_vm;

    % 5. Vibration & Dynamic Deflection (Rayleigh Stiffness + 3P Square-Wave Pulse)
    % Unit top load deflection -> equivalent fundamental stiffness k_eq
    M_unit     = 1.0 * (z_tower(end) - z_tower);
    defl_unit  = cumtrapz(z_tower, cumtrapz(z_tower, M_unit ./ (E_steel .* I_sec)));
    mode_shape = defl_unit / defl_unit(end); % Normalized 1st bending mode shape
    k_eq       = 1.0 / defl_unit(end);

    m_tower_tot = trapz(z_tower, rho_steel .* A_sec);
    m_eff       = m_top + 0.23 * m_tower_tot; % Effective modal mass
    omega_1     = sqrt(k_eq / m_eff);
    f_1         = omega_1 / (2 * pi);         % 1st natural frequency [Hz]
    f_1P        = rpm / 60;
    f_3P        = 3 * f_1P;

    % Assumptions 6, 7, & 8: Binary shadow square-wave pulse with constant tower OD
    D_avg      = mean(OD);
    r_eff      = 0.75 * max(bladeData.DistanceFromCenterOfRotation / 1000);
    pulseWidth = D_avg / (omega * r_eff);
    delta_T    = (1/3) * T_final;
    [mag_harm, phase_harm] = pulseFFT(f_3P, pulseWidth, delta_T, 5);

    % SDOF steady-state harmonic response along mode shape (zeta = 1% damping)
    zeta    = 0.01;
    t_cycle = linspace(0, 2 / f_3P, 200);
    u_top_t = zeros(size(t_cycle));
    for h = 2:5
        Omega_h = 2 * pi * (h - 1) * f_3P;
        H_freq  = 1 / ((k_eq - m_eff * Omega_h^2) + 1i * (2 * zeta * omega_1 * m_eff * Omega_h));
        u_top_t = u_top_t + real(mag_harm(h) * H_freq .* exp(1i * (Omega_h * t_cycle + phase_harm(h))));
    end
    dyn_env_max = max(abs(u_top_t)) .* mode_shape;
    dyn_amp_1st = abs(mag_harm(2) / ((k_eq - m_eff * (2*pi*f_3P)^2) + 1i * (2*zeta*omega_1*m_eff*2*pi*f_3P))) .* mode_shape;

    % 6. Buckling Analysis (Rayleigh-Ritz Euler Column + Local Shell)
    L_tot    = z_tower(end);
    dy       = (pi / (2 * L_tot)) .* sin(pi .* z_tower ./ (2 * L_tot));
    d2y      = (pi / (2 * L_tot))^2 .* cos(pi .* z_tower ./ (2 * L_tot));
    SF_euler = trapz(z_tower, E_steel .* I_sec .* d2y.^2) / trapz(z_tower, W_above .* dy.^2);

    R_mid    = (OD - t_wall) / 2;
    sigma_cl = (E_steel / sqrt(3 * (1 - nu_steel^2))) .* (t_wall ./ R_mid);
    alpha_kd = 0.83 ./ sqrt(1 + 0.01 .* (R_mid ./ t_wall));
    sigma_cr = alpha_kd .* sigma_cl;
    [SF_shell, idx_buck] = min(sigma_cr ./ sigma_comp);

    % 7. Fatigue Analysis (Modified Goodman + S-N Curve)
    Sn_prime   = 0.5 * Su;
    Sn         = Sn_prime * 1.0 * 0.7 * 0.72 * 1.0 * 0.814; % CL * CG * CS * CT * CR
    Kf         = 1 + 1.0 * (1.5 - 1);                       % q = 1.0, Kt = 1.5
    sigma_m    = Kf * sigma_comp(idx_vm);
    sigma_a    = Kf * (mag_harm(2) * (H_hub - z_tower(idx_vm)) * c_out(idx_vm) / I_sec(idx_vm));
    SF_fatigue = 1 / ((sigma_a / Sn) + (sigma_m / Su));

    Sf_1e3 = 0.9 * Su;
    b_exp  = -log10(Sf_1e3 / Sn) / 3;
    A_coef = Sf_1e3 / (10^3)^b_exp;

    % Print Summary
    fprintf('=== Deliverable 5: Simplified Tower Structural Results ===\n');
    fprintf('Max Static Deflection:        %.1f mm\n', deflection(end)*1000);
    fprintf('Max Bending Stress:           %.2f MPa (Max VM = %.2f MPa, Yield SF = %.2f)\n', ...
            max_bend/1e6, max_vm/1e6, SF_yield);
    fprintf('Natural Frequency (f1):   %.4f Hz (1P = %.4f Hz, 3P = %.4f Hz)\n', f_1, f_1P, f_3P);
    fprintf('   Peak Dynamic Deflection (3P): +/- %.2f mm\n', dyn_env_max(end)*1000);
    fprintf('Buckling SF:     %.2f | Local Shell Buckling SF: %.2f\n', SF_euler, SF_shell);
    fprintf('Endurance Limit Sn: %.2f MPa | Goodman Fatigue SF: %.2f\n\n', Sn/1e6, SF_fatigue);

    % Plots (Figures 4-7)
    figure(4);
    plot(deflection * 1000, z_tower, 'b-', 'LineWidth', 2.0); grid on;
    xlabel('Static Lateral Deflection (mm)'); ylabel('Tower Height z (m)');

    figure(5);
    plot(dyn_env_max * 1000, z_tower, 'r-', 'LineWidth', 2.0); hold on;
    plot(-dyn_env_max * 1000, z_tower, 'r-', 'LineWidth', 2.0, 'HandleVisibility', 'off');
    plot(dyn_amp_1st * 1000, z_tower, 'k--', 'LineWidth', 1.5);
    plot(-dyn_amp_1st * 1000, z_tower, 'k--', 'LineWidth', 1.5, 'HandleVisibility', 'off');
    grid on;
    xlabel('Dynamic Deflection Amplitude \Delta x (mm)'); ylabel('Tower Height z (m)');
    legend('Multi-Harmonic Peak Envelope', 'Fundamental 3P Harmonic', 'Location', 'southeast');

    figure(6);
    p1=plot(sigma_bend / 1e6, z_tower, 'b-', 'LineWidth', 2.0); hold on;
    p2=xline(Sy / 1e6, 'r--', sprintf('Yield Strength S_y = %.0f MPa', Sy/1e6), 'LineWidth', 1.5);
    grid on;
    xlabel('Stress (MPa)'); ylabel('Tower Height z (m)');
  legend([p1, p2], {'Bending Stress (\sigma_b)', 'Yield Strength (S_y)'}, 'Location', 'southeast');
    xlim([0, 380]);

    N_curve = [logspace(0, 3, 50), logspace(3, 6, 100), logspace(6, 9, 50)];
    S_curve = [linspace(Su, Sf_1e3, 50), A_coef .* (logspace(3, 6, 100) .^ b_exp), Sn .* ones(1, 50)];
    figure(7);
    loglog(N_curve, S_curve / 1e6, 'b-', 'LineWidth', 2.2); hold on;
    loglog([1e3, 1e6], [Sf_1e3, Sn] / 1e6, 'bo', 'MarkerFaceColor', 'b', 'MarkerSize', 7);
    yline(Sn / 1e6, 'k--', sprintf('Endurance Limit S_n = %.1f MPa', Sn/1e6), 'LineWidth', 1.2);
    loglog(1e7, sigma_a / 1e6, 'rp', 'MarkerFaceColor', 'r', 'MarkerSize', 11);
    grid on;
    xlabel('Cycles to Failure, N (Log Scale)'); ylabel('Fatigue Strength / Stress Amplitude S (MPa)');
    legend('Corrected S-N Curve', '10^3 & 10^6 Anchor Points', 'Endurance Limit (S_n)', 'Tower Design Point', 'Location', 'southwest');
    xlim([1e2, 1e9]); ylim([1, 1000]);

    results = struct('deflection', deflection, 'sigma_bend', sigma_bend, 'sigma_vm', sigma_vm, ...
                     'SF_yield', SF_yield, 'f_1', f_1, 'SF_euler', SF_euler, 'SF_shell', SF_shell, 'SF_fatigue', SF_fatigue);
end