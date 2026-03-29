%% =========================================================================
%  LEO Conjunction Avoidance Simulation
%  Runs a synthetic satellite constellation through a debris field,
%  detects close approaches, fires avoidance burns, and reports fuel usage.
%% =========================================================================

clc; clear; close all;
addpath(genpath(pwd));

fprintf('========================================================\n');
fprintf('   LEO Satellite Conjunction Avoidance Simulation\n');
fprintf('========================================================\n\n');

%% Simulation setup
params.mu       = 3.986004418e14;   % Earth GM [m^3/s^2]
params.Re       = 6371e3;           % Earth radius [m]
params.J2       = 1.08263e-3;       % J2 oblateness coefficient
params.dt       = 10;               % time step [s]
params.T_days   = 5;
params.T_total  = params.T_days * 86400;
params.t_vec    = 0 : params.dt : params.T_total;
params.N_steps  = length(params.t_vec);

% Thresholds for conjunction alerts
params.warn_dist_km   = 50;    % [km] - flag for monitoring
params.danger_dist_km = 10;    % [km] - triggers a burn
params.Pc_threshold   = 1e-4;  % Pc above this also triggers a burn

fprintf('Simulation: %.0f days | dt = %d s | Steps = %d\n\n', ...
    params.T_days, params.dt, params.N_steps);

%% Step 1: Build the constellation
fprintf('[1/6] Generating satellite constellation...\n');
satellites = generate_satellites(params);
fprintf('      %d satellites initialized.\n\n', length(satellites));

%% Step 2: Seed the debris field
fprintf('[2/6] Generating synthetic debris field...\n');
debris = generate_debris(params);
fprintf('      %d debris objects initialized.\n\n', length(debris));

%% Step 3: Propagate everything
fprintf('[3/6] Propagating orbits (Two-Body + J2 perturbation)...\n');
[sat_states, debris_states] = propagate_all_orbits(satellites, debris, params);
fprintf('      Propagation complete.\n\n');

%% Step 4: Find close approaches and dodge them
fprintf('[4/6] Running conjunction detection...\n');
[conjunctions, maneuver_log] = detect_and_avoid(sat_states, debris_states, ...
                                                  satellites, params);
fprintf('      %d conjunction events detected.\n', length(conjunctions));
fprintf('      %d avoidance maneuvers executed.\n\n', length(maneuver_log));

%% Step 5: Fuel and lifetime accounting
fprintf('[5/6] Computing fuel consumption and mission lifetime...\n');
fuel_report = fuel_analysis(maneuver_log, satellites, params);

%% Step 6: Plots
fprintf('[6/6] Generating visualizations...\n');
visualize_simulation(sat_states, debris_states, conjunctions, ...
                     maneuver_log, fuel_report, params);

fprintf('\n========================================================\n');
fprintf('   Simulation Complete.\n');
fprintf('========================================================\n');
