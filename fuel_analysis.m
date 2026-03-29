function fuel_report = fuel_analysis(maneuver_log, satellites, params)
% Summarize how much fuel each satellite burned and estimate how long they'll last.
%
% Uses the Tsiolkovsky rocket equation to work out remaining delta-V capacity,
% then divides by a typical LEO stationkeeping budget to get years of life left.

g0     = 9.80665;
n_sats = length(satellites);

fuel_report.satellites  = satellites;
fuel_report.total_dv    = 0;
fuel_report.total_fuel  = 0;

fprintf('\n--- FUEL & MISSION LIFETIME REPORT ---\n');
fprintf('%-10s %-12s %-12s %-12s %-15s %-15s\n', ...
    'Satellite','Maneuvers','dV used(m/s)','Fuel used(kg)','Fuel left(kg)','Est. Life(yr)');
fprintf('%s\n', repmat('-',1,80));

% Rough stationkeeping budget at 550 km altitude
yearly_dv = 5.0;   % m/s/year

for s = 1:n_sats
    sat = satellites(s);

    if ~isempty(maneuver_log) && isstruct(maneuver_log) && isfield(maneuver_log,'sat_idx')
        rows  = ([maneuver_log.sat_idx] == s);
        dv    = sum([maneuver_log(rows).delta_v_ms]);
        spent = sum([maneuver_log(rows).fuel_used_kg]);
        burns = sum(rows);
    else
        dv = 0; spent = 0; burns = 0;
    end

    init_fuel  = sat.fuel_mass;
    fuel_left  = max(0, init_fuel - spent);

    % Remaining delta-V from rocket equation
    wet_left   = sat.mass + max(fuel_left, 0);
    dv_left    = sat.Isp * g0 * log(wet_left / sat.mass);

    % Years of life based on what's left after maneuvers
    life_yrs   = max(0, (dv_left - dv) / yearly_dv);

    fprintf('%-10s %-12d %-12.3f %-12.3f %-15.3f %-15.2f\n', ...
        sat.name, burns, dv, spent, fuel_left, life_yrs);

    fuel_report.total_dv   = fuel_report.total_dv   + dv;
    fuel_report.total_fuel = fuel_report.total_fuel + spent;

    fuel_report.per_sat(s).name        = sat.name;
    fuel_report.per_sat(s).n_maneuvers = burns;
    fuel_report.per_sat(s).dv_used     = dv;
    fuel_report.per_sat(s).fuel_used   = spent;
    fuel_report.per_sat(s).fuel_left   = fuel_left;
    fuel_report.per_sat(s).est_life_yr = life_yrs;
end

fprintf('%s\n', repmat('-',1,80));
fprintf('TOTAL: dV = %.3f m/s | Fuel consumed = %.3f kg\n\n', ...
    fuel_report.total_dv, fuel_report.total_fuel);
end
