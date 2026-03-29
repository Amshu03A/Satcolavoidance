function [conjunctions, maneuver_log] = detect_and_avoid(sat_states, debris_states, satellites, params)
% Check for close approaches and dodge them if needed.
%
% Steps:
%  1. For each sat-debris pair, find the closest point in time.
%  2. If close enough to warrant a warning, log it.
%  3. If dangerously close or Pc is too high, plan a burn.
%  4. Burn: prograde + slight out-of-plane nudge ~30 min before closest pass.
%  5. Update the satellite trajectory from that point forward.

mu         = params.mu;
warn_km    = params.warn_dist_km * 1e3;
danger_km  = params.danger_dist_km * 1e3;
pc_thresh  = params.Pc_threshold;
dt         = params.dt;
t_vec      = params.t_vec;
num_steps  = params.N_steps;
num_sats   = size(sat_states, 1);
num_debris = size(debris_states, 1);

% Cell arrays avoid MATLAB's struct-assignment quirks
close_calls = {};
burns       = {};
cc_count    = 0;
burn_count  = 0;

% Don't let a sat burn twice in quick succession
cooldown = zeros(num_sats, 1);

for s = 1:num_sats
    sat_name = satellites(s).name;
    isp      = satellites(s).Isp;
    dry_mass = satellites(s).mass;
    fuel     = satellites(s).fuel_mass;

    for d = 1:num_debris
        % Find closest approach time for this pair
        gap   = squeeze(sat_states(s,:,1:3)) - squeeze(debris_states(d,:,1:3));
        range = sqrt(sum(gap.^2, 2));

        [min_range, tca] = min(range);

        if min_range > warn_km
            continue;
        end

        % Log this close approach
        cc_count = cc_count + 1;

        body_r   = satellites(s).hard_body_r;
        debris_r = 0.5;
        rel_spd  = norm(squeeze(sat_states(s,tca,4:6)) - ...
                        squeeze(debris_states(d,tca,4:6)));
        pc = compute_pc(min_range, body_r, debris_r, rel_spd);

        rec.sat_name    = sat_name;
        rec.sat_idx     = s;
        rec.deb_idx     = d;
        rec.tca_time_s  = t_vec(tca);
        rec.min_dist_km = min_range / 1e3;
        rec.risk_level  = tag_risk(min_range, warn_km, danger_km);
        rec.Pc          = pc;
        close_calls{end+1} = rec;

        % Decide whether to burn
        too_close = (min_range < danger_km) || (pc > pc_thresh);

        if too_close && (tca > cooldown(s)) && (fuel > 0.1)
            % Burn 30 min before closest approach
            burn_step = max(1, tca - round(1800/dt));

            v_sat = squeeze(sat_states(s, burn_step, 4:6))';
            r_pos = squeeze(sat_states(s, burn_step, 1:3))';

            r_hat = r_pos / norm(r_pos);
            v_hat = v_sat / norm(v_sat);
            h_hat = cross(r_pos, v_sat);
            h_hat = h_hat / norm(h_hat);

            dv_mag = 0.5 + rand() * 0.5;   % 0.5-1.0 m/s
            dv_vec = dv_mag * (0.8*v_hat + 0.2*h_hat);

            % How much fuel does this cost?
            g0   = 9.80665;
            wet  = dry_mass + fuel;
            dm   = wet * (1 - exp(-norm(dv_vec) / (isp*g0)));
            dm   = min(dm, fuel);
            fuel = fuel - dm;

            % Shift all subsequent states by the burn
            for k = burn_step:num_steps
                sat_states(s, k, 4) = sat_states(s, k, 4) + dv_vec(1);
                sat_states(s, k, 5) = sat_states(s, k, 5) + dv_vec(2);
                sat_states(s, k, 6) = sat_states(s, k, 6) + dv_vec(3);
                if k > burn_step
                    for ax = 1:3
                        sat_states(s,k,ax) = sat_states(s,k,ax) + dv_vec(ax)*dt*(k-burn_step);
                    end
                end
            end

            % Record the burn
            burn_count = burn_count + 1;
            brec.sat_name           = sat_name;
            brec.sat_idx            = s;
            brec.deb_idx            = d;
            brec.man_time_s         = t_vec(burn_step);
            brec.delta_v_ms         = norm(dv_vec);
            brec.fuel_used_kg       = dm;
            brec.fuel_remain_kg     = fuel;
            brec.Pc_before          = pc;
            brec.min_dist_km_before = min_range / 1e3;
            burns{end+1}            = brec;

            satellites(s).fuel_mass = fuel;
            satellites(s).fuel_used = satellites(s).fuel_used + dm;
            satellites(s).maneuvers = satellites(s).maneuvers + 1;

            % 2-hour cooldown before this sat can burn again
            cooldown(s) = tca + round(7200/dt);

            fprintf('      MANEUVER: %s avoids DEB-%03d | dV=%.3f m/s | dm=%.3f kg | Fuel left=%.2f kg\n', ...
                sat_name, d, norm(dv_vec), dm, fuel);
        end
    end
end

if isempty(close_calls)
    conjunctions = struct('sat_name',{},'sat_idx',{},'deb_idx',{},...
        'tca_time_s',{},'min_dist_km',{},'risk_level',{},'Pc',{});
else
    conjunctions = [close_calls{:}];
end

if isempty(burns)
    maneuver_log = struct('sat_name',{},'sat_idx',{},'deb_idx',{},...
        'man_time_s',{},'delta_v_ms',{},'fuel_used_kg',{},...
        'fuel_remain_kg',{},'Pc_before',{},'min_dist_km_before',{});
else
    maneuver_log = [burns{:}];
end
end

% -------------------------------------------------------------------------
function level = tag_risk(dist_m, warn_m, danger_m)
if dist_m < danger_m
    level = 'CRITICAL';
elseif dist_m < warn_m
    level = 'WARNING';
else
    level = 'NOMINAL';
end
end

% -------------------------------------------------------------------------
function pc = compute_pc(miss_dist, r_sat, r_deb, rel_vel)
% Chan simplified Pc estimate.
% Assumes 100 m 1-sigma position uncertainty on each axis.
sigma = 100;
area  = pi * (r_sat + r_deb)^2;
pc    = (area / (2*pi*sigma^2)) * exp(-miss_dist^2 / (2*sigma^2));
pc    = min(pc, 1.0);
end
