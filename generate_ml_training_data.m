function [X, Y, label_names] = generate_ml_training_data(num_samples)
% Generate labelled satellite-debris pairs for ML training.
%
% Each sample is a random sat + debris with orbital elements.
% We propagate both, measure their closest approach, and label:
%   - min approach distance [km]
%   - probability of collision
%   - risk class: 0=Safe, 1=Warning, 2=Critical
%   - suggested burn size [m/s]
%
% Features X (12 cols per row):
%   [sat_a, sat_e, sat_i, sat_RAAN, sat_w, sat_M0,
%    deb_a, deb_e, deb_i, deb_RAAN, deb_w, deb_M0]
%   (axes in km, angles in degrees)
%
% Labels Y (4 cols):
%   [min_dist_km, Pc, risk_class, rec_dv_ms]

if nargin < 1, num_samples = 5000; end

fprintf('Generating %d training samples...\n', num_samples);

p.mu      = 3.986004418e14;
p.Re      = 6371e3;
p.J2      = 1.08263e-3;
p.dt      = 30;
p.T_total = 1 * 86400;
p.t_vec   = 0:p.dt:p.T_total;
p.N_steps = length(p.t_vec);
p.warn_dist_km   = 50;
p.danger_dist_km = 10;
p.Pc_threshold   = 1e-4;

X = zeros(num_samples, 12);
Y = zeros(num_samples, 4);

rng(123);

for k = 1:num_samples
    if mod(k, 500) == 0
        fprintf('  Sample %d / %d\n', k, num_samples);
    end

    % Cycle through classes so we don't end up with almost no Critical samples
    target = mod(k-1, 3);   % 0=Safe, 1=Warning, 2=Critical

    % Random satellite in LEO imaging band (500-600 km)
    sat_alt   = 500 + rand()*100;
    sat.a     = p.Re + sat_alt*1e3;
    sat.e     = rand()*0.02;
    sat.i     = deg2rad(80 + rand()*20);
    sat.RAAN  = deg2rad(rand()*360);
    sat.omega = deg2rad(rand()*360);
    sat.M0    = deg2rad(rand()*360);

    % Debris placement tuned to hit the desired class
    switch target
        case 2   % Critical: nearly same alt + RAAN, tiny inclination offset
            deb_alt  = sat_alt + (rand()-0.5)*5;
            deb_incl = rad2deg(sat.i) + (rand()-0.5)*2;
            deb_RAAN = rad2deg(sat.RAAN) + (rand()-0.5)*3;
            deb_M0   = rad2deg(sat.M0) + 170 + rand()*20;

        case 1   % Warning: moderate separation
            deb_alt  = sat_alt + (rand()-0.5)*20;
            deb_incl = rad2deg(sat.i) + (rand()-0.5)*8;
            deb_RAAN = rad2deg(sat.RAAN) + (rand()-0.5)*15;
            deb_M0   = rad2deg(sat.M0) + 150 + rand()*60;

        otherwise % Safe: push debris well away
            deb_alt  = sat_alt + 30 + rand()*80;
            deb_incl = rad2deg(sat.i) + 15 + rand()*40;
            deb_RAAN = rad2deg(sat.RAAN) + rand()*360;
            deb_M0   = rad2deg(sat.M0) + rand()*360;
    end

    deb.a           = p.Re + deb_alt*1e3;
    deb.e           = rand()*0.03;
    deb.i           = deg2rad(mod(deb_incl, 180));
    deb.RAAN        = deg2rad(mod(deb_RAAN, 360));
    deb.omega       = deg2rad(rand()*360);
    deb.M0          = deg2rad(mod(deb_M0, 360));
    deb.hard_body_r = 0.1 + rand()*0.9;

    X(k,:) = [sat.a/1e3, sat.e, rad2deg(sat.i), ...
               rad2deg(sat.RAAN), rad2deg(sat.omega), rad2deg(sat.M0), ...
               deb.a/1e3, deb.e, rad2deg(deb.i), ...
               rad2deg(deb.RAAN), rad2deg(deb.omega), rad2deg(deb.M0)];

    sat_traj = propagate_one(sat, p);
    deb_traj = propagate_one(deb, p);

    gap      = sat_traj(:,1:3) - deb_traj(:,1:3);
    dist_km  = sqrt(sum(gap.^2, 2)) / 1e3;
    [closest, tca] = min(dist_km);

    rel_spd  = norm(sat_traj(tca,4:6) - deb_traj(tca,4:6));
    pc       = chan_pc(closest*1e3, 1.5, deb.hard_body_r, rel_spd);

    if closest < p.danger_dist_km || pc > p.Pc_threshold
        risk  = 2;
        rec_dv = 0.5 + closest/p.danger_dist_km * 0.5;
    elseif closest < p.warn_dist_km
        risk  = 1;
        rec_dv = 0.1 + (p.warn_dist_km - closest)/p.warn_dist_km * 0.4;
    else
        risk  = 0;
        rec_dv = 0;
    end

    Y(k,:) = [closest, pc, risk, rec_dv];
end

label_names = {'min_dist_km','Pc','risk_class','recommended_dv_ms'};
fprintf('Done. Class distribution: Safe=%d | Warning=%d | Critical=%d\n',...
    sum(Y(:,3)==0), sum(Y(:,3)==1), sum(Y(:,3)==2));
end

% -----------------------------------------------------------------------
function traj = propagate_one(obj, p)
% Propagate one object and return [N_steps x 6] ECI states.
mu = p.mu; Re = p.Re; J2 = p.J2;
a = obj.a; e = obj.e; i = obj.i;
n  = sqrt(mu/a^3);
lp = a*(1-e^2);
dR = -1.5*n*J2*(Re/lp)^2*cos(i);
dw =  0.75*n*J2*(Re/lp)^2*(5*cos(i)^2-1);
N  = p.N_steps;
traj = zeros(N,6);
for k = 1:N
    t    = p.t_vec(k);
    RAAN = obj.RAAN + dR*t;
    omg  = obj.omega + dw*t;
    M    = obj.M0 + n*t;
    E    = solve_kepler(M, e);
    nu   = 2*atan2(sqrt(1+e)*sin(E/2), sqrt(1-e)*cos(E/2));
    r    = a*(1-e*cos(E));
    rp   = [r*cos(nu); r*sin(nu); 0];
    vp   = sqrt(mu/lp)*[-sin(nu); e+cos(nu); 0];
    R    = rot3(-RAAN)*rot1(-i)*rot3(-omg);
    traj(k,:) = [R*rp; R*vp]';
end
end

function E = solve_kepler(M, e)
M = mod(M, 2*pi); E = M;
for it = 1:50
    dE = (M - E + e*sin(E)) / (1 - e*cos(E));
    E  = E + dE;
    if abs(dE) < 1e-12, break; end
end
end

function R = rot1(a); R = [1 0 0; 0 cos(a) sin(a); 0 -sin(a) cos(a)]; end
function R = rot3(a); R = [cos(a) sin(a) 0; -sin(a) cos(a) 0; 0 0 1]; end

function pc = chan_pc(d, rs, rd, ~)
sig = 100; A = pi*(rs+rd)^2;
pc  = min((A/(2*pi*sig^2))*exp(-d^2/(2*sig^2)), 1);
end
