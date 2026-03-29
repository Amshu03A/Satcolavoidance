%% =========================================================================
%  Validation & Test Suite
%  Checks the Kepler solver, OE->ECI, J2 drift, rocket equation, Pc, and
%  a quick end-to-end smoke test.
%% =========================================================================

clc; clear; close all;
fprintf('==============================================\n');
fprintf('   Validation & Test Suite\n');
fprintf('==============================================\n\n');

p.mu = 3.986004418e14;
p.Re = 6371e3;
p.J2 = 1.08263e-3;
p.dt = 10;

passed = 0;
failed = 0;

%% Test 1: Kepler solver round-trip
fprintf('[TEST 1] Kepler Equation Solver...\n');
test_M = [0, pi/4, pi/2, pi, 1.5*pi, 2*pi-0.01];
test_e = [0, 0.1,  0.3,  0.5, 0.7,   0.01];
for k = 1:length(test_M)
    E      = kepler_local(test_M(k), test_e(k));
    M_back = E - test_e(k)*sin(E);
    err    = abs(mod(M_back, 2*pi) - mod(test_M(k), 2*pi));
    if err < 1e-10
        fprintf('  PASS: M=%.4f, e=%.2f -> E=%.6f (residual=%.2e)\n', test_M(k), test_e(k), E, err);
        passed = passed + 1;
    else
        fprintf('  FAIL: residual %.2e exceeds tolerance\n', err);
        failed = failed + 1;
    end
end

%% Test 2: Orbital radius at epoch
fprintf('\n[TEST 2] Orbital radius from OE2ECI...\n');
oe.a = 6778e3; oe.e = 0.001; oe.i = deg2rad(51.6);
oe.RAAN = deg2rad(40); oe.omega = deg2rad(0); oe.M = deg2rad(0);
st        = oe2eci_local(oe, p);
r_got     = norm(st(1:3));
r_want    = oe.a * (1 - oe.e^2) / (1 + oe.e);
err       = abs(r_got - r_want) / r_want;
if err < 0.01
    fprintf('  PASS: r=%.1f km (expected ~%.1f km, err=%.4f%%)\n',...
        r_got/1e3, r_want/1e3, err*100);
    passed = passed + 1;
else
    fprintf('  FAIL: r error = %.4f%%\n', err*100);
    failed = failed + 1;
end

%% Test 3: Energy conservation around a circular orbit
fprintf('\n[TEST 3] Orbital Energy Conservation...\n');
oe.a = 6928e3; oe.e = 0.0; oe.i = deg2rad(45);
oe.RAAN = 0; oe.omega = 0;
energies = zeros(1,360);
for deg = 1:360
    oe.M = deg2rad(deg);
    st   = oe2eci_local(oe, p);
    r    = norm(st(1:3));
    v    = norm(st(4:6));
    energies(deg) = 0.5*v^2 - p.mu/r;
end
spread = (max(energies) - min(energies)) / abs(mean(energies));
if spread < 1e-6
    fprintf('  PASS: Energy variation = %.2e (< 1e-6)\n', spread);
    passed = passed + 1;
else
    fprintf('  FAIL: Energy variation = %.2e\n', spread);
    failed = failed + 1;
end

%% Test 4: J2 RAAN drift goes westward for prograde orbit
fprintf('\n[TEST 4] J2 RAAN drift direction (prograde orbit -> westward)...\n');
oe.a = p.Re + 550e3; oe.e = 0.001; oe.i = deg2rad(51.6);
oe.RAAN = 0; oe.omega = 0; oe.M0 = 0;
lp         = oe.a*(1 - oe.e^2);
n          = sqrt(p.mu/oe.a^3);
raan_rate  = -1.5 * n * p.J2 * (p.Re/lp)^2 * cos(oe.i);
if raan_rate < 0
    fprintf('  PASS: dRAAN/dt = %.4e rad/s (negative = westward)\n', raan_rate);
    fprintf('        = %.4f deg/day\n', rad2deg(raan_rate)*86400);
    passed = passed + 1;
else
    fprintf('  FAIL: Expected negative dRAAN/dt for prograde orbit\n');
    failed = failed + 1;
end

%% Test 5: Tsiolkovsky rocket equation
fprintf('\n[TEST 5] Rocket equation...\n');
m_dry = 150; fuel = 20; isp = 220; g0 = 9.80665;
dv_calc = isp * g0 * log((m_dry+fuel)/m_dry);
dv_ref  = 220 * 9.80665 * log(170/150);
err     = abs(dv_calc - dv_ref);
if err < 0.01
    fprintf('  PASS: Total dV = %.3f m/s (expected %.3f m/s, err=%.4f m/s)\n',...
        dv_calc, dv_ref, err);
    passed = passed + 1;
else
    fprintf('  FAIL: dV error = %.4f m/s\n', err);
    failed = failed + 1;
end

%% Test 6: Pc decreases with miss distance
fprintf('\n[TEST 6] Pc monotonically decreases with miss distance...\n');
dists  = 100:100:5000;
pc_arr = arrayfun(@(d) pc_local(d, 1.5, 0.5, 7800), dists);
if all(diff(pc_arr) <= 0)
    fprintf('  PASS: Pc decreases monotonically. Pc(100m)=%.2e, Pc(5km)=%.2e\n',...
        pc_arr(1), pc_arr(end));
    passed = passed + 1;
else
    fprintf('  FAIL: Pc not monotonically decreasing\n');
    failed = failed + 1;
end

%% Test 7: Quick end-to-end smoke test
fprintf('\n[TEST 7] End-to-end smoke test (mini run)...\n');
try
    mini = p;
    mini.dt = 60; mini.T_days = 0.1;
    mini.T_total = mini.T_days*86400;
    mini.t_vec   = 0:mini.dt:mini.T_total;
    mini.N_steps = length(mini.t_vec);
    mini.warn_dist_km   = 10;
    mini.danger_dist_km = 2;
    mini.Pc_threshold   = 1e-4;

    sats = generate_satellites(mini);
    debs = generate_debris(mini);
    [ss, ds] = propagate_all_orbits(sats(1:2), debs(1:5), mini);

    assert(size(ss,1)==2, 'Wrong sat count');
    assert(size(ds,1)==5, 'Wrong debris count');
    assert(size(ss,3)==6, 'Wrong state dimension');

    fprintf('  PASS: Mini run done. ss=[%dx%dx%d], ds=[%dx%dx%d]\n',...
        size(ss,1),size(ss,2),size(ss,3),...
        size(ds,1),size(ds,2),size(ds,3));
    passed = passed + 1;
catch ME
    fprintf('  FAIL: %s\n', ME.message);
    failed = failed + 1;
end

%% Results
fprintf('\n==============================================\n');
fprintf('  Results: %d passed | %d failed\n', passed, failed);
if failed == 0
    fprintf('  ALL TESTS PASSED\n');
else
    fprintf('  SOME TESTS FAILED - review output above\n');
end
fprintf('==============================================\n');

%% Local helpers

function E = kepler_local(M, e)
M = mod(M, 2*pi); E = M;
for it = 1:100
    dE = (M - E + e*sin(E))/(1 - e*cos(E));
    E  = E + dE;
    if abs(dE) < 1e-12, break; end
end
end

function state = oe2eci_local(oe, p)
mu = p.mu;
a  = oe.a; e = oe.e; i = oe.i;
W  = oe.RAAN; w = oe.omega;
if isfield(oe,'M'), M = oe.M; else, M = oe.M0; end
E  = kepler_local(M, e);
nu = 2*atan2(sqrt(1+e)*sin(E/2), sqrt(1-e)*cos(E/2));
r  = a*(1 - e*cos(E));
lp = a*(1 - e^2);
r_pf = [r*cos(nu); r*sin(nu); 0];
v_pf = sqrt(mu/lp)*[-sin(nu); e+cos(nu); 0];
R1  = [1 0 0; 0 cos(-i) sin(-i); 0 -sin(-i) cos(-i)];
R3W = [cos(-W) sin(-W) 0; -sin(-W) cos(-W) 0; 0 0 1];
R3w = [cos(-w) sin(-w) 0; -sin(-w) cos(-w) 0; 0 0 1];
R   = R3W * R1 * R3w;
state = [R*r_pf; R*v_pf];
end

function pc = pc_local(miss_dist, r_sat, r_deb, ~)
sig  = 100;
area = pi*(r_sat+r_deb)^2;
pc   = (area/(2*pi*sig^2)) * exp(-miss_dist^2/(2*sig^2));
end
