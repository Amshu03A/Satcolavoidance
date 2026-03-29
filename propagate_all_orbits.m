function [sat_states, debris_states] = propagate_all_orbits(satellites, debris, params)
% Propagate all satellites and debris objects through the simulation window.
%
% Uses Keplerian two-body propagation with J2 secular corrections
% for RAAN regression and perigee rotation.
%
% Returns:
%   sat_states    [N_sat x N_steps x 6]  ECI (x,y,z [m], vx,vy,vz [m/s])
%   debris_states [N_deb x N_steps x 6]

n_sats  = length(satellites);
n_deb   = length(debris);
n_steps = params.N_steps;
t_vec   = params.t_vec;

sat_states    = zeros(n_sats, n_steps, 6);
debris_states = zeros(n_deb,  n_steps, 6);

for s = 1:n_sats
    elems = pull_elements(satellites(s));
    for k = 1:n_steps
        elems_t = advance_j2(elems, t_vec(k), params);
        sat_states(s, k, :) = elements_to_eci(elems_t, params);
    end
end

for d = 1:n_deb
    elems = pull_elements(debris(d));
    for k = 1:n_steps
        elems_t = advance_j2(elems, t_vec(k), params);
        debris_states(d, k, :) = elements_to_eci(elems_t, params);
    end
end
end

% -------------------------------------------------------------------------
function elems = pull_elements(obj)
% Pull orbital elements out of a satellite or debris struct.
elems.a     = obj.a;
elems.e     = obj.e;
elems.i     = obj.i;
elems.RAAN  = obj.RAAN;
elems.omega = obj.omega;
elems.M0    = obj.M0;
end

% -------------------------------------------------------------------------
function elems_t = advance_j2(elems, t, params)
% Advance orbital elements to time t, accounting for J2 secular drift.
%
% Drift rates:
%   dRAAN/dt  = -3/2 * n * J2 * (Re/p)^2 * cos(i)
%   domega/dt = +3/4 * n * J2 * (Re/p)^2 * (5*cos^2(i) - 1)

mu = params.mu;
J2 = params.J2;
Re = params.Re;

a = elems.a;
e = elems.e;
i = elems.i;

n = sqrt(mu / a^3);
p = a * (1 - e^2);

raan_dot  = -1.5  * n * J2 * (Re/p)^2 * cos(i);
omega_dot =  0.75 * n * J2 * (Re/p)^2 * (5*cos(i)^2 - 1);

elems_t.a     = a;
elems_t.e     = e;
elems_t.i     = i;
elems_t.RAAN  = elems.RAAN  + raan_dot  * t;
elems_t.omega = elems.omega + omega_dot * t;
elems_t.M     = elems.M0    + n         * t;
end

% -------------------------------------------------------------------------
function state = elements_to_eci(elems, params)
% Convert orbital elements to ECI Cartesian state [x,y,z,vx,vy,vz].

mu = params.mu;
a  = elems.a;
e  = elems.e;
i  = elems.i;
W  = elems.RAAN;
w  = elems.omega;
M  = elems.M;

E  = solve_kepler(M, e);
nu = 2 * atan2(sqrt(1+e)*sin(E/2), sqrt(1-e)*cos(E/2));
r  = a * (1 - e*cos(E));
p  = a * (1 - e^2);

r_pf = [r*cos(nu); r*sin(nu); 0];
v_pf = sqrt(mu/p) * [-sin(nu); e+cos(nu); 0];

R = rot3(-W) * rot1(-i) * rot3(-w);

state = [R*r_pf; R*v_pf]';
end

% -------------------------------------------------------------------------
function E = solve_kepler(M, e)
% Newton-Raphson solve for Kepler's equation.
M = mod(M, 2*pi);
E = M;
for it = 1:50
    dE = (M - E + e*sin(E)) / (1 - e*cos(E));
    E  = E + dE;
    if abs(dE) < 1e-12, break; end
end
end

function R = rot1(a)
R = [1 0 0; 0 cos(a) sin(a); 0 -sin(a) cos(a)];
end

function R = rot3(a)
R = [cos(a) sin(a) 0; -sin(a) cos(a) 0; 0 0 1];
end
