function satellites = generate_satellites(params)
% Set up a 12-satellite LEO imaging constellation.
%
% Layout: 3 orbital planes, 4 birds per plane.
% All at ~550 km, sun-synchronous-ish inclination (97.6 deg).
% Planes spaced 120 deg apart in RAAN; satellites evenly spaced within each plane.

Re          = params.Re;
n_planes    = 3;
per_plane   = 4;
total       = n_planes * per_plane;

alt_km   = 550;
a_base   = Re + alt_km*1e3;
ecc      = 0.001;
incl_deg = 97.6;

raan_offsets = linspace(0, 120, n_planes);   % 120 deg between planes

% Pre-fill all fields to avoid MATLAB "dissimilar structures" errors
satellites = repmat(struct('name','','a',0,'e',0,'i',0,'RAAN',0,...
                           'omega',0,'M0',0,'mass',0,'fuel_mass',0,...
                           'Isp',0,'hard_body_r',0,'fuel_used',0,...
                           'maneuvers',0), 1, total);

idx = 0;
for p = 1:n_planes
    for s = 1:per_plane
        idx    = idx + 1;
        m0_deg = (s-1) * (360/per_plane);   % even spacing within plane

        satellites(idx).name        = sprintf('IMG-%02d', idx);
        satellites(idx).a           = a_base + (p-1)*5e3;
        satellites(idx).e           = ecc;
        satellites(idx).i           = deg2rad(incl_deg);
        satellites(idx).RAAN        = deg2rad(raan_offsets(p));
        satellites(idx).omega       = deg2rad(0);
        satellites(idx).M0          = deg2rad(m0_deg);
        satellites(idx).mass        = 150;    % dry mass [kg]
        satellites(idx).fuel_mass   = 20;     % [kg]
        satellites(idx).Isp         = 220;    % [s]
        satellites(idx).hard_body_r = 1.5;    % [m]
        satellites(idx).fuel_used   = 0;
        satellites(idx).maneuvers   = 0;
    end
end
end
