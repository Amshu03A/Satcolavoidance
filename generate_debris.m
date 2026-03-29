function debris = generate_debris(params)
% Build a synthetic debris field in LEO.
%
% Three rough size buckets:
%   Large  - rocket body fragments (~30 cm to 1 m)
%   Medium - collision chips (9-30 cm)
%   Small  - paint flecks and tiny bits (< 9 cm)
%
% Spread across 530-580 km altitude so debris actually shares
% the orbital shell with the satellites.

rng(42);   % reproducible

Re         = params.Re;
num_pieces = 150;

alt_lo = 530;
alt_hi = 580;

% 60% small, 30% medium, 10% large
size_pdf    = [0.6, 0.3, 0.1];
size_ranges = {[0.01, 0.09], [0.09, 0.30], [0.30, 1.0]};   % diameter [m]
rough_mass  = [0.2, 5, 50];   % kg per bucket

debris = repmat(struct('name','','a',0,'e',0,'i',0,'RAAN',0,...
                       'omega',0,'M0',0,'size_m',0,'mass_kg',0,...
                       'hard_body_r',0), 1, num_pieces);

for d = 1:num_pieces
    alt_km = alt_lo + rand() * (alt_hi - alt_lo);

    % Mostly near-circular, a few stragglers with some eccentricity
    if rand() < 0.85
        ecc = rand() * 0.01;
    else
        ecc = 0.01 + rand() * 0.08;
    end

    % Inclination: bias toward sun-sync band so debris crosses our sats often
    r_draw = rand();
    if r_draw < 0.50
        incl = 85 + rand() * 15;   % sun-sync band
    elseif r_draw < 0.80
        incl = 70 + rand() * 15;   % high inclination
    else
        incl = 20 + rand() * 50;   % lower inclination
    end

    % Pick size bucket
    r = rand();
    if r < size_pdf(1)
        bucket = 1;
    elseif r < size_pdf(1) + size_pdf(2)
        bucket = 2;
    else
        bucket = 3;
    end
    lo = size_ranges{bucket}(1);
    hi = size_ranges{bucket}(2);
    sz = lo + rand() * (hi - lo);

    debris(d).name        = sprintf('DEB-%03d', d);
    debris(d).a           = Re + alt_km*1e3;
    debris(d).e           = ecc;
    debris(d).i           = deg2rad(incl);
    debris(d).RAAN        = deg2rad(rand()*360);
    debris(d).omega       = deg2rad(rand()*360);
    debris(d).M0          = deg2rad(rand()*360);
    debris(d).size_m      = sz;
    debris(d).mass_kg     = rough_mass(bucket) * (1 + rand());
    debris(d).hard_body_r = sz / 2;
end

fprintf('      Debris breakdown: %d small | %d medium | %d large\n', ...
    sum(arrayfun(@(d) d.size_m < 0.09, debris)), ...
    sum(arrayfun(@(d) d.size_m >= 0.09 && d.size_m < 0.30, debris)), ...
    sum(arrayfun(@(d) d.size_m >= 0.30, debris)));
end
