function visualize_simulation(sat_states, debris_states, conjunctions, ...
                              maneuver_log, fuel_report, params)
% Generate all plots for the conjunction avoidance sim.
%
% Fig 1 - 3D orbital view (Earth + sat + debris trajectories)
% Fig 2 - Conjunction timeline (top 20 closest approaches)
% Fig 3 - Risk heatmap (sat vs debris min distance)
% Fig 4 - Fuel and delta-V budget per satellite
% Fig 5 - Probability of collision distribution
% Fig 6 - Maneuver timeline

Re     = params.Re / 1e3;   % km
t_days = params.t_vec / 86400;
n_sats = size(sat_states, 1);
n_deb  = size(debris_states, 1);

%% Figure 1: 3D Orbital View
figure('Name','3D Orbital View','NumberTitle','off','Color','k','Position',[50 50 900 750]);
ax = axes; hold on; grid on; axis equal;
ax.Color = 'k'; ax.GridColor = [0.3 0.3 0.3];
ax.XColor = 'w'; ax.YColor = 'w'; ax.ZColor = 'w';

[xe,ye,ze] = sphere(40);
surf(Re*xe, Re*ye, Re*ze, 'FaceColor',[0.1 0.3 0.6],...
    'EdgeColor','none','FaceAlpha',0.8,'DisplayName','Earth');

pal = lines(n_sats);
for s = 1:min(n_sats, 12)
    x = squeeze(sat_states(s,:,1))/1e3;
    y = squeeze(sat_states(s,:,2))/1e3;
    z = squeeze(sat_states(s,:,3))/1e3;
    plot3(x, y, z, '-', 'Color', pal(s,:), 'LineWidth', 1.2, ...
        'DisplayName', sprintf('SAT %d', s));
    plot3(x(1), y(1), z(1), 'o', 'Color', pal(s,:), 'MarkerSize', 6, ...
        'MarkerFaceColor', pal(s,:), 'HandleVisibility','off');
end

n_show_deb = min(n_deb, 30);
for d = 1:n_show_deb
    x = squeeze(debris_states(d,:,1))/1e3;
    y = squeeze(debris_states(d,:,2))/1e3;
    z = squeeze(debris_states(d,:,3))/1e3;
    plot3(x, y, z, '--', 'Color',[0.8 0.3 0.1 0.4], 'LineWidth', 0.5,...
        'HandleVisibility','off');
end
plot3(NaN,NaN,NaN,'--','Color',[0.8 0.3 0.1],'DisplayName','Debris');

title('LEO Constellation + Debris Field','Color','w','FontSize',14);
xlabel('X [km]','Color','w'); ylabel('Y [km]','Color','w'); zlabel('Z [km]','Color','w');
legend('show','TextColor','w','Color','k','Location','northeast');
view(30, 25);

%% Figure 2: Conjunction Timeline
if ~isempty(conjunctions) && isstruct(conjunctions) && isfield(conjunctions,'min_dist_km')
    figure('Name','Conjunction Events','NumberTitle','off','Position',[100 50 900 500]);

    n_show  = min(length(conjunctions), 20);
    [~, si] = sort([conjunctions.min_dist_km]);
    si      = si(1:n_show);

    close_dists = [conjunctions(si).min_dist_km];
    tca_hrs     = [conjunctions(si).tca_time_s] / 3600;

    subplot(1,2,1);
    barh(1:n_show, close_dists, 'FaceColor',[0.2 0.6 0.8]);
    hold on;
    xline(params.danger_dist_km, 'r--', 'LineWidth', 2, 'Label','Danger');
    xline(params.warn_dist_km,   'y--', 'LineWidth', 1.5, 'Label','Warning');
    xlabel('Min Separation [km]'); ylabel('Conjunction Event #');
    title('Closest Approach Distances (Top 20)');
    grid on; set(gca,'YDir','reverse');

    subplot(1,2,2);
    scatter(tca_hrs, close_dists, 60, close_dists, 'filled');
    colormap(flipud(hot)); colorbar;
    xlabel('Time of Closest Approach [hours]');
    ylabel('Miss Distance [km]');
    title('Conjunction TCA Timeline');
    yline(params.danger_dist_km, 'r--', 'LineWidth',2);
    yline(params.warn_dist_km,   'y--', 'LineWidth',1.5);
    grid on;
end

%% Figure 3: Satellite-Debris Risk Heatmap
figure('Name','Risk Heatmap','NumberTitle','off','Position',[150 50 800 500]);
risk_mat = inf(n_sats, n_deb);
if ~isempty(conjunctions) && isfield(conjunctions,'sat_idx')
    for c = 1:length(conjunctions)
        s = conjunctions(c).sat_idx;
        d = conjunctions(c).deb_idx;
        risk_mat(s,d) = min(risk_mat(s,d), conjunctions(c).min_dist_km);
    end
end
risk_plot = risk_mat;
risk_plot(isinf(risk_plot)) = params.warn_dist_km;
imagesc(risk_plot, [0 params.warn_dist_km]);
colormap(flipud(hot)); colorbar;
xlabel('Debris Object Index'); ylabel('Satellite Index');
title('Min Separation Heatmap [km]  (red = closest)');
yticks(1:n_sats);
yticklabels(arrayfun(@(x) sprintf('IMG-%02d',x), 1:n_sats,'UniformOutput',false));
grid off;

%% Figure 4: Fuel Budget
figure('Name','Fuel Budget','NumberTitle','off','Position',[100 80 1100 750]);

if isfield(fuel_report,'per_sat')
    ps    = fuel_report.per_sat;
    n     = length(ps);
    names = {ps.name};
    dvs   = [ps.dv_used];
    spent = [ps.fuel_used];
    left  = [ps.fuel_left];
    burns = [ps.n_maneuvers];
    lives = [ps.est_life_yr];

    ax1 = subplot(2,2,1);
    c_dv = repmat([0.75 0.75 0.75], n, 1);
    c_dv(dvs > 0, :) = repmat([0.2 0.5 0.9], sum(dvs>0), 1);
    b1 = bar(dvs, 'FaceColor','flat');
    b1.CData = c_dv;
    ylabel('\DeltaV Used [m/s]');
    title('\DeltaV Consumed per Satellite','FontSize',11,'FontWeight','bold');
    xticks(1:n); xticklabels(names); xtickangle(40);
    ylim([0, max(dvs)*1.3 + 0.1]);
    for i = 1:n
        if dvs(i) > 0
            text(i, dvs(i) + max(dvs)*0.04, sprintf('%.2f',dvs(i)), ...
                'HorizontalAlignment','center','FontSize',7.5,'Color',[0.1 0.1 0.6]);
        end
    end
    grid on; box off;

    ax2 = subplot(2,2,2);
    c_burns = repmat([0.75 0.75 0.75], n, 1);
    c_burns(burns > 0, :) = repmat([0.85 0.33 0.1], sum(burns>0), 1);
    b2 = bar(burns, 'FaceColor','flat');
    b2.CData = c_burns;
    ylabel('Number of Maneuvers');
    title('Avoidance Maneuvers per Satellite','FontSize',11,'FontWeight','bold');
    xticks(1:n); xticklabels(names); xtickangle(40);
    ylim([0, max(burns)*1.4 + 0.5]);
    yticks(0:1:max(burns)+1);
    for i = 1:n
        if burns(i) > 0
            text(i, burns(i) + 0.05, sprintf('%d',burns(i)), ...
                'HorizontalAlignment','center','FontSize',9,'FontWeight','bold');
        end
    end
    grid on; box off;

    ax3 = subplot(2,2,3);
    c_fuel = repmat([0.75 0.75 0.75], n, 1);
    c_fuel(spent > 0, :) = repmat([0.47 0.67 0.19], sum(spent>0), 1);
    b3 = bar(spent*1000, 'FaceColor','flat');   % grams for readability
    b3.CData = c_fuel;
    ylabel('Fuel Consumed [g]');
    title('Fuel Used per Satellite  (grams)','FontSize',11,'FontWeight','bold');
    xticks(1:n); xticklabels(names); xtickangle(40);
    ylim([0, max(spent*1000)*1.4 + 1]);
    for i = 1:n
        if spent(i) > 0
            text(i, spent(i)*1000 + max(spent*1000)*0.04, ...
                sprintf('%.1fg', spent(i)*1000), ...
                'HorizontalAlignment','center','FontSize',7.5,'Color',[0.1 0.4 0.1]);
        end
    end
    grid on; box off;

    ax4 = subplot(2,2,4);
    full_tank = left + spent;
    pct_left  = (left ./ full_tank) * 100;
    c_life = repmat([0.2 0.6 0.8], n, 1);
    c_life(spent > 0, :) = repmat([0.9 0.6 0.1], sum(spent>0), 1);
    b4 = bar(pct_left, 'FaceColor','flat');
    b4.CData = c_life;
    ylabel('Fuel Remaining [%]');
    title('Fuel Remaining (% of Initial)','FontSize',11,'FontWeight','bold');
    xticks(1:n); xticklabels(names); xtickangle(40);
    ylim([98.5, 100.5]);
    yline(100,'k--','LineWidth',0.8);
    for i = 1:n
        text(i, pct_left(i) - 0.08, sprintf('%.2f%%', pct_left(i)), ...
            'HorizontalAlignment','center','FontSize',7,'Color',[0.2 0.2 0.2]);
    end
    grid on; box off;

    sgtitle('Satellite Fuel Budget & Maneuver Summary','FontSize',13,'FontWeight','bold');
end

%% Figure 5: Pc Distribution
figure('Name','Pc Distribution','NumberTitle','off','Position',[250 50 700 450]);
if ~isempty(conjunctions) && isfield(conjunctions,'Pc')
    pc_vals = [conjunctions.Pc];
    pc_vals = pc_vals(pc_vals > 0);
    histogram(log10(pc_vals + 1e-12), 20, 'FaceColor',[0.8 0.2 0.2],'EdgeColor','w');
    xlabel('log_{10}(Probability of Collision)');
    ylabel('Count');
    title('Distribution of Conjunction Pc Values');
    xline(log10(params.Pc_threshold),'b--','LineWidth',2,...
        'Label',sprintf('Threshold (%.0e)',params.Pc_threshold));
    grid on;
end

%% Figure 6: Maneuver Timeline
if ~isempty(maneuver_log) && isstruct(maneuver_log) && isfield(maneuver_log,'man_time_s')
    figure('Name','Maneuver Timeline','NumberTitle','off','Position',[300 50 900 400]);

    burn_times = [maneuver_log.man_time_s] / 3600;
    burn_dv    = [maneuver_log.delta_v_ms];
    burn_sat   = [maneuver_log.sat_idx];

    scatter(burn_times, burn_sat, burn_dv*80, burn_dv, 'filled');
    colorbar; colormap(hot);
    xlabel('Time [hours]'); ylabel('Satellite Index');
    title('Avoidance Maneuver Timeline (bubble size = \DeltaV magnitude)');
    yticks(1:n_sats);
    yticklabels(arrayfun(@(x) sprintf('IMG-%02d',x), 1:n_sats,'UniformOutput',false));
    grid on;
end

fprintf('      All plots generated.\n');
end
