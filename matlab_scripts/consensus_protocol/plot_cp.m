function [] = plot_cp()




%% With Consensus Protocol
cd consensus_protocol/
    LL0_cp = parseData("agent_0_cp.txt");
    LL1_cp = parseData("agent_1_cp.txt");
cd ..

x0 = LL0_cp.pos(:, 1); % x position of 0th agent
y0 = LL0_cp.pos(:, 2); % y position of 0th agent

x1 = LL1_cp.pos(:, 1); % x position of 1st agent
y1 = LL1_cp.pos(:, 2); % y position of 1st agent

time = LL0_cp.time(:); % time

idxx = 1:10000;

x0 = x0(idxx); y0 = y0(idxx);
x1 = x1(idxx); y1 = y1(idxx);
% Calculate the distance between agent 0 and agent 1


% Calculate the LJ potential
epsilon = 0.5;
sigma = 0.8;
dist_cp = sqrt((x1 - x0).^2 + (y1 - y0).^2);
LJ_potential_cp = 1e9 * 4 * epsilon * ( (sigma ./ dist_cp).^12 - (sigma ./ dist_cp).^6 );

cd consensus_protocol/
    LL0_wocp = parseData("agent_0_wocp.txt");
    LL1_wocp = parseData("agent_1_wocp.txt");
cd ..

x0 = LL0_wocp.pos(:, 1); % x position of 0th agent
y0 = LL0_wocp.pos(:, 2); % y position of 0th agent

x1 = LL1_wocp.pos(:, 1); % x position of 1st agent
y1 = LL1_wocp.pos(:, 2); % y position of 1st agent

time = LL0_wocp.time(:); % time

idxx = 1:10000;

x0 = x0(idxx); y0 = y0(idxx);
x1 = x1(idxx); y1 = y1(idxx);
% Calculate the distance between agent 0 and agent 1


% Calculate the LJ potential
epsilon = 0.5;
sigma = 0.8;

dist_wocp = sqrt((x1 - x0).^2 + (y1 - y0).^2);
LJ_potential_wocp = 0 * 4 * epsilon * ( (sigma ./ dist_wocp).^12 - (sigma ./ dist_wocp).^6 );



% Plot the LJ potential over time
figure(8); clf;
subplot(1, 2, 1);
plot(time(idxx), LJ_potential_cp, 'LineWidth', 2.5); hold on;
plot(time(idxx), LJ_potential_wocp, 'LineWidth', 2.5); grid on;
lgd = legend('$\kappa = 1e9$', '$\kappa = 0$');
lgd.FontSize = 18; lgd.Interpreter = 'latex';
ax = gca; ax.FontSize = 18;
ax.XLabel.Interpreter = 'latex';
ax.YLabel.Interpreter = 'latex';
ax.ZLabel.Interpreter = 'latex';  % If you have a 3D plot
ax.Title.Interpreter = 'latex';
ax.TickLabelInterpreter = 'latex';
ax.LineWidth = 1.5;

hold on;
xlabel('Time', 'FontSize', 18, 'Interpreter', 'latex');
ylabel('Consensus cost $\kappa  \Gamma(x^i, x^j)$', 'FontSize', 18, 'Interpreter', 'latex');
title('Consensus Cost between Agent 1 and Agent 2', 'FontSize', 18, 'Interpreter', 'latex');
subplot(1, 2, 2);
plot(time(idxx), (dist_cp - 0.9) + 1, 'LineWidth', 4); hold on;
plot(time(idxx), (dist_wocp - 0.9)*1.0 +1, 'LineWidth', 4)
lgd = legend('$\kappa = 1e9$', '$\kappa = 0$');
lgd.FontSize = 18; lgd.Interpreter = 'latex';
ax = gca; ax.FontSize = 18;
ax.XLabel.Interpreter = 'latex';
ax.YLabel.Interpreter = 'latex';
ax.ZLabel.Interpreter = 'latex';  % If you have a 3D plot
ax.Title.Interpreter = 'latex';
ax.TickLabelInterpreter = 'latex';
ax.LineWidth = 1.5;

ylim([0.5 1.5])
xlabel('Time', 'FontSize', 18, 'Interpreter', 'latex');
ylabel('Distance (m)', 'FontSize', 18, 'Interpreter', 'latex');
title('Distance between Agent 1 and Agent 2', 'FontSize', 18, 'Interpreter', 'latex');
grid on;

gcf;
print('consensus_plot', '-djpeg', '-r600');

end
