function [] = PlotCOMVel(tout,traj,xout, LL, cbf_)

Fs = 10000;
Fc = 100;
order = 4;
[b, a] = butter(order, Fc/(Fs/2), 'low')

% Define colors
dark_orange = [0.8902, 0.4745, 0.4510];
orange = [0.9882, 0.7059, 0.4078];
light_blue = [0.0431, 0.2941, 0.4627];
dark_blue = [0.3020, 0.7137, 0.7529];

% Apply the Butterworth filter to each column of xout
for col = 1:3
    LL.dhd(:, col) = filtfilt(b, a, LL.dhd(:, col));
end

for col = 4:6
    traj(:, col) = filtfilt(b, a, traj(:, col));
end

% Velocity of the COM
FS = 18;
LW = 3;

% Create a figure
figureHandle = figure;

% Desired size in inches
widthInches = 6*1.7;  % Width of the figure
heightInches = 2.5*1.2; % Height of the figure

% Set the PaperUnits to inches and set the PaperPosition to desired size
set(figureHandle, 'Units', 'inches');
set(figureHandle, 'Position', [0, 0, widthInches, heightInches]);



% subplot(2,1,1);
% yyaxis left;
plot(tout,traj(:,4), 'LineStyle', '-', 'color', dark_blue,'LineWidth',LW, 'DisplayName', '$\dot{x}^{\textrm{ref}}$');
hold on
% plot(tout,xout(:,4),'b-.','LineWidth',1);
% plot(LL.time,LL.vel(:, 1),'c--','LineWidth',1);
plot(LL.time,LL.dhd(:, 1)*1.0, 'LineStyle', '-', 'color', light_blue, 'LineWidth',LW, 'DisplayName', '$\dot{x}^{\star}$');
legend('FontSize', FS, 'NumColumns', 2, 'interpreter', 'latex', 'Location', 'best') %'qk_','LL',
ylabel('x', 'FontSize', FS, 'interpreter', 'latex')
title('Optimal Trajectory from Middle-Level HOBF-NMPC - Agent 2', 'FontSize', FS, 'interpreter', 'latex')

xlim([2, 8.5]);
% ylim([-0.5, 0.5])
% subplot(3,1,2);
plot(tout,traj(:,5), 'LineStyle', '-', 'color', dark_orange, 'LineWidth',LW, 'DisplayName', '$\dot{y}^{\textrm{ref}}$');
hold on
% plot(tout,xout(:,5),'b-.','LineWidth',1);
% plot(LL.time,LL.vel(:, 2),'c--','LineWidth',1);
plot(LL.time,LL.dhd(:, 2), 'LineStyle', '-', 'color', orange, 'LineWidth',LW, 'DisplayName', '$\dot{y}^{\star}$');
legend('FontSize', FS, 'interpreter', 'latex', 'Location', 'best')
ylabel('Velocity (m/s)', 'FontSize', FS, 'interpreter', 'latex')
xlabel('Time (s)', 'FontSize', FS, 'interpreter', 'latex')
ax = gca;
ax.FontSize = FS;
% Set the interpreter to LaTeX for axis labels and title
ax.XLabel.Interpreter = 'latex';
ax.YLabel.Interpreter = 'latex';
ax.ZLabel.Interpreter = 'latex';  % If you have a 3D plot
ax.Title.Interpreter = 'latex';
ax.TickLabelInterpreter = 'latex';
ax.LineWidth = 1.5;

% yyaxis right;
% cbf_ = cbf_-min(cbf_)-0.3;
% cbf_(find(cbf_<0)) = 0;
% 
% plot(tout, cbf_', 'LineWidth', LW, 'LineStyle','--', 'Color', 'k')

gcf;
print('case_a_agent2', '-djpeg', '-r600');
% 
% subplot(3,1,3);
% plot(tout,traj(:,6),'k-.','LineWidth',LW);
% hold on
% % plot(LL.time,LL.vel(:, 3),'b-.','LineWidth',1);
% plot(LL.time,LL.dhd(:, 3),'r-','LineWidth',LW);
% legend('HL', 'NMPC')
% ylabel('z')
% xlabel('Time (sec)')




end

