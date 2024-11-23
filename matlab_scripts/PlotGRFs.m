function [] = PlotGRFs(tout,uout)

% Ground reaction force
fplt = figure;
subplot(2,2,1);
plot(tout,uout(:,6),'b'); hold on; plot(tout,uout(:,5),'g'); plot(tout,uout(:,4),'r');
legend('f_z');  % 'f_y', f_x
title('Leg 2 GRFs')
ylabel('f_2 (N)')
xlabel('Time (sec)')

subplot(2,2,2);
plot(tout,uout(:,3),'b'); hold on; plot(tout,uout(:,2),'g'); plot(tout,uout(:,1),'r');
legend('f_z','f_y','f_x')
title('Leg 1 GRFs')
ylabel('f_1 (N)')
xlabel('Time (sec)')

subplot(2,2,3);
plot(tout,uout(:,12),'b'); hold on; plot(tout,uout(:,11),'g'); plot(tout,uout(:,10),'r');
legend('f_z','f_y','f_x')
title('Leg 4 GRFs')
ylabel('f_4 (N)')
xlabel('Time (sec)')

subplot(2,2,4);
plot(tout,uout(:,9),'b'); hold on; plot(tout,uout(:,8),'g'); plot(tout,uout(:,7),'r');
legend('f_z','f_y','f_x')
title('Leg 3 GRFs')
ylabel('f_3 (N)')
xlabel('Time (sec)')

% Velocity of the COM
grfplt = figure;
FS = 18;
LW = 3;

Fs = 10000;
Fc = 70;
order = 4;
[b, a] = butter(order, Fc/(Fs/2), 'low')

% Desired size in inches
widthInches = 6*1.7;  % Width of the figure
heightInches = 2.5*1.2; % Height of the figure
% Define colors
dark_orange = [0.8902, 0.4745, 0.4510];
orange = [0.9882, 0.7059, 0.4078];
light_blue = [0.0431, 0.2941, 0.4627];
dark_blue = [0.3020, 0.7137, 0.7529];

% Set the PaperUnits to inches and set the PaperPosition to desired size
set(grfplt, 'Units', 'inches');
set(grfplt, 'Position', [0, 0, widthInches, heightInches]);


plot(tout-0.14,uout(:,3)*1.6,"Color",light_blue, 'LineWidth',2);
title('Optimal Vertical GRF on Rough Terrain', 'FontSize', FS)

xlim([3 8]);
ylim([0 80]);
% legend('FontSize', FS, 'interpreter', 'latex', 'Location', 'best')
ylabel('Newton (N)', 'FontSize', FS, 'interpreter', 'latex')
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

gcf;
print('rough_GRFs', '-djpeg', '-r600');



end

