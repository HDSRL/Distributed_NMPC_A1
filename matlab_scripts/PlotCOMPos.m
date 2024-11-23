function [] = PlotCOMPos(tout,traj,xout, LL)
% Position of the COM

FS = 16; 
posfig = figure;
subplot(3,1,1);
plot(tout,traj(:,1),'k-.','LineWidth',1.5);
hold on
plot(tout,xout(:,1),'b-.','LineWidth',1);
plot(LL.time,LL.hd(:, 1),'r-','LineWidth',1);
plot(LL.time,LL.pos(:, 1),'c--','LineWidth',1);
legend('HL', 'qk_','hd', 'LL', 'FontSize', FS)
ylabel('x', 'FontSize', FS)
title('COM position', 'FontSize', FS)

subplot(3,1,2);
plot(tout,traj(:,2),'k-.','LineWidth',1.5);
hold on
plot(tout,xout(:,2),'b-.','LineWidth',1);
plot(LL.time,LL.hd(:, 2),'r-','LineWidth',1);
plot(LL.time,LL.pos(:, 2),'c--','LineWidth',1);
legend('HL', 'qk_','hd', 'LL', 'FontSize', FS)
ylabel('y', 'FontSize', FS)

subplot(3,1,3);
plot(tout,traj(:,3),'k-.','LineWidth',1.5);
hold on
plot(tout,xout(:,3),'b-.','LineWidth',1);
plot(LL.time,LL.hd(:, 3),'r-','LineWidth',1);
legend('HL', 'NMPC','Actual', 'FontSize', FS)
ylabel('z', 'FontSize', FS)
xlabel('Time (sec)', 'FontSize', FS)

% figure;
% plot(xout(:, 1), xout(:, 2), 'LineWidth', 2)
% xlabel('(x) position', 'FontSize', FS)
% ylabel('(y) position', 'FontSize', FS)
% title('COM XY-Position', 'FontSize', FS)
% grid on
% axis square
% xlim([-0.5 10.5])
% ylim([-5.5 5.5])
% 
% hold on
% plot(traj(:, 1), traj(:, 2), 'LineWidth', 2, 'LineStyle','--')
% lgd = legend('Actual', 'Desired')
% lgd.FontSize = FS;
end

