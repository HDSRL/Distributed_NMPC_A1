function [] = PlotEulerAngls(tout,traj,xout)

% Pitch and roll orientation
FS = 16;
angfig = figure;
subplot(3,1,1);
plot(tout,xout(:,7),'b-','LineWidth',1);
hold on
plot(tout,traj(:,7),'k-.','LineWidth',1.5);
legend('Desired','Actual', 'FontSize', FS)
ylabel('\phi', 'FontSize', FS)

subplot(3,1,2);
plot(tout,xout(:,8),'b-','LineWidth',1);
hold on
plot(tout,traj(:,8),'k-.','LineWidth',1.5);
legend('Desired','Actual', 'FontSize', FS)
ylabel('\theta', 'FontSize', FS)
xlabel('Time (sec)', 'FontSize', FS)

subplot(3,1,3);
plot(tout,xout(:,9),'b-','LineWidth',1);
hold on
plot(tout,traj(:,9),'k-.','LineWidth',1.5);
legend('Desired','Actual', 'FontSize', FS)
ylabel('\psi', 'FontSize', FS)
xlabel('Time (sec)', 'FontSize', FS)

end

