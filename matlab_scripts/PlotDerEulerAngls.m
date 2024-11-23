function [] = PlotDerEulerAngls(tout,traj,xout)

% Pitch and roll orientation
angfig = figure;
subplot(3,1,1);
plot(tout,xout(:,10),'b-','LineWidth',1);
hold on
plot(tout,traj(:,10),'k-.','LineWidth',1.5);
legend('Desired','Actual')
ylabel('d\phi')

subplot(3,1,2);
plot(tout,xout(:,11),'b-','LineWidth',1);
hold on
plot(tout,traj(:,11),'k-.','LineWidth',1.5);
legend('Desired','Actual')
ylabel('d\theta')
xlabel('Time (sec)')

subplot(3,1,3);
plot(tout,xout(:,12),'b-','LineWidth',1);
hold on
plot(tout,traj(:,12),'k-.','LineWidth',1.5);
legend('Desired','Actual')
ylabel('d\psi')
xlabel('Time (sec)')

end

