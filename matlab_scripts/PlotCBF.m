function PlotCBF(tout, cbf_)
FS = 16;

% cbf_(find(cbf_<0)) = 0; 
cbf_ = cbf_-min(cbf_)-0.3;
cbf_(find(cbf_<0)) = 0;

figure;
plot(tout, cbf_', 'LineWidth', 2.5)
xlabel('Time', 'FontSize', FS)
ylabel('CBF \Psi_2(x_k)', 'FontSize', FS)
ax = gca;
ax.FontSize = FS;
end