function fighand = CH1D_Plot_Evolution(cvecs)
T = size(cvecs,2);
N = size(cvecs,1) - 1;
x_grid = linspace(0,1,N+1);
fighand = figure();
    for n = 1:T
        plot(x_grid, cvecs(:,n), 'LineWidth', 2);
        ylim([-1.1, 1.1]); % Set fixed y-axis limits
        title(['Time-step: ', num2str(n)]);
        xlabel('Spatial domain');
        ylabel('Concentration');
        grid on;
        drawnow;
        pause(0.05); 
    end
end