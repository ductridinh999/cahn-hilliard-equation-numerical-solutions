N = 128; 
T = 2000; 
ep = 0.01; 
mu = 1; 
seed = 10; 
k = 2; 
[cvecs, time] = CH1D12(N,T,ep,mu,seed,k);

v = VideoWriter('CH_FDM_1D.mp4', 'MPEG-4');
open(v);

figure();
y_limits = [min(cvecs(:)) max(cvecs(:))]; 

for i = 1:size(cvecs, 2) 
    plot(cvecs(:, i), 'LineWidth', 2);
    ylim(y_limits);      
    xlim([1 N]);         
    title(['Step: ' num2str(i)]);
    grid on;
    
    frame = getframe(gcf);
    writeVideo(v, frame);
end

close(v);
disp('Video saved.');