% Example Script
addpath('../')
ca = readmatrix('Calcium.csv'); %load calcium file
time = ca(:,2);%extract time as it's own variable and remove from the matrix
ca(:,1:2) = [];


%identify oscillation end of first phase
figure, plot(time, mean(ca,2))
title('Select end of first phase. Press enter when complete')
starttime =  ginput()  %here you put the time in seconds that you want to start the analysis
starttime = starttime(:,1);

%only use second phase
ca = ca(1:starttime,:);
time = time(1:starttime);

[time_half, cells_sorted]  = RunFirstResponder(ca);

% plot: 
     figure,
     %find top 3 high phase cells from first oscillation- 
   

     plot(time, ca, 'color',[0.9,0.9,0.9])
     hold on, line1 = plot(time, ca(:,cells_sorted(1:3)), 'linewidth',1, 'color', 'blue')
     hold on, line2= plot(time, ca(:,cells_sorted(end-3:end)), 'linewidth',1, 'color', 'red')

   
    legend([line1(1), line2(1)], {'First Responder','Last Responder'})
    set(gcf, 'color','white')
    xlabel('Time (s)')
    ylabel('Normalized Ca^{2+} Fluoresence')
    set(gca, 'box','off')

     
