% Example Script
addpath('../')
ca = readmatrix('Calcium.csv'); %load calcium file
time = ca(:,2);%extract time as it's own variable and remove from the matrix
ca(:,1:2) = [];


%identify oscillation start of second phase
figure, plot(time, mean(ca,2))
title('Select beginning of second phase. Press enter when complete')
starttime =  ginput()  %here you put the time in seconds that you want to start the analysis
starttime = starttime(:,1);

%only use second phase
ca = ca(starttime:end,:);
time = time(starttime:end);

[averagephase, sorted_highphase]  = RunPhaseAnalysis_allsecondphase(ca);

% plot: 
     figure,
     %find top 3 high phase cells from first oscillation- 
     cell_high= sorted_highphase(1,1:3);
     cell_low= sorted_highphase(1,end-3:end);

     plot(time, ca, 'color',[0.9,0.9,0.9])
     hold on, line1 = plot(time, ca(:,cell_high), 'linewidth',1, 'color', 'blue')
     hold on, line2= plot(time, ca(:,cell_low), 'linewidth',1, 'color', 'red')

   
    legend([line1(1), line2(1)], {'High Phase','Low Phase'})
    set(gcf, 'color','white')
    xlabel('Time (s)')
    ylabel('Normalized Ca^{2+} Fluoresence')
    set(gca, 'box','off')

     
