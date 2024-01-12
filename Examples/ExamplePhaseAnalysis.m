% Example Script
addpath('../')
ca = readmatrix('Calcium.csv'); %load calcium file
time = ca(:,1);%extract time as it's own variable and remove from the matrix
ca(:,1) = [];

%set presets: 
opts.thresholdsetting = 'Degree';
opts.avDegree = 4;
opts.figs = 1; %I do want figures

%identify oscillation start and end times
figure(7), plot(time, mean(calcium')) %plot calcium
[start_indx, end_indx] = identify_oscillations(ca, time, 0) %0 if you want to manually select beginning and end
xline(start_indx, 'label','This is where we start')
xline(end_indx, 'label','This is where we end')



[averagephase, sorted_highphase, finalphase]  = RunPhaseAnalysis_individual(ca, start_indx, end_indx, opts);

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

     end
