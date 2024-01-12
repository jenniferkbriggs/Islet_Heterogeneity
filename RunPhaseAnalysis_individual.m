function [averagephase, sorted_highphase, finalphase]  = RunPhaseAnalysis_individual(calcium, start_indx, end_indx, opts);
%% by Jennifer Briggs 02/28/2022
% This function will run phase analysis given a calcium array (time x
% calcium for each cell). For this version of phase analysis, we only look
% at the depolarization phase of the calcium oscillation. To do this we
% must first input start and end times (note that these times are the
% data indices, not actual times). 

% !!!! If you want to automatically locate  these, run
% "identify_oscillations.m" !!!

% The final result will give an array with the phase lag for each cell at
% each oscillation

% INPUTS: 
% -- ca: an array containing the calcium timeseries for the islet.
%               With the rows corresponding to individual time points and the columns
%               corresponding to different cells (see example if needed)
% -- start_indx: array with the index for the start time of each
% oscillation
% -- end_indx : array with the index for the peak of each oscillation
% ----:  opts.figs: boolean - set to 1 if you want matlab to make figures or 0
%       if not. Default is 1

% OUTPUTS:
% -- averagephase: average phase (number are indices not time) for each
% cell
% -- sorted_highpase: cells sorted from high to low phase for each
% oscillation
% -- final_phase: phase value for each cell for each oscillation. 


%% Wave origin analysis
    
%in order to make the analysis more accurate, we linearly interpolate between the
%points to artificially increase resolution. We may need to play with this
%as the outputs are not getting very good differentiation between phases.

%Demean the calcium  (this isn't actually necessary becasue the cross
%correlation calculates demeaned but it can be helpful for visualization)
    calstore = calcium;
    calcium_demeaned = (calcium-min(calcium))./(max(calcium)-min(calcium)); 
    %Lets normalize all calcium ranges, assuming the average flouresence value for a cell is a
    %reflection on the staining rather than the actual cell's properties
    cashort = [];

        
          
   for i = 1:length(end_indx)
    % run phase analysis for each oscillation
    cashort =  [calcium_demeaned(start_indx(i):end_indx(i),:)];
    step = .005; %this gives how much to interpolate by
    xq = 1:step:size(cashort,1)+1;
    vq1 = interp1(1:size(cashort,1),cashort,xq); %may need to investigate a better way to do this
    %vq2 = spline(1:size(cashort,1),cashort',xq);
    
    numcells=size(cashort,2);

    calciumT = (vq1);                           % new calcium time course
    [row,col] = find(isnan(calciumT)); %remove NaN's
    calciumT = calciumT(1:row(1)-1,:);
   % 
   % 
   %  % 2. MAKING THE REFERENCE SIGNAL TO COMPARE THE SIGNAL OF INDIVIDUAL CELL'S CROSSCORRELATION WITH THIS REFERENCE
    clear vq1 cashort xq camax

    % 3. OBTAINING CROSS-CORRELATION OF THE REFERENCE SIGNAL (MEANISLET) WITH EACH INDIVIDUAL CELL
    tic

    %   HERE WE CALCULATE PHASE USING CROSS CORRELATIONS
        MeanIslet= mean(calciumT,2);       % reference signal. Index (i-(st-1)) is here to account for times when st is not 0, otherwise indexing is wrong
        st = size(calciumT,1);
        for j=1:numcells % itterative index for cells
           [c(:,j)]=xcov(calciumT(1:round(4*st/5),j),MeanIslet,'none');      % cross-covariance  measures the similarity between currentcell and shifted (lagged) copies of MeanIslet as a function of the lag      % cross-covariance  measures the similarity between currentcell and shifted (lagged) copies of MeanIslet as a function of the lag.
        end
        toc
    
        [maxCV, maxCL]=max(c);
        while length(unique(maxCL))<2;
            c(unique(maxCL), :) = [];
            [maxCV, maxCL]=max(c);
        end
        clear c

        %maxCL is the max cross correlation - here we just demean
    newmaxCLvec = maxCL-mean(maxCL);

   %  clear maxCL
    if isempty(nonzeros(newmaxCLvec))
        keyboard
        error('Something wrong with calcium signals')
    end

    %this is where you get the final output. phasevecsort gives you the
    %sorted vector of the phase lag compared to the islet mean.
    %cells_sorted is what you are really interested in. This gives the cell
    %index in order from high phase (start earlier) to low
    %phase (start later)
    [phasevecsort_init, cells_sortedinit] = sort(newmaxCLvec); 
    sorted_highphase(i,:)= cells_sortedinit;
    phasevecsort(i,:) = phasevecsort_init;
    finalphase(i,:) = newmaxCLvec;
    
   end
    %% find the average phase for all oscillations
    averagephase = mean(finalphase);

end

