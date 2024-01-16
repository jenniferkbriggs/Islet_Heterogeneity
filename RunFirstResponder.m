function [time_sorted, cells_sorted]  = RunFirstResponder(calcium)
%% by Jennifer Briggs 02/28/2022
% This function will run phase analysis given a calcium array (time x
% calcium for each cell). For this version of phase analysis, we only look
% at the entire second phase of the calcium oscillation.

% The final result will give an array with the phase lag for each cell at
% each oscillation

% INPUTS: 
% -- ca: an array containing the calcium timeseries for the islet.
%               With the rows corresponding to individual time points and the columns


% OUTPUTS:
% -- time_half: indexed time that each cell reaches half max
% -- cells_sorted: cells sorted from first to last responders

    calciumT = normalize(calcium, "range"); %normalize calcium signals]

    %sort through signals and calculate when the cell first reached half
    %max.
   for j = 1:size(calcium, 2)
       timehalf_init = find(calciumT(:,j)>0.5);
       timehalf(j) = timehalf_init(1);
   end
  
    %sort
    [time_sorted, cells_sorted] = sort(timehalf); 
end

