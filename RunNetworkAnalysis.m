function out = RunNetworkAnalysis(ca, opts)
% Run Network analysis: 
% Jennifer Briggs 2022
% This function will run network analysis given a calcium array (time x calcium for each cell). 
% The final result will give a list of hub cells, the correlation coefficient matrix, and the adjacency matrix in a matlab structure. 
% INPUTS: 
% -- ca: an array containing the calcium timeseries for the islet.
%               With the rows corresponding to individual time points and the columns
%               corresponding to different cells (see example if needed)
% -- opts:   optional options structure. 
%      opts.thresholdsetting: Structure to define how the correlation threshold is set
%           there are a few options (see findoptRth.m).  If nothing is given here, 
%           the default is to set a threshold of 0.9
%           1) You can preset a threshold, in which case set opts.thresholdsetting = 'preset' 
%               and opts.Threshold = x. Where x is a float with the value of the threshold
%           2) You can define the average degree of the islet, in which
%               case set opts.thresholdsetting = 'Degree' and opts.avDegree =  x.
%               Where x is an int with the average degree. 
%           3) You can find a threshold such that the degree distribution
%               roughly matches a scale free distribution (See Stozer2013 or
%               Briggs2023). In this case, set opts.thresholdsetting =
%               'Scale-Free'
%     opts.figs: boolean - set to 1 if you want matlab to make figures or 0
%       if not. Default is 1
%     opts.hubs: Set how hubs are defined: 
%           1) 'sixty' - 60% of max degree Johnston2016 -- THIS IS DEFAULT
%           2) 'percentile' - Top x% of the islet (e.g. top 10% of the
%           islet. If this is set, the top percentile must also be set -
%           e.g. opts.hubspercentile = .1


if ~exist('opts')
    opts = struct()
end

if ~isfield(opts, 'figs')
    Opts.figs = 1 %Set 1 if you want figures, 0 if not
else
    Opts.figs = opts.figs; %Opts is the option structure for the NetworkAnalysis script
end

%1) Find threshold
if ~isfield(opts, 'thresholdsetting')
    Threshold = 0.9 %Default threshold
else
    switch opts.thresholdsetting
        case 'Degree'
            Opts.Method = 'Degree';
            try
            Opts.avDeg = opts.avDegree;
            catch
                error('Average Degree not defined. Input opts.avDegree to set the average degree')
            end
        case 'Scale-Free'
            Opts.Method = 'Scale-Free' 
            %%Set the bounds for max and min average degree for scale free: 
            Opts.Max = 20
            Opts.Min = 2
    end

    %NOW RUN TO FIND THRESHOLD
    Threshold = findoptRth(ca, Opts);
end


%2) Run network analysis!

[degree, Adj, kpercent, histArrayPercShort,pval,Rij,s] = NetworkAnalysis(ca, Threshold, Opts);


if ~isfield(opts, 'hubs')
    opts.hubs = 'sixty'; %default is to define hubs as top 60% of maximum degree
end
switch opts.hubs
    case 'sixty'
        a = find(kpercent > 60); hubthreshold = (a(1)); %Find degree threshold for hubs 
        Hubs = find(degree>hubthreshold) %List each hub
    case 'percentile'   
        try
        [sorteddeg, sortedcell] = sort(degree, 'descend');
        Hubs = sortedcell(1:length(Adj)*opts.hubspercentile)
        catch
            error('Percentile not defined. Input opts.hubspercentile')
        end
end


%create output
out.Hubs = Hubs;
out.degree = degree;
out.Adj = Adj;
out.Rij = Rij;
out.Threshold = Threshold;