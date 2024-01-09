% Example Script
addpath('../')
ca = readmatrix('Calcium.csv'); %load calcium file

%set presets: 
opts.thresholdsetting = 'Degree';
opts.avDegree = 4;
opts.figs = 1; %I do want figures

% I want to define hubs as top 10% of high degree cells
opts.hubs = 'percentile'
opts.hubspercentile = .1;

out = RunNetworkAnalysis(ca, opts)
