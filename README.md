# Islet_Heterogeneity
This repository contains sample code for studying islet heterogeneity. Main scripts include network/hub cell analysis, wave/phase analysis, and first responder analysis. Code is available in Matlab. All code was written by Jennifer Briggs (2024) and Richard Benninger's lab at the University of Colorado Anschutz Medical Campus.

## Examples 
Examples to run all scripts and an example calcium file are found under the folder "examples"

## Network Analysis in Matlab
To run general network analysis, you will use the function *RunNetworkAnalysis.m* 

If you want to determine the optimal threshold for getting a log-log plot, input calcium into *findoptRth*.

## Phase Analysis in Matlab
To run phase analysis, you must decide which method to use. *RunPhaseAnalysis_indiviual.m* looks at each individual oscillation and determines the phase of each cell during the depolarization stage only. *RunPhaseAnalysis_allsecondphase.m* looks at the phase cells over the entire second phase. 

To run *RunPhaseAnalysis_indiviual.m*, you must input start and end times of the depolarization of each oscillation. Use *identify_oscillations.m* to create these arrays (see examples/*ExamplePhaseAnalysis.m* for help). 

## First Responder Analysis 
First responder analysis is very simple - you are just finding when cells reach half of their max-height during the first phase of insulin/calcium secretion. The cells that reach their half-max first are the first responders. To run first responder analysis correctly, you want to make sure that only the first phase is fed into the function. Then just run *RunFirstResponder.m*. For an example, see *examples/ExampleFirstResponder.m*

