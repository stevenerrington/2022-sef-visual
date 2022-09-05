%{
Properties of visually responsive neurons in the Supplementary Eye Field
2022-08-23

Dependencies:
gramm toolbox
donut function (in src folder)

%}

%% Setup workspace %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all; clc
% Get directories %%%%%%%%%%%%%%%%%%%%%%%%
user = 'erringsp_WH633'; % erringsp
dirs = setup_dirs(user);

% Define windows %%%%%%%%%%%%%%%%%%%%%%%%%%
%   Data windows
timewins.sdf = [-999:1500];
timewins.zero = abs(min(timewins.sdf));
%   Analysis epochs
timewins.baseline_target = [-250:0];
timewins.visual_target = [50:200];

% Load session information
load(fullfile(dirs.data_proc,'behavior\','executiveBeh.mat'))
load(fullfile(dirs.data_proc,'behavior\','bayesianSSRT.mat'))
load(fullfile(dirs.data_proc,'behavior\','FileNames.mat'))

% Extract matched trial indices
get_matchedtrials

%% Analysis codes %%%%%%%%%%%%%%%%%%%%%%%%%%%
%{
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Identify visual neuron activity
- This script will:
    x Print target-aligned SDF
    x Look for neurons with a significant change in modulation between
      baseline and the visual period
    x Determine the onset and offset of these neurons (relative to target)
    x Plot individual SDFs for ID'd visual neurons
    x Output the index of visual/non-visual neurons, no-stop SDFs, and
      broad information about visually responsive neurons
%}

identify_visualactivity

%{
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Explore the functional features of visual related neurons
- This script will:
    x Find neurons that are facilitated/suppressed relative to target onset
    x Find the proportion of these classes that are broad/narrow spikes
    x Plot the population average SDF for these classes
    x Plot the cumulative density function of onset times for these
      neurons.
    x Plot example neurons in each class (fac and suppressed)
    x Output the updated neuron-index with fac and sup visual neurons
%}

functional_visualactivity
%{ 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Determine whether these neurons have motor features
- This script will:
    x Generate SDFs for no-stop trials, aligned on target and saccade
    x Calculate the visuomotor index for each neuron in our sample
    x Plot the population average for facilitated neurons, aligned on
      target and saccade
    x Plot the distribution of calculated VMI
    x Find the onset/offset of activity relative to saccade
%}

identify_motoractivity

%{
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Compare the latency of our neuronal sample with those in other cortical and
 subcortical areas.
- This script will:
    x Import latencies derived from other neurons (Schmolesky et al., 1998)
    x Calculate CDF of onset latencies.
    x Plot CDF and boxplot of these latencies
%}

compare_latency;

%{
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Determine whether these neurons are moderated by value context 
- This script will:
    x Generate SDFs for no-stop trials, split by low and high reward
      contexts
    x Calculate the value-sensitivity index (VSI) for each facilitated visual
      neuron
    x Plot the population average for facilitated neurons, aligned on
      target, for high and low reward contexts
    x Plot the distribution of calculated VSI
%}

value_visual;

%{
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Determine whether these neurons are moderated by laterality of target 
- This script will:
    x Generate SDFs for no-stop trials, split by left and right target
      presentations
    x Calculate the laterality-sensitivity index (LSI) for each facilitated visual
      neuron
    x Plot the population average for facilitated neurons, aligned on
      target, for left and right target presentations
    x Plot the distribution of calculated LSI
%}

lateral_visual;

%{
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Explore the laminar distribution of visually-responsive neurons in SEF
    x ID visual facilitated neurons in perpendicular SEF penetrations
    x Plot time-depth plot of visual neurons aligned on target and saccade.
    x Look at the distriubtion of LSI through depth.
%}

laminar_visual;











