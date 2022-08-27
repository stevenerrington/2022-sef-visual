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
user = 'erringsp';
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

%% Analysis codes %%%%%%%%%%%%%%%%%%%%%%%%%%%
% Identify: This script will extract X, Y, and Z.
identify_visualactivity

% Functional: This script will extract X, Y, and Z.
functional_visualactivity




%%% IN DEVLEOPMENT %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Laterality: This script will extract X, Y, and Z.
get_matchedtrials;
value_visual;
lateral_visual;

laminar_visual