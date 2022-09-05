% Get matched trials for fairer comparison

 
%% Analyse: Extract trial indices across all sessions
for session = 1:29

% Define high/low value, and right/left no-stop trials
clear GOrL GOrH GOlL GOlH
GOrL = executiveBeh.ttx.GO_Right_L{session}; % No-stop, right target, low value
GOrH = executiveBeh.ttx.GO_Right_H{session}; % No-stop, right target, high value
GOlL = executiveBeh.ttx.GO_Left_L{session};  % No-stop, left target, low value
GOlH = executiveBeh.ttx.GO_Left_H{session};  % No-stop, left target, low value

% Get response times for whole session, on which trial matching will be
% dervied from
clear RTX
RTX = executiveBeh.TrialEventTimes_Overall{session}(:,4)-...
    executiveBeh.TrialEventTimes_Overall{session}(:,2);

% Laterality matching: %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Here, there is an equal number of high and low reward trials in for
%   each direction, and RT's within condition are matched.
clear GOrL_matched GOlL_matched GOrH_matched GOlH_matched
[GOrL_matched, GOlL_matched,  ~,  ~] = rtmatchTrials_laterality([GOrL RTX(GOrL)], [GOlL RTX(GOlL)], 10);
[GOrH_matched, GOlH_matched,  ~,  ~] = rtmatchTrials_laterality([GOrH RTX(GOrH)], [GOlH RTX(GOlH)], 10);

trial_index_matched.lateral{session}.right = [GOrL_matched(:,1); GOrH_matched(:,1)];
trial_index_matched.lateral{session}.left = [GOlL_matched(:,1); GOlH_matched(:,1)];

% Value matching: %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%&&&&&
%   Here, there is an equal number of left and right direction trials in for
%   each value context, and RT's within condition are matched.

clear GOL GOH
GOL = [GOrL; GOlL]; GOH = [GOrH; GOlH];
[GOL_matched,  GOH_matched] = rtmatchTrials_value(GOL, GOH, executiveBeh.TrialEventTimes_Overall{session}, executiveBeh.SessionInfo{session});
trial_index_matched.value{session}.low = GOL_matched;
trial_index_matched.value{session}.high = GOH_matched;

end

%% Organize: output data
save(fullfile(dirs.root,'results','mat_files','trial_index_matched.mat'),'trial_index_matched')
