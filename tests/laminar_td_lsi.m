%{
Examine Laminar Profile of Visuomotor Activity in SEF
2022-09-02 | S P Errington

! Insert description here
%}

%% Setup: Load previously extracted data
load(fullfile(dirs.root,'results','mat_files','neuron_index.mat'))
load(fullfile(dirs.root,'results','mat_files','visual_info.mat'))
load(fullfile(dirs.root,'data','samplingBias.mat'))
load(fullfile(dirs.root,'data','td_visual_colormap.mat'))
load(fullfile(dirs.root,'results','mat_files','visual_sdf_1.mat'))
load(fullfile(dirs.root,'results','mat_files','visual_lateral.mat'))
load(fullfile(dirs.root,'results','mat_files','visual_value.mat'))

%% Organize: Find information about neurons within laminar subsample
%  Note: perpendicular sessions start at session 14; the first neuron in
%  this session is neuron 283. Thus, we take any neuron above (& including) this index
neuron_index.visual_fac_perp = neuron_index.visual_pos(neuron_index.visual_pos > 282);

% Define neurons to examine (here, we are running with facilitated, as
% there are only a few suppressive neurons in our sample)
neuronlist = []; neuronlist = neuron_index.visual_fac_perp;

clear timedepth_table
% For each neuron in our sample, get information about biophysical
% properties of the neuron (i.e. depth and spk width), and visuomotor
% properties (i.e. onset/offset). Save this information to a table for
% future use.

for neuronlist_i = 1:length(neuronlist)
    neuron_i = neuronlist(neuronlist_i);
    session = executiveBeh.neuronMatPosit(neuron_i,1);
    monkey = executiveBeh.nhpSessions.monkeyNameLabel(session);
    site = executiveBeh.bioInfo.sessionSite(session);
    spkwidth = executiveBeh.bioInfo.spkWidth(neuron_i)*25;
    depth_ch = executiveBeh.bioInfo.depthInfo(neuron_i);
    depth_um = depth_ch*150;
    mod_onset_target = visual_info.cd_detect_onset(neuron_i);
    mod_offset_target = visual_info.cd_detect_offset(neuron_i);
    mod_onset_saccade = visual_info.saccade_onset(neuron_i);
    mod_offset_saccade = visual_info.saccade_offset(neuron_i); 
    vs_flag = visual_value.value_flag_onset(visual_value.neuron_i == neuron_i);
    vs_index = visual_value.VSI_onset(visual_value.neuron_i == neuron_i);
    ls_flag = visual_lateral.lateral_flag_onset(visual_lateral.neuron_i == neuron_i); 
    ls_index = visual_lateral.LSI_onset(visual_lateral.neuron_i == neuron_i); 
    
    
    timedepth_table (neuronlist_i,:) =...
        table(neuron_i,session,monkey,site,spkwidth,depth_ch,...
        depth_um,mod_onset_target,mod_offset_target,mod_onset_saccade,mod_offset_saccade,...
        vs_flag,vs_index,ls_flag,ls_index);

end

timedepth_table(timedepth_table.ls_index > 0,:) = [];
%% Organize: Setup time-depth plot information
% Define how many depths we have within our sample (constant)
n_depths = 19;

% Define windows, and initialise an array to add to for
%   - target
td_plot_win_target = [-500:1000]; td_plot_zero_target = abs(td_plot_win_target(1));
td_plot_array_target = zeros(n_depths,length(td_plot_win_target));
%   - saccade
td_plot_win_saccade = [-1000:2000]; td_plot_zero_saccade = abs(td_plot_win_saccade(1));
td_plot_array_saccade = zeros(n_depths,length(td_plot_win_saccade));

% Then we will add to this array, marking each active period for the depth of
% each neuron with a one. This will add, so that each neuron at a depth and
% time adds a "1" value.

for neuronlist_i = 1:size(timedepth_table,1) 
    % Get the depth
    depth = timedepth_table.depth_ch(neuronlist_i);
    
    % Onset and offset times of the modulation aligned on target and
    % saccade
    onset_target = timedepth_table.mod_onset_target(neuronlist_i);
    offset_target = timedepth_table.mod_offset_target(neuronlist_i);
    onset_saccade = timedepth_table.mod_onset_saccade(neuronlist_i);
    offset_saccade = timedepth_table.mod_offset_saccade(neuronlist_i);    
    
    % But if this time goes beyond the size of our array, we cut it short.
    if offset_target > max(td_plot_win_target); offset_target = max(td_plot_win_target); end
    if offset_saccade > max(td_plot_win_saccade); offset_saccade = max(td_plot_win_saccade); end
    
    % We then add one to the array at the depth and times at which a neuron
    % is active. We do this for the target...
    td_plot_array_target(depth,onset_target+td_plot_zero_target:offset_target+td_plot_zero_target) = ...
        td_plot_array_target(depth,onset_target+td_plot_zero_target:offset_target+td_plot_zero_target) + 1;
    
    % ... and saccade, assuming we observe saccadic activity (i.e. it's not
    % a pure visual neuron).
    if ~isnan(onset_saccade)
    td_plot_array_saccade(depth,onset_saccade+td_plot_zero_saccade:offset_saccade+td_plot_zero_saccade) = ...
        td_plot_array_saccade(depth,onset_saccade+td_plot_zero_saccade:offset_saccade+td_plot_zero_saccade) + 1;

    end
end

% We then normalise the array across depths to account for the sampling
% bias. This gives us a p(active|depth) measure, rather than a
% N(active|depth).
for depth_i = 1:n_depths
    % For target:
    td_plot_array_target(depth_i,:) = td_plot_array_target(depth_i,:)./samplingBias.countPerDepth(depth_i);
    % For saccade:
    td_plot_array_saccade(depth_i,:) = td_plot_array_saccade(depth_i,:)./samplingBias.countPerDepth(depth_i);
end


%% Figure: Plot laminar profile of visuomotor activity:
% Initialise figure
timedepth_fig = figure('Renderer', 'painters', 'Position', [100 100 1000 300]);

% Plot 1: Target aligned activity ---------------------
subplot(1,3,1);hold on
% Plot the time depth activity heatmap
imagesc('XData',td_plot_win_target,'YData',1:n_depths,'CData',td_plot_array_target)
% Add markers for broad-spike neurons (those 250+ us spk width)
scatter(timedepth_table.mod_onset_target(timedepth_table.spkwidth > 249),...
    timedepth_table.depth_ch(timedepth_table.spkwidth > 249),'^','filled','MarkerEdgeColor','black','MarkerFacecolor','white')
% Add markers for narrow-spike neurons (those <250 us spk width)
scatter(timedepth_table.mod_onset_target(timedepth_table.spkwidth < 249),...
    timedepth_table.depth_ch(timedepth_table.spkwidth < 249),60,'h','filled','MarkerEdgeColor','black','MarkerFacecolor','white')
% Tidy the axes limits (as imagesc goes a little weird sometimes...)
set(gca,'XLim',[-200 500],'YLim',[0.5 n_depths+0.5],...
    'YDir','Reverse')
% Add axes labels
xlabel('Time from Target (ms)'); ylabel('Depth');
% Update the color map to match the color class of neurons (previously
% extracted)
colormap(td_visual_colormap)
% Add lines to mark onset and laminar boundaries
vline(0,'k'); hline(8.5,'k--')

% Plot 2: Saccade aligned activity ---------------------
subplot(1,3,[2 3]);hold on
% Plot the time depth activity heatmap
imagesc('XData',td_plot_win_saccade,'YData',1:n_depths,'CData',td_plot_array_saccade)
% Tidy the axes limits (as imagesc goes a little weird sometimes...)
set(gca,'XLim',[-200 1200],'YLim',[0.5 n_depths+0.5],...
    'YDir','Reverse')
% Add axes labels
xlabel('Time from Saccade (ms)'); ylabel('Depth');
% Add a colorbar to get the scale
colorbar
% Update the color map to match the color class of neurons (previously
% extracted)
colormap(td_visual_colormap)
% Add lines to mark onset (0 ms), tone onset (600 ms, fixed), and laminar boundaries
vline(0,'k'); vline(600,'k:'); hline(8.5,'k--');



