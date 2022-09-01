%{
Value context encoding by of visually responsive neurons.
2022-08-27 | S P Errington

In this script, we will X, Y, and Z

Dependencies: identify_visualactivity (neuron_index and visual_info
variables), get_matchedtrials (trial_index_matched).
%}
%% Organize: Load in processed data
load(fullfile(dirs.root,'results','mat_files','neuron_index.mat'))
load(fullfile(dirs.root,'results','mat_files','visual_info.mat'))
load(fullfile(dirs.root,'results','mat_files','trial_index_matched.mat'))

%% Extract: Get average spiking activity around target onset, split by value.

neuronlist = []; neuronlist = neuron_index.visual_pos; n_neurons = length(neuronlist);
clear value_*

parfor neuronlist_i = 1:n_neurons
    neuron_i = neuronlist(neuronlist_i)
    fprintf('Extracting target aligned SDF for neuron %i of %i. \n',neuronlist_i, length(neuronlist))
    
    session_i = executiveBeh.neuronMatPosit(neuron_i,1);
    
    spk_in = load(fullfile(dirs.data_spike,['SDF_' int2str(neuron_i)]));
    
    % Define trials for comparison:
    low_value_trl = []; high_value_trl = [];
    low_value_trl = trial_index_matched.value{session_i}.low
    high_value_trl = trial_index_matched.value{session_i}.high
    
    % Define epochs of comparison
    baseline_win = []; baseline_win = timewins.baseline_target + timewins.zero;
    target_win = []; target_win = timewins.visual_target + timewins.zero;
    onset_win = []; onset_win = [0:200] + visual_info.cd_detect_onset(neuron_i) + timewins.zero;
    
    
    % Save averaged SDF for future use:
    % (1) Raw
    value_sdf_low_raw(neuronlist_i,:) = nanmean(spk_in.SDF.target(low_value_trl,:));
    value_sdf_high_raw(neuronlist_i,:) = nanmean(spk_in.SDF.target(high_value_trl,:));
    
    % (2) Z-scored
    fr_visual_baseline_mean_low = nanmean(nanmean(spk_in.SDF.target(low_value_trl,baseline_win)));
    fr_visual_baseline_mean_high = nanmean(nanmean(spk_in.SDF.target(high_value_trl,baseline_win)));
    fr_visual_baseline_std_low = nanstd(nanmean(spk_in.SDF.target(low_value_trl,baseline_win)));
    fr_visual_baseline_std_high = nanstd(nanmean(spk_in.SDF.target(high_value_trl,baseline_win)));
    value_sdf_low_zscore(neuronlist_i,:) = (nanmean(spk_in.SDF.target(low_value_trl,:))-fr_visual_baseline_mean_low)./fr_visual_baseline_std_low;
    value_sdf_high_zscore(neuronlist_i,:) = (nanmean(spk_in.SDF.target(high_value_trl,:))-fr_visual_baseline_mean_high)./fr_visual_baseline_std_high;   
    
    
    % Get mean FR between trial types for future use.
    %   for a time period, relative to target onset
    value_fr_low_target{neuronlist_i} = nanmean(spk_in.SDF.target(low_value_trl,target_win),2)
    value_fr_high_target{neuronlist_i} = nanmean(spk_in.SDF.target(high_value_trl,target_win),2)
    %   for a time period, relative to modulation onset
    value_fr_low_onset{neuronlist_i} = nanmean(spk_in.SDF.target(low_value_trl,onset_win),2)
    value_fr_high_onset{neuronlist_i} = nanmean(spk_in.SDF.target(high_value_trl,onset_win),2)    
    
end

%% Analyse: calculate the value-sensitivity for these neurons

clear visual_value
for neuronlist_i = 1:n_neurons
    
    neuron_i = neuronlist(neuronlist_i);
    session_i = executiveBeh.neuronMatPosit(neuron_i);
    monkey = executiveBeh.nhpSessions.monkeyNameLabel(session_i);
    
    [value_flag_target,value_p_target,~,value_stat_target] =...
        ttest2(value_fr_low_target{neuronlist_i},value_fr_high_target{neuronlist_i});
    [value_flag_onset,value_p_onset,~,value_stat_onset] =...
        ttest2(value_fr_low_onset{neuronlist_i},value_fr_high_onset{neuronlist_i});
    
    % Values of 1 = greater activity on high reward trials
    VSI_target = (mean(value_fr_high_target{neuronlist_i}) - mean(value_fr_low_target{neuronlist_i}))./...
        (mean(value_fr_high_target{neuronlist_i}) + mean(value_fr_low_target{neuronlist_i}));
    VSI_onset = (mean(value_fr_high_onset{neuronlist_i}) - mean(value_fr_low_onset{neuronlist_i}))./...
        (mean(value_fr_high_onset{neuronlist_i}) + mean(value_fr_low_onset{neuronlist_i}));
    
    visual_value(neuronlist_i,:) = ...
        table(neuron_i,session_i,monkey,value_flag_target,value_p_target,value_stat_target,...
        value_flag_onset,value_p_onset,value_stat_onset,VSI_target,VSI_onset);
    
end


sum(visual_value.value_flag_onset == 1 & visual_value.VSI_onset > 0)
sum(visual_value.value_flag_onset == 1 & visual_value.VSI_onset < 0)

%% Figure: Histogram of value-sensitivity index
vsi_histogram = figure('Renderer', 'painters', 'Position', [100 100 300 300]); hold on
histogram(visual_value.VSI_onset(visual_value.value_flag_onset == 0),-0.6:0.05:0.6,'LineStyle','None')
histogram(visual_value.VSI_onset(visual_value.value_flag_onset == 1),-0.6:0.05:0.6,'LineStyle','None')
xlabel('Value-sensitivity index'); ylabel('Frequency')


% Once we're done with a page, save it and close it.
filename = fullfile(dirs.root,'results','gen_figures','vsi_histogram_visual.pdf');
set(vsi_histogram,'PaperSize',[20 10]); %set the paper size to what you want
print(vsi_histogram,filename,'-dpdf') % then print it
close(vsi_histogram)

%% Figure: Population SDF for high/low value context
clear input_sdf_low input_sdf_high value_population_sdf
input_sdf_low = num2cell(value_sdf_low_zscore, 2);
input_sdf_high = num2cell(value_sdf_high_zscore, 2);

labels_value = [repmat({'1_Low'},length(input_sdf_low),1);repmat({'2_High'},length(input_sdf_high),1)];
labels_monkey = [repmat(visual_info.monkey_label(neuronlist),2,1)];

% Produce the figure, collapsed across all monkeys
value_population_sdf(1,1)=gramm('x',timewins.sdf,'y',[input_sdf_low;input_sdf_high],'color',labels_value);
value_population_sdf(1,1).stat_summary();
value_population_sdf(1,1).axe_property('XLim',[-200 500],'YLim',[0 10]);
value_population_sdf(1,1).set_names('x','Time from Target (ms)','y','FR (Z-score)');

% and create a subplot, split by monkey for transparency
value_population_sdf(2,1)=gramm('x',timewins.sdf,'y',[input_sdf_low;input_sdf_high],'color',labels_value);
value_population_sdf(2,1).stat_summary();
value_population_sdf(2,1).axe_property('XLim',[-200 500],'YLim',[-15 15]);
value_population_sdf(2,1).set_names('x','Time from Target (ms)','y','FR (Z-score)');
value_population_sdf(2,1).facet_grid([],labels_monkey); % <- performs monkey split

value_pop_figure = figure('Renderer', 'painters', 'Position', [100 100 600 600]);
value_population_sdf.draw();

% Once we're done with a page, save it and close it.
filename = fullfile(dirs.root,'results','gen_figures','valuepop_figure_visual.pdf');
set(value_pop_figure,'PaperSize',[20 10]); %set the paper size to what you want
print(value_pop_figure,filename,'-dpdf') % then print it
close(value_pop_figure)

%% Organize: Save visual value data for future use.
save(fullfile(dirs.root,'results','mat_files','visual_value.mat'),'visual_value')


