%{
Extract spiking activity and identify neurons with marked visual modulation
2022-08-26 | S P Errington

In this script, we will loop through our neurons and identify neurons that
show periods of significant deviation of activity from baseline, in
response to a visual target. We will do this through identifying neurons
that show a significant difference in mean FR between these two periods,
and finding when this activity starts/ends.

This code will produce figures of these spike density functions for manual
curation (in the results folder), and will give the index of visual neurons
(neuron_index and visual_info variables).
%}


%% Extract: Get average spiking activity around target onset.
parfor neuron_i = 1:575
    fprintf('Extracting target aligned SDF for neuron %i of %i. \n',neuron_i, 575)
    
    session_i = executiveBeh.neuronMatPosit(neuron_i,1);
    
    spk_in = load(fullfile(dirs.data_spike,['SDF_' int2str(neuron_i)]));
    
    trial_in_all = []; trial_in_all = executiveBeh.ttx.GO{session_i};
    
    foreperiod_dist = [];
    foreperiod_dist = ...
        executiveBeh.TrialEventTimes_Overall{session_i}(trial_in_all,2)-...
        executiveBeh.TrialEventTimes_Overall{session_i}(trial_in_all,1);
    
    Y = [];
    [Y,~] = discretize(foreperiod_dist,3);

    trial_fp_long = []; trial_fp_mid = []; trial_fp_short = [];
    trial_fp_long = trial_in_all(Y == 3);
    trial_fp_mid = trial_in_all(Y == 2);
    trial_fp_short = trial_in_all(Y == 1);
    
    
    baseline_win = []; baseline_win = [-250:-50] + timewins.zero;
    
    fr_visual_baseline_mean = nanmean(nanmean(spk_in.SDF.trialStart(trial_in_all,baseline_win)));
    fr_visual_baseline_std = nanstd(nanmean(spk_in.SDF.trialStart(trial_in_all,baseline_win)));
    
    
    if fr_visual_baseline_std == 0
        sdf_visual_raw_short(neuron_i,:) = nan(1,length(timewins.sdf));
        sdf_visual_zscore_short(neuron_i,:) = nan(1,length(timewins.sdf));
        
        sdf_visual_raw_mid(neuron_i,:) = nan(1,length(timewins.sdf));
        sdf_visual_zscore_mid(neuron_i,:) = nan(1,length(timewins.sdf));
        
        sdf_visual_raw_long(neuron_i,:) = nan(1,length(timewins.sdf));
        sdf_visual_zscore_long(neuron_i,:) = nan(1,length(timewins.sdf));
        
    else
        sdf_visual_raw_short(neuron_i,:) = nanmean(spk_in.SDF.trialStart(trial_fp_short,:));
        sdf_visual_zscore_short(neuron_i,:) = (sdf_visual_raw_short(neuron_i,:) - fr_visual_baseline_mean)./...
            fr_visual_baseline_std;
        
        sdf_visual_raw_mid(neuron_i,:) = nanmean(spk_in.SDF.trialStart(trial_fp_mid,:));
        sdf_visual_zscore_mid(neuron_i,:) = (sdf_visual_raw_mid(neuron_i,:) - fr_visual_baseline_mean)./...
            fr_visual_baseline_std;
        
        sdf_visual_raw_long(neuron_i,:) = nanmean(spk_in.SDF.trialStart(trial_fp_long,:));
        sdf_visual_zscore_long(neuron_i,:) = (sdf_visual_raw_long(neuron_i,:) - fr_visual_baseline_mean)./...
            fr_visual_baseline_std;
    end
end





%% Figure: Population SDF for right/left value context
clear input_sdf_short input_sdf_mid foreperiod_sdf 
input_sdf_short = num2cell(sdf_visual_zscore_short, 2);
input_sdf_mid = num2cell(sdf_visual_zscore_mid, 2);
input_sdf_long = num2cell(sdf_visual_zscore_long, 2);

labels_value = [repmat({'1_short'},length(input_sdf_short),1);repmat({'2_mid'},length(input_sdf_mid),1);repmat({'3:long'},length(input_sdf_long),1)];
labels_monkey = [repmat(visual_info.monkey_label(1:575),3,1)];

% Produce the figure, collapsed across all monkeys
foreperiod_sdf(1,1)=gramm('x',timewins.sdf,'y',[input_sdf_short;input_sdf_mid;input_sdf_long],'color',labels_value);
foreperiod_sdf(1,1).stat_summary();
foreperiod_sdf(1,1).axe_property('XLim',[-800 100],'YLim',[-1 2]);
foreperiod_sdf(1,1).set_names('x','Time from Target (ms)','y','FR (Z-score)');

% and create a subplot, split by monkey for transparency
foreperiod_sdf(2,1)=gramm('x',timewins.sdf,'y',[input_sdf_short;input_sdf_mid;input_sdf_long],'color',labels_value);
foreperiod_sdf(2,1).stat_summary();
foreperiod_sdf(2,1).axe_property('XLim',[-800 100],'YLim',[-1 2]);
foreperiod_sdf(2,1).set_names('x','Time from Target (ms)','y','FR (Z-score)');
foreperiod_sdf(2,1).facet_grid([],labels_monkey); % <- performs monkey split

foreperiod_pop_figure = figure('Renderer', 'painters', 'Position', [100 100 600 600]);
foreperiod_sdf.draw();

