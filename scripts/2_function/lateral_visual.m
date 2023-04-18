%{
Laterality encoding by of visually responsive neurons.
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
clear lateral_*

parfor neuronlist_i = 1:n_neurons
    neuron_i = neuronlist(neuronlist_i)
    fprintf('Extracting target aligned SDF for neuron %i of %i. \n',neuronlist_i, length(neuronlist))
    
    session_i = executiveBeh.neuronMatPosit(neuron_i,1);
    
    spk_in = load(fullfile(dirs.data_spike,['SDF_' int2str(neuron_i)]));
    
    % Define trials for comparison:
    left_target_trl = []; right_target_trl = [];
    left_target_trl = trial_index_matched.lateral{session_i}.left
    right_target_trl = trial_index_matched.lateral{session_i}.right
    
    % Define epochs of comparison
    baseline_win = []; baseline_win = timewins.baseline_target + timewins.zero;
    target_win = []; target_win = timewins.visual_target + timewins.zero;
    onset_win = []; onset_win = [0:200] + visual_info.cd_detect_onset(neuron_i) + timewins.zero;
    
    
    % Save averaged SDF for future use:
    % (1) Raw
    lateral_sdf_left_raw(neuronlist_i,:) = nanmean(spk_in.SDF.target(left_target_trl,:));
    lateral_sdf_right_raw(neuronlist_i,:) = nanmean(spk_in.SDF.target(right_target_trl,:));
    
    % (2) Z-scored
    fr_visual_baseline_mean = nanmean(nanmean(spk_in.SDF.target([left_target_trl;right_target_trl],baseline_win)));
    fr_visual_baseline_std = nanstd(nanmean(spk_in.SDF.target([left_target_trl;right_target_trl],baseline_win)));
    lateral_sdf_left_zscore(neuronlist_i,:) = (nanmean(spk_in.SDF.target(left_target_trl,:))-fr_visual_baseline_mean)./fr_visual_baseline_std;
    lateral_sdf_right_zscore(neuronlist_i,:) = (nanmean(spk_in.SDF.target(right_target_trl,:))-fr_visual_baseline_mean)./fr_visual_baseline_std;   
    
    lateral_sdf_left_zscore_saccade(neuronlist_i,:) = (nanmean(spk_in.SDF.saccade(left_target_trl,:))-fr_visual_baseline_mean)./fr_visual_baseline_std;
    lateral_sdf_right_zscore_saccade(neuronlist_i,:) = (nanmean(spk_in.SDF.saccade(right_target_trl,:))-fr_visual_baseline_mean)./fr_visual_baseline_std;   
       
    
    % Get mean FR between trial types for future use.
    %   for a time period, relative to target onset
    lateral_fr_left_target{neuronlist_i} = nanmean(spk_in.SDF.target(left_target_trl,target_win),2)
    lateral_fr_right_target{neuronlist_i} = nanmean(spk_in.SDF.target(right_target_trl,target_win),2)
    %   for a time period, relative to modulation onset
    lateral_fr_left_onset{neuronlist_i} = nanmean(spk_in.SDF.target(left_target_trl,onset_win),2)
    lateral_fr_right_onset{neuronlist_i} = nanmean(spk_in.SDF.target(right_target_trl,onset_win),2)    
    
end

%% Analyse: calculate the laterality-sensitivity for these neurons

clear visual_lateral
for neuronlist_i = 1:n_neurons
    
    neuron_i = neuronlist(neuronlist_i);
    session_i = executiveBeh.neuronMatPosit(neuron_i);
    monkey = executiveBeh.nhpSessions.monkeyNameLabel(session_i);
    
    [lateral_flag_target,lateral_p_target,~,lateral_stat_target] =...
        ttest2(lateral_fr_left_target{neuronlist_i},lateral_fr_right_target{neuronlist_i});
    [lateral_flag_onset,lateral_p_onset,~,lateral_stat_onset] =...
        ttest2(lateral_fr_left_onset{neuronlist_i},lateral_fr_right_onset{neuronlist_i});
    
    % Values of 1 = greater activity on right side trials
    LSI_target = (mean(lateral_fr_right_target{neuronlist_i}) - mean(lateral_fr_left_target{neuronlist_i}))./...
        (mean(lateral_fr_right_target{neuronlist_i}) + mean(lateral_fr_left_target{neuronlist_i}));
    LSI_onset = (mean(lateral_fr_right_onset{neuronlist_i}) - mean(lateral_fr_left_onset{neuronlist_i}))./...
        (mean(lateral_fr_right_onset{neuronlist_i}) + mean(lateral_fr_left_onset{neuronlist_i}));
    
    visual_lateral(neuronlist_i,:) = ...
        table(neuron_i,session_i,monkey,lateral_flag_target,lateral_p_target,lateral_stat_target,...
        lateral_flag_onset,lateral_p_onset,lateral_stat_onset,LSI_target,LSI_onset);
    
end

sum(visual_lateral.lateral_flag_onset == 1 & visual_lateral.LSI_onset > 0)
sum(visual_lateral.lateral_flag_onset == 1 & visual_lateral.LSI_onset < 0)

%% Figure: Histogram of value-sensitivity index
laterality_histogram = figure('Renderer', 'painters', 'Position', [100 100 300 300]); hold on
histogram(visual_lateral.LSI_onset(visual_lateral.lateral_flag_onset == 0),-1:0.05:1,'LineStyle','None')
histogram(visual_lateral.LSI_onset(visual_lateral.lateral_flag_onset == 1),-1:0.05:1,'LineStyle','None')
xlabel('Laterality-sensitivity index'); ylabel('Frequency')


% Once we're done with a page, save it and close it.
filename = fullfile(dirs.root,'results','gen_figures','lsi_histogram_lateral.pdf');
set(laterality_histogram,'PaperSize',[20 10]); %set the paper size to what you want
print(laterality_histogram,filename,'-dpdf') % then print it
close(laterality_histogram)

ttest(visual_lateral.LSI_onset)
mean(visual_lateral.LSI_onset)

%% Figure: Population SDF for right/left value context
clear input_sdf_left input_sdf_right lateral_population_sdf
input_sdf_left = num2cell(lateral_sdf_left_zscore, 2);
input_sdf_right = num2cell(lateral_sdf_right_zscore, 2);

labels_value = [repmat({'1_left'},length(input_sdf_left),1);repmat({'2_right'},length(input_sdf_right),1)];
labels_monkey = [repmat(visual_info.monkey_label(neuronlist),2,1)];

% Produce the figure, collapsed across all monkeys
lateral_population_sdf(1,1)=gramm('x',timewins.sdf,'y',[input_sdf_left;input_sdf_right],'color',labels_value);
lateral_population_sdf(1,1).stat_summary();
lateral_population_sdf(1,1).axe_property('XLim',[-200 500],'YLim',[0 10]);
lateral_population_sdf(1,1).set_names('x','Time from Target (ms)','y','FR (Z-score)');

% and create a subplot, split by monkey for transparency
lateral_population_sdf(2,1)=gramm('x',timewins.sdf,'y',[input_sdf_left;input_sdf_right],'color',labels_value);
lateral_population_sdf(2,1).stat_summary();
lateral_population_sdf(2,1).axe_property('XLim',[-200 500],'YLim',[-15 15]);
lateral_population_sdf(2,1).set_names('x','Time from Target (ms)','y','FR (Z-score)');
lateral_population_sdf(2,1).facet_grid([],labels_monkey); % <- performs monkey split

lateral_pop_figure = figure('Renderer', 'painters', 'Position', [100 100 600 600]);
lateral_population_sdf.draw();

% Once we're done with a page, save it and close it.
filename = fullfile(dirs.root,'results','gen_figures','lateralpop_figure_visual.pdf');
set(lateral_pop_figure,'PaperSize',[20 10]); %set the paper size to what you want
print(lateral_pop_figure,filename,'-dpdf') % then print it
close(lateral_pop_figure)


%%




%% Organize: Save output
save(fullfile(dirs.root,'results','mat_files','visual_lateral.mat'),'visual_lateral')