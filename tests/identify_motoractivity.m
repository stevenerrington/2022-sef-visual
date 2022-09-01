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
    fprintf('Extracting saccade aligned SDF for neuron %i of %i. \n',neuron_i, 575)
    
    session_i = executiveBeh.neuronMatPosit(neuron_i,1);
    ssd_n = length(executiveBeh.inh_SSD{session_i}(executiveBeh.inh_SSD{session_i} > 50 & executiveBeh.inh_SSD{session_i} < 500));
    spk_in = load(fullfile(dirs.data_spike,['SDF_' int2str(neuron_i)]));
    
    % Initialise arrays for parfor loop
    sdf_target_ns_fast_ssd_raw = []; sdf_target_ns_slow_ssd_raw = [];
    sdf_target_nc_ssd_raw = []; sdf_target_canc_ssd_raw = [];
    sdf_target_ns_fast_ssd_zscore = []; sdf_target_ns_slow_ssd_zscore = [];
    sdf_target_nc_ssd_zscore = []; sdf_target_canc_ssd_zscore = [];
    sdf_saccade_ns_all_raw = []; sdf_saccade_ns_all_zscore = [];
    
    % As we are looking at canceled/non-canc trials here, we will latency
    % match.
    for ssd_i = 1:ssd_n
        
        % Define trials
        trial_ns_fast = []; trial_ns_fast = executiveBeh.ttm_c.GO_NC{session_i,ssd_i}.all;
        trial_ns_slow = []; trial_ns_slow = executiveBeh.ttm_c.GO_C{session_i,ssd_i}.all;
        trial_nc = [];      trial_nc = executiveBeh.ttm_c.NC{session_i,ssd_i}.all;
        trial_c = [];       trial_c =  executiveBeh.ttm_c.C{session_i,ssd_i}.all;
        trial_all = [];     trial_all =  [trial_ns_fast; trial_ns_slow; trial_nc; trial_c] ;
        
        % Define time windows
        baseline_win = []; baseline_win = timewins.baseline_target + timewins.zero;
        saccade_win = []; saccade_win = [-100:100] + timewins.zero;
        
        % Get baseline firing rates (across all trial types)
        fr_target_baseline_mean = nanmean(nanmean(spk_in.SDF.target(trial_all,baseline_win)));
        fr_target_baseline_std = nanstd(nanmean(spk_in.SDF.target(trial_all,baseline_win)));
        
        % Get RT for matching alignment on canceled trials
        RT = []; RT = executiveBeh.TrialEventTimes_Overall{session_i}(:,4) - ...
            executiveBeh.TrialEventTimes_Overall{session_i}(:,2);
        mean_match_RT = round(nanmean(RT(trial_ns_slow)));
        alignment_index_canc = timewins.zero + mean_match_RT + [-600:600];
        
        if fr_target_baseline_std == 0 | executiveBeh.inh_trcount_SSD{session_i}(ssd_i) < 2
            sdf_target_ns_fast_ssd_raw(ssd_i,:) = nan(1,length(timewins.sdf));
            sdf_target_ns_slow_ssd_raw(ssd_i,:) = nan(1,length(timewins.sdf));
            sdf_target_nc_ssd_raw(ssd_i,:) = nan(1,length(timewins.sdf));
            sdf_target_canc_ssd_raw(ssd_i,:) = nan(1,length(timewins.sdf));
            sdf_saccade_ns_all_raw(ssd_i,:) = nan(1,length(timewins.sdf));
            
            sdf_target_ns_fast_ssd_zscore(ssd_i,:) = nan(1,length(timewins.sdf));
            sdf_target_ns_slow_ssd_zscore(ssd_i,:) = nan(1,length(timewins.sdf));
            sdf_target_nc_ssd_zscore(ssd_i,:) = nan(1,length(timewins.sdf));
            sdf_target_canc_ssd_zscore(ssd_i,:) = nan(1,length(timewins.sdf));         
        else
            sdf_target_ns_fast_ssd_raw(ssd_i,:) = nanmean(spk_in.SDF.target(trial_ns_fast,:));
            sdf_target_ns_slow_ssd_raw(ssd_i,:) = nanmean(spk_in.SDF.target(trial_ns_slow,:));
            sdf_target_nc_ssd_raw(ssd_i,:) = nanmean(spk_in.SDF.target(trial_nc,:));
            sdf_target_canc_ssd_raw(ssd_i,:) = nanmean(spk_in.SDF.target(trial_c,:));
            
            sdf_target_ns_fast_ssd_zscore(ssd_i,:) = (sdf_target_ns_fast_ssd_raw(ssd_i,:) - fr_target_baseline_mean)./...
                fr_target_baseline_std;
            sdf_target_ns_slow_ssd_zscore(ssd_i,:) = (sdf_target_ns_slow_ssd_raw(ssd_i,:) - fr_target_baseline_mean)./...
                fr_target_baseline_std;
            sdf_target_nc_ssd_zscore(ssd_i,:) = (sdf_target_nc_ssd_raw(ssd_i,:) - fr_target_baseline_mean)./...
                fr_target_baseline_std;
            sdf_target_canc_ssd_zscore(ssd_i,:) = (sdf_target_canc_ssd_raw(ssd_i,:) - fr_target_baseline_mean)./...
                fr_target_baseline_std;
            
            
            
        end
    end
    
    sdf_target_ns_raw(neuron_i,:) = nanmean(spk_in.SDF.target(executiveBeh.ttx.GO{session_i},:));
    sdf_saccade_ns_raw(neuron_i,:) = nanmean(spk_in.SDF.saccade(executiveBeh.ttx.GO{session_i},:));
    
    sdf_target_ns_zscore(neuron_i,:) = (sdf_target_ns_raw(neuron_i,:) - fr_target_baseline_mean)./...
        fr_target_baseline_std;
    sdf_saccade_ns_zscore(neuron_i,:) = (sdf_saccade_ns_raw(neuron_i,:) - fr_target_baseline_mean)./...
        fr_target_baseline_std;
    
    sdf_target_ns_fast_raw(neuron_i,:) = nanmean(sdf_target_ns_fast_ssd_raw);
    sdf_target_ns_slow_raw(neuron_i,:) = nanmean(sdf_target_ns_slow_ssd_raw);
    sdf_target_nc_raw(neuron_i,:)      = nanmean(sdf_target_nc_ssd_raw);
    sdf_target_canc_raw(neuron_i,:)    = nanmean(sdf_target_canc_ssd_raw);
    
    sdf_target_ns_fast_zscore(neuron_i,:) = nanmean(sdf_target_ns_fast_ssd_zscore);
    sdf_target_ns_slow_zscore(neuron_i,:) = nanmean(sdf_target_ns_slow_ssd_zscore);
    sdf_target_nc_zscore(neuron_i,:)     = nanmean(sdf_target_nc_ssd_zscore);
    sdf_target_canc_zscore(neuron_i,:)   = nanmean(sdf_target_canc_ssd_zscore);
end




%% Figure: Population spike-density function for visual neurons
clear input_sdf* population_sdf
input_sdf_target = num2cell(sdf_target_ns_zscore(neuron_index.visual_pos,:), 2);
input_sdf_saccade = num2cell(sdf_saccade_ns_zscore(neuron_index.visual_pos,:), 2);

% Produce the figure, collapsed across all monkeys
population_sdf(1,1)=gramm('x',timewins.sdf,'y',input_sdf_target);
population_sdf(1,2)=gramm('x',timewins.sdf,'y',input_sdf_saccade);
population_sdf(1,1).stat_summary(); population_sdf(1,2).stat_summary();
population_sdf(1,1).axe_property('XLim',[-200 800],'YLim',[-2 15]);
population_sdf(1,2).axe_property('XLim',[-400 400],'YLim',[-2 15]);
population_sdf(1,1).set_names('x','Time from Target (ms)','y','FR (Z-score)');
population_sdf(1,2).set_names('x','Time from Saccade (ms)','y','FR (Z-score)');

pop_figure = figure('Renderer', 'painters', 'Position', [100 100 800 300]);
population_sdf.draw();

% Once we're done with a page, save it and close it.
filename = fullfile(dirs.root,'results','gen_figures','population_figure_movement.pdf');
set(pop_figure,'PaperSize',[20 10]); %set the paper size to what you want
print(pop_figure,filename,'-dpdf') % then print it
close(pop_figure)




%% Analysis: Calculate visuomotor index
neuronlist = []; neuronlist = neuron_index.visual_pos; n_neurons = length(neuronlist);
clear visuomotor_*

baseline_win = []; baseline_win = timewins.baseline_target + timewins.zero;
target_win = []; target_win = timewins.visual_target + timewins.zero;
saccade_win = []; saccade_win = [-100:100] + timewins.zero;

    
for neuronlist_i = 1:n_neurons
    neuron_i = neuronlist(neuronlist_i);
    fprintf('Extracting target aligned SDF for neuron %i of %i. \n',neuronlist_i, length(neuronlist))
    
    session_i = executiveBeh.neuronMatPosit(neuron_i,1);

    % Target aligned activity
    vmi_activity_target = nanmean(sdf_target_ns_raw(neuron_i,target_win));
    % Target aligned activity
    vmi_activity_saccade = nanmean(sdf_saccade_ns_raw(neuron_i,saccade_win));
    
    vmi_index(neuronlist_i) = (vmi_activity_saccade-vmi_activity_target)/(vmi_activity_saccade+vmi_activity_target);
end


visuomotor_histogram = figure('Renderer', 'painters', 'Position', [100 100 300 300]); hold on
histogram(vmi_index,-1:0.05:1,'LineStyle','None')
xlabel('Visuomotor index'); ylabel('Frequency')

