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
    
    sdf_target_ns_raw(neuron_i,:) = nanmean(spk_in.SDF.target(trial_all,:));
    sdf_saccade_ns_raw(neuron_i,:) = nanmean(spk_in.SDF.saccade(trial_all,:));
    
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



%% Extract: Produce summary sheet figure for mean SDF

n_plot_x = 4; n_plot_y = 3; n_plot_sheet = n_plot_x*n_plot_y;
n_batches = round(575/n_plot_sheet,-1)+1;
getColors

neuron_i = 0;
for page_i = 1:n_batches
    fig_out = figure('Renderer', 'painters', 'Position', [100 100 1200 800]);
    
    for plot_i = 1:n_plot_sheet
        neuron_i = neuron_i+1;
        
        
        try
            session_i = executiveBeh.neuronMatPosit(neuron_i,1);
            ssd_mean = executiveBeh.inh_SSD{session_i}(executiveBeh.midSSDindex(session_i));
            
            subplot(n_plot_x, n_plot_y, plot_i); hold on
            plot(timewins.sdf, sdf_target_canc_raw(neuron_i,:),'color',colors.canceled,'LineWidth',0.5)
            plot(timewins.sdf, sdf_target_ns_slow_raw(neuron_i,:),'color',colors.nostop,'LineWidth',0.5)
            
            xlim([-100 800]); vline(0,'k--'); hline(0,'k--');
            vline(ssd_mean,'r--'); vline(ssd_mean + executiveBeh.SSRT_integrationWeighted_all(session_i),'r:')
            xlabel('Time from Target (ms)')
            ylabel('Firing rate (z-score)')
            title(['Neuron: ' int2str(neuron_i) ])
            
        catch
            continue
        end
        
    end
    
    filename = fullfile(dirs.root,'results','sdf_overview_figs',['sdf_saccade_overview_pg' int2str(page_i) '.pdf']);
    set(fig_out,'PaperSize',[20 10]); %set the paper size to what you want
    print(fig_out,filename,'-dpdf') % then print it
    close(fig_out)
end

%% Analysis: Probe deviation from a baseline period during target onset period

%% Figure: Population spike-density function for visual neurons
clear input_sdf population_sdf
input_sdf_canc = num2cell(sdf_target_canc_zscore(neuron_index.visual_pos,:), 2);
input_sdf_nostop = num2cell(sdf_saccade_ns_zscore(neuron_index.visual_pos,:), 2);
class_label = [repmat({'1_Canceled'},length(input_sdf_canc),1); repmat({'2_Nostop'},length(input_sdf_nostop),1)];
% Produce the figure, collapsed across all monkeys
population_sdf(1,1)=gramm('x',timewins.sdf,'y',[input_sdf_canc;input_sdf_nostop],'color',class_label);
population_sdf(1,1).stat_summary();
population_sdf(1,1).axe_property('XLim',[-200 800],'YLim',[-2 15]);
population_sdf(1,1).set_names('x','Time from Target (ms)','y','FR (Z-score)');

pop_figure = figure('Renderer', 'painters', 'Position', [100 100 400 300]);
population_sdf.draw();

% Once we're done with a page, save it and close it.
filename = fullfile(dirs.root,'results','gen_figures','population_figure_movement.pdf');
set(pop_figure,'PaperSize',[20 10]); %set the paper size to what you want
print(pop_figure,filename,'-dpdf') % then print it
close(pop_figure)






















