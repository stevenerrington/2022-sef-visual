%{
Identify visuomotor activity
2022-08-31 | S P Errington
%}

%% Extract: Get average spiking activity around target onset.

baseline_win = []; baseline_win = timewins.baseline_target + timewins.zero;
target_win = []; target_win = timewins.visual_target + timewins.zero;
saccade_win = []; saccade_win = [-100:100] + timewins.zero;

clear vmi*
parfor neuron_i = 1:575
    fprintf('Extracting saccade aligned SDF for neuron %i of %i. \n',neuron_i, 575)
    
    session_i = executiveBeh.neuronMatPosit(neuron_i,1);
    spk_in = load(fullfile(dirs.data_spike,['SDF_' int2str(neuron_i)]));
    
    fr_target_baseline_mean = nanmean(nanmean(spk_in.SDF.target(executiveBeh.ttx.GO{session_i},baseline_win)));
    fr_target_baseline_std = nanstd(nanmean(spk_in.SDF.target(executiveBeh.ttx.GO{session_i},baseline_win)));
    
    sdf_target_ns_raw(neuron_i,:) = nanmean(spk_in.SDF.target(executiveBeh.ttx.GO{session_i},:));
    sdf_saccade_ns_raw(neuron_i,:) = nanmean(spk_in.SDF.saccade(executiveBeh.ttx.GO{session_i},:));
    sdf_tone_ns_raw(neuron_i,:) = nanmean(spk_in.SDF.tone(executiveBeh.ttx.GO{session_i},:));
    
    sdf_target_ns_zscore(neuron_i,:) = (sdf_target_ns_raw(neuron_i,:) - fr_target_baseline_mean)./...
        fr_target_baseline_std;
    sdf_saccade_ns_zscore(neuron_i,:) = (sdf_saccade_ns_raw(neuron_i,:) - fr_target_baseline_mean)./...
        fr_target_baseline_std;
    sdf_tone_ns_zscore(neuron_i,:) = (sdf_tone_ns_raw(neuron_i,:) - fr_target_baseline_mean)./...
        fr_target_baseline_std;
    
    % Target aligned activity
    vmi_activity_target = []; vmi_activity_saccade = [];
    vmi_activity_target = nanmean(spk_in.SDF.target(executiveBeh.ttx.GO{session_i},target_win),2);
    % Saccade aligned activity
    vmi_activity_saccade = nanmean(spk_in.SDF.saccade(executiveBeh.ttx.GO{session_i},saccade_win),2);
    
    [vmi_flag(neuron_i,1),vmi_p(neuron_i,1),~,vmi_stat{neuron_i,1}] =...
        ttest2(vmi_activity_target,vmi_activity_saccade);
    
    vmi_dir(neuron_i,1) = nanmean(vmi_activity_target) > nanmean(vmi_activity_saccade);
    
    vmi_index(neuron_i,1) = (nanmean(vmi_activity_saccade)-nanmean(vmi_activity_target))./...
        (nanmean(vmi_activity_saccade)+nanmean(vmi_activity_target));
end

%% Figure: Population spike-density function for visual neurons
clear input_sdf* population_sdf
input_sdf_target = num2cell(sdf_target_ns_zscore(neuron_index.visual_neg,:), 2);
input_sdf_saccade = num2cell(sdf_saccade_ns_zscore(neuron_index.visual_neg,:), 2);

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

%% Figure: Plot VMI
visuomotor_histogram = figure('Renderer', 'painters', 'Position', [100 100 300 300]); hold on
histogram(vmi_index(neuron_index.visual_pos),-1:0.05:1,'LineStyle','None')
xlabel('Visuomotor index'); ylabel('Frequency')
vline([-0.4 0.4],'k:'); vline(0,'k')
filename = fullfile(dirs.root,'results','gen_figures','vmi_histogram_visualneuron_pos.pdf');
set(visuomotor_histogram,'PaperSize',[20 10]); %set the paper size to what you want
print(visuomotor_histogram,filename,'-dpdf') % then print it
close(visuomotor_histogram)


%% Analysis: Probe deviation from a baseline period during target onset period
cd_1 = 3; cd_1t = 75; % Signal must be above 3 sd for 100 ms,
cd_2 = 6; cd_2t = 25;  % ...and above 6 sd for 25 ms.

clear onset_cd_thresh offset_cd_thresh
% Get the neuron and z-scored activity
for neuron_i = 1:575
    neuron_sdf = []; neuron_sdf = sdf_saccade_ns_zscore(neuron_i,:);
    
    % Find when all points at which activity rises above thresholds
    clear start_cd1 start_cd2 len_cd1 len_cd2
    [start_cd1, len_cd1, ~] = ZeroOnesCount(abs(neuron_sdf) > cd_1);
    [start_cd2, len_cd2, ~] = ZeroOnesCount(abs(neuron_sdf) > cd_2);
    
    % Then limit this to points in which it is above threshold for the defined
    % period of time
    clear modOnset_cd1 modOnset_cd2
    modOnset_cd1 = start_cd1(find(len_cd1 > cd_1t)) - timewins.zero;
    modOnset_cd2 = start_cd2(find(len_cd2 > cd_2t)) - timewins.zero;
    
    % We will then just limit the periods of modulation that occur within 0 and
    % 200 ms following target onset.
    modOnset_cd1 = modOnset_cd1(modOnset_cd1 > -500 & modOnset_cd1 < 0);
    
    if length(modOnset_cd1) > 1
        [minValue,closestIndex] = min(abs(modOnset_cd1-0));
        modOnset_cd1 = modOnset_cd1(closestIndex);
    end
    
    % After determining the onset falls within this criteria, we can then "walk
    % out" from this point and stop when we fall back below the threshold.
    trig_flag = 0; duration = 0;
    while trig_flag == 0
        duration = duration + 1;
        try
            trig_flag = abs(neuron_sdf(modOnset_cd1 + timewins.zero + duration)) < cd_1;
        catch
            trig_flag = 1;
        end
    end
    
    if ~isempty(modOnset_cd1)
        onset_cd_thresh(neuron_i,1) = modOnset_cd1;
        offset_cd_thresh(neuron_i,1) = modOnset_cd1+duration;
    else
        onset_cd_thresh(neuron_i,1) = NaN;
        offset_cd_thresh(neuron_i,1) = NaN;
    end
end

%% Figure: Plot visually responsive neurons
% Setup figure dimensions and calculate number of batches to show multiple
% figures on one sheet
n_plot_x = 4; n_plot_y = 3; n_plot_sheet = n_plot_x*n_plot_y;
n_batches = round(length(neuron_index.visual_pos)/n_plot_sheet,-1)+1;
neuron_list_i = 0;

% Then for each page
for page_i = 1:n_batches
    % Create the figure canvas
    fig_out = figure('Renderer', 'painters', 'Position', [100 100 1200 800]);
    
    % Then for each plot space
    for plot_i = 1:n_plot_sheet
        % Find the neuron id number
        neuron_list_i = neuron_list_i+1;
        try
            neuron_i = neuron_index.visual_pos(neuron_list_i);
            % Plot the mean activity between -200 and 500 ms relative to
            % target
            subplot(n_plot_x, n_plot_y, plot_i)
            plot(timewins.sdf, sdf_saccade_ns_zscore(neuron_i,:),'k');
            xlim([-500 1200]); vline(0,'k--'); vline(600,'k:'); hline(0,'k--')%, hline([-3 3],'r:'); hline([-6 6],'r:');
            % Draw the detected onset/offset in red (dotted)
            vline([onset_cd_thresh(neuron_i,1) offset_cd_thresh(neuron_i,1)],'r:')
            % and label our figures
            xlabel('Time from Saccade (ms)')
            ylabel('Firing rate (z-score)')
            title(['Neuron: ' int2str(neuron_i) ])
            
        catch
            continue
        end
        
    end
    
    % Once we're done with a page, save it and close it.
    filename = fullfile(dirs.root,'results','sdf_saccade_figs',['sdf_saccade_overview_pg' int2str(page_i) '.pdf']);
    set(fig_out,'PaperSize',[20 10]); %set the paper size to what you want
    print(fig_out,filename,'-dpdf') % then print it
    close(fig_out)
end

%% Organize: Manual curation adjustments
%  Note: these adjustments were made after visual inspection of SDF's with
%  ***describe properties here

% ONSET %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% (1): Adjust offset latency %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%      Plot figures and use pointer to find onset ----------------
for neuron_i = [9,10,11,45,47,53,84,85,93,103,106,113,114,116,120,150,160,222,289,298,318,331,...
        335,336,337,353,354,356,369,370,379,380,389,410,436,463]
    fig_out = figure('Renderer', 'painters', 'Position', [100 100 400 400]);
    
    plot(timewins.sdf, sdf_saccade_ns_zscore(neuron_i,:),'k');
    xlim([-750 250]); vline(0,'k--'); hline(0,'k--')%, hline([-3 3],'r:'); hline([-6 6],'r:');
    % Draw the detected onset in red (dotted)
    vline(onset_cd_thresh(neuron_i),'r:')
    xlabel('Time from Saccade (ms)')
    ylabel('Firing rate (z-score)')
    title(['Neuron: ' int2str(neuron_i) ])
end

% Update timings with manual markings (note: most were made at the
% 2sd mark).
onset_cd_thresh([9,10,11,45,47,53,84,85,93,103,106,113,114,116,120,150,160,222,289,298,318,331,...
    335,336,337,353,354,356,369,370,379,380,389,410,436,463],1)=...
    [-345,-305,-321,-324,-377,-371,-320,-388,-230,-218,-330,-229,-355,-358,-334,-164,...
    -242,-164,-302,-306,-272,-161,-344,-372,-314,-353,-350,-289,-270,-315,-351,...
    -343,-270,-225,-190,-242];

% OFFSET %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% (1): Adjust offset latency %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%      Plot figures and use pointer to find onset ----------------
for neuron_i = [53, 57, 62, 84, 85, 103, 106, 114, 120, 129, 150, 167, 295,...
        298, 337, 353, 354, 363, 368, 377, 378, 379, 392, 553]
    fig_out = figure('Renderer', 'painters', 'Position', [100 100 400 400]);
    
    plot(timewins.sdf, sdf_saccade_ns_zscore(neuron_i,:),'k');
    xlim([-200 500]); vline(0,'k--'); hline(0,'k--')%, hline([-3 3],'r:'); hline([-6 6],'r:');
    % Draw the detected onset in red (dotted)
    vline(offset_cd_thresh(neuron_i),'r:')
    xlabel('Time from Saccade (ms)')
    ylabel('Firing rate (z-score)')
    title(['Neuron: ' int2str(neuron_i) ])
end

% Update timings with manual markings (note: most were made at the
% 2sd mark).
offset_cd_thresh([53, 57, 62, 84, 85, 103, 106, 114, 120, 129, 150, 167, 295,...
    298, 337, 353, 354, 363, 368, 377, 378, 379, 392, 553],1)=...
    [110, 315, 332, 55, 14,136,89,273,328,119,22,344,351,361,302,...
    266,196,383,204,313,331,226,-31,47];

%% Organize: Save updated information to the main analysis table, and save it.

visual_info.saccade_onset = onset_cd_thresh;
visual_info.saccade_offset = offset_cd_thresh;
save(fullfile(dirs.root,'results','mat_files','visual_info.mat'),'visual_info')