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
    
    % Get the corresponding session for the neuron
    session_i = executiveBeh.neuronMatPosit(neuron_i,1);
    
    % Load in pre-processed SDF
    spk_in = load(fullfile(dirs.data_spike,['SDF_' int2str(neuron_i)]));
    
    % Define input trials and windows
    trial_in = []; trial_in = executiveBeh.ttx.GO{session_i};
    baseline_win = []; baseline_win = timewins.baseline_target + timewins.zero;
    target_win = []; target_win = timewins.visual_target + timewins.zero;
    
    % Get baseline mean and std for z-score
    fr_visual_baseline_mean = nanmean(nanmean(spk_in.SDF.target(trial_in,baseline_win)));
    fr_visual_baseline_std = nanstd(nanmean(spk_in.SDF.target(trial_in,baseline_win)));
    
    % Get the mean SDF on no-stop trials.
    if fr_visual_baseline_std == 0
        sdf_visual_raw(neuron_i,:) = nan(1,length(timewins.sdf));
        sdf_visual_zscore(neuron_i,:) = nan(1,length(timewins.sdf));
    else
        sdf_visual_raw(neuron_i,:) = nanmean(spk_in.SDF.target(trial_in,:));
        sdf_visual_zscore(neuron_i,:) = (sdf_visual_raw(neuron_i,:) - fr_visual_baseline_mean)./...
            fr_visual_baseline_std;
    end
    
    % Analysis: Compare mean FR between baseline and visual period
    fr_comp_visual_baseline_mean = []; fr_comp_visual_target_mean = [];
    fr_comp_visual_baseline_mean = nanmean(spk_in.SDF.target(trial_in,baseline_win),2);
    fr_comp_visual_target_mean = nanmean(spk_in.SDF.target(trial_in,target_win),2);

    fr_comp_visual_sig_h(neuron_i,1) = ttest2(fr_comp_visual_baseline_mean, fr_comp_visual_target_mean);
    fr_comp_visual_sig_dir(neuron_i,1) = mean(fr_comp_visual_baseline_mean) < mean(fr_comp_visual_target_mean);
    
end

% Get a normalised (max z-score) density function for heatmap
for neuron_i = 1:575
    target_win = []; target_win = timewins.visual_target + timewins.zero;
    
    fr_visual_max = max(abs(sdf_visual_zscore(neuron_i,target_win)));
    sdf_visual_maxnorm(neuron_i,:) = sdf_visual_zscore(neuron_i,:)./fr_visual_max;
    
end


%% Extract: Produce summary sheet figure for mean SDF

n_plot_x = 4; n_plot_y = 3; n_plot_sheet = n_plot_x*n_plot_y;
n_batches = round(575/n_plot_sheet,-1)+1;

neuron_i = 0;
for page_i = 1:n_batches
    fig_out = figure('Renderer', 'painters', 'Position', [100 100 1200 800]);
    
    for plot_i = 1:n_plot_sheet
        neuron_i = neuron_i+1;
        try
            subplot(n_plot_x, n_plot_y, plot_i)
            plot(timewins.sdf, sdf_visual_zscore(neuron_i,:),'k')
            xlim([-200 500]); vline(0,'k--'); hline(0,'k--')%, hline([-3 3],'r:'); hline([-6 6],'r:'); 
            xlabel('Time from Target (ms)')
            ylabel('Firing rate (z-score)')
            title(['Neuron: ' int2str(neuron_i) ])
            
        catch
            continue
        end
        
    end
    
    filename = fullfile(dirs.root,'results','sdf_overview_figs',['sdf_target_overview_pg' int2str(page_i) '.pdf']);
    set(fig_out,'PaperSize',[20 10]); %set the paper size to what you want
    print(fig_out,filename,'-dpdf') % then print it
    close(fig_out)
end

%% Analysis: Probe deviation from a baseline period during target onset period
cd_1 = 3; cd_1t = 75; % Signal must be above 3 sd for 100 ms,
cd_2 = 6; cd_2t = 25;  % ...and above 6 sd for 25 ms.

clear onset_cd_thresh offset_cd_thresh
% Get the neuron and z-scored activity
for neuron_i = 1:575
    neuron_sdf = []; neuron_sdf = sdf_visual_zscore(neuron_i,:);
    
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
    modOnset_cd1 = modOnset_cd1(modOnset_cd1 > 0 & modOnset_cd1 < 150);
    
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

%% Analysis: Find peaks and onset/offset
clear onset_peak_thresh peak_peak_thresh offset_peak_thresh
for neuron_i = 1:575
    % Get the z-scored activity for the neuron of interest
    neuron_sdf = []; neuron_sdf = sdf_visual_zscore(neuron_i,:);
    
    % Identify peak periods, and find the absolute peak activity
    [pks,locs] = findpeaks(neuron_sdf);
    peak_time = locs(find(pks == max(abs(pks))))- timewins.zero;
    peak_time = peak_time(peak_time > 0 & peak_time < 150);
    
    if ~isempty(peak_time)
        % Once we've identified the peak, then we can "walk out" forwards and
        % backwards to find the onset and offset, when it falls below the
        % initial threshold.
        trig_flag_start = 0; trig_flag_end = 0; duration_start = 0; duration_end = 0;
        
        % For the onset, we will walk backwards from the identified peak
        while trig_flag_start == 0
            duration_start = duration_start - 1;
            try
                trig_flag_start = abs(neuron_sdf(peak_time + timewins.zero + duration_start)) < cd_1;
            catch % If we go beyond the range of the SDF, we will bail out of the while loop.
                trig_flag_start = 1;
            end
        end
        
        % For the offset, we will walk forwards from the identified peak
        while trig_flag_end == 0
            duration_end = duration_end + 1;
            try
                trig_flag_end = abs(neuron_sdf(peak_time + timewins.zero + duration_end)) < cd_1;
            catch % If we go beyond the range of the SDF, we will bail out of the while loop.
                trig_flag_end = 1;
            end
        end
        
        onset_peak_thresh(neuron_i,1) = duration_start+peak_time;
        peak_peak_thresh(neuron_i,1) = peak_time;
        offset_peak_thresh(neuron_i,1) = duration_end+peak_time;
    else
        onset_peak_thresh(neuron_i,1) = NaN;
        peak_peak_thresh(neuron_i,1) = NaN;
        offset_peak_thresh(neuron_i,1) = NaN;
    end
    
end


%% Organize: Compile analyses into one table to look for convergent evidence of visual activty
% Generate table
clear visual_info
for neuron_i = 1:575
    
    session = executiveBeh.neuronMatPosit(neuron_i,1);
    session_label = FileNames(session);
    monkey_label = executiveBeh.nhpSessions.monkeyNameLabel(session);
    spkWidth = executiveBeh.bioInfo.spkWidth(neuron_i)*25;
    site = executiveBeh.bioInfo.sessionSite(session);
    
    visual_flag = fr_comp_visual_sig_h(neuron_i);
    visual_fac = fr_comp_visual_sig_dir(neuron_i) == 1;
    visual_sup = fr_comp_visual_sig_dir(neuron_i) == 0;
    
    cd_detect_onset = onset_cd_thresh(neuron_i);
    cd_detect_offset = offset_cd_thresh(neuron_i);
    
    peak_detect_onset = onset_peak_thresh(neuron_i);
    peak_detect_offset = offset_peak_thresh(neuron_i);
    peak_detect_peak = peak_peak_thresh(neuron_i);
    
    visual_info(neuron_i,:) = table(neuron_i,session,session_label,...
        monkey_label,spkWidth,site,visual_flag,visual_fac,visual_sup,...
        cd_detect_onset,cd_detect_offset,...
        peak_detect_onset,peak_detect_offset,peak_detect_peak);
    
end


%% Organize: Get indices of visually responsive neurons
% Find neurons that meet a defined criteria
%   Here, we are looking at neurons that have a significantly different
%   firing rate between a baseline and visual period, and also hit the
%   criteria to allow for time extraction (i.e activity > 2SD / 6SD).
neuron_index.visual = find(visual_info.visual_flag == 1 & ...
    ~isnan(visual_info.cd_detect_onset));

neuron_index.nonvisual = find(~ismember(1:575,neuron_index.visual));

%% Figure: Plot visually responsive neurons
% Setup figure dimensions and calculate number of batches to show multiple
% figures on one sheet
n_plot_x = 4; n_plot_y = 3; n_plot_sheet = n_plot_x*n_plot_y;
n_batches = round(length(neuron_index.visual)/n_plot_sheet,-1)+1;
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
            neuron_i = neuron_index.visual(neuron_list_i);
            % Plot the mean activity between -200 and 500 ms relative to
            % target
            subplot(n_plot_x, n_plot_y, plot_i)
            plot(timewins.sdf, sdf_visual_zscore(neuron_i,:),'k');
            xlim([-200 500]); vline(0,'k--'); hline(0,'k--')%, hline([-3 3],'r:'); hline([-6 6],'r:'); 
            % Draw the detected onset in red (dotted)
            vline([visual_info.cd_detect_onset(neuron_i) visual_info.cd_detect_offset(neuron_i)],'r:')
            % and label our figures
            xlabel('Time from Target (ms)')
            ylabel('Firing rate (z-score)')
            title(['Neuron: ' int2str(neuron_i) ])
            
        catch
            continue
        end
        
    end
    
    % Once we're done with a page, save it and close it.
    filename = fullfile(dirs.root,'results','sdf_visual_figs',['sdf_target_overview_pg' int2str(page_i) '.pdf']);
    set(fig_out,'PaperSize',[20 10]); %set the paper size to what you want
    print(fig_out,filename,'-dpdf') % then print it
    close(fig_out)
end

%% Organize: Manual curation adjustments
%  Note: these adjustments were made after visual inspection of SDF's with
%  properties derived from a 3/6 sd cutoff, at 75/25 ms thresholds, with a
%  latest onset of 150 ms.

% (1): Incorrectly identified visual-neurons %%%%%%%%%%%%%%%%%%%%%%%%%
visual_info.visual_flag([49,52,68,133,382,402],:) = 0;

% (2): Adjust onset latency %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%       Plot figures and use pointer to find onset ----------------
for neuron_i = [80,83,103,125,160,162,163,202,269,270,276,331]
    fig_out = figure('Renderer', 'painters', 'Position', [100 100 400 400]);
    
    plot(timewins.sdf, sdf_visual_zscore(neuron_i,:),'k');
    xlim([-200 500]); vline(0,'k--'); hline(0,'k--')%, hline([-3 3],'r:'); hline([-6 6],'r:');
    % Draw the detected onset in red (dotted)
    vline(visual_info.cd_detect_onset(neuron_i),'r:')
    xlabel('Time from Target (ms)')
    ylabel('Firing rate (z-score)')
    title(['Neuron: ' int2str(neuron_i) ])
end
%       Update timings with manual markings (note: most were made at the
%       2sd mark).
visual_info.cd_detect_onset...
    ([80,83,103,125,160,162,163,202,269,270,276,331]) =...
     [85,74, 75, 38, 52, 36,125, 71, 82, 92, 91,113];

% (3): Update visual neuron index accordingly %%%%%%%%%%%%%%%%%%%%%%%%
neuron_index.visual = []; neuron_index.nonvisual = [];

neuron_index.visual = find(visual_info.visual_flag == 1 & ...
    ~isnan(visual_info.cd_detect_onset));

neuron_index.nonvisual = find(~ismember(1:575,neuron_index.visual));


%% Figure: Plot histogram of spike widths
spkwidth_figure = figure('Renderer', 'painters', 'Position', [100 100 400 200]);
histogram(executiveBeh.bioInfo.spkWidth*25,0:25:800,'LineStyle','None')
vline(250,'k--')
% Once we're done with a page, save it and close it.
filename = fullfile(dirs.root,'results','gen_figures','spkwidth_histogram.pdf');
set(spkwidth_figure,'PaperSize',[20 10]); %set the paper size to what you want
print(spkwidth_figure,filename,'-dpdf') % then print it
close(spkwidth_figure)

%% Organize: output data
save(fullfile(dirs.root,'results','mat_files','neuron_index.mat'),'neuron_index')
save(fullfile(dirs.root,'results','mat_files','visual_info.mat'),'visual_info')
save(fullfile(dirs.root,'results','mat_files','visual_sdf_1.mat'),'sdf_visual_raw','sdf_visual_zscore','sdf_visual_maxnorm')

