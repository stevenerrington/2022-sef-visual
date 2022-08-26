%% Extract: Get average spiking activity around target onset.

parfor neuron_i = 1:575
    fprintf('Extracting target aligned SDF for neuron %i of %i. \n',neuron_i, 575)
    
    session_i = executiveBeh.neuronMatPosit(neuron_i,1);
    
    spk_in = load(fullfile(dirs.data_spike,['SDF_' int2str(neuron_i)]));
    
    trial_in = []; trial_in = executiveBeh.ttx.GO{session_i};
    baseline_win = []; baseline_win = timewins.baseline_target + timewins.zero;
    target_win = []; target_win = timewins.visual_target + timewins.zero;
    
    fr_visual_baseline_mean = nanmean(nanmean(spk_in.SDF.target(trial_in,baseline_win)));
    fr_visual_baseline_std = nanstd(nanmean(spk_in.SDF.target(trial_in,baseline_win)));
    
    
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
cd_1 = 3; cd_1t = 100; % Signal must be above 2 sd for 100 ms,
cd_2 = 6; cd_2t = 25;  % ...and above 6 sd for 10 ms.

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
    modOnset_cd1 = modOnset_cd1(modOnset_cd1 > 0 & modOnset_cd1 < 200);
    
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
    peak_time = peak_time(peak_time > 0 & peak_time < 200);
    
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
clear visual_extract
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
    
    visual_extract(neuron_i,:) = table(neuron_i,session,session_label,...
        monkey_label,spkWidth,site,visual_flag,visual_fac,visual_sup,...
        cd_detect_onset,cd_detect_offset,...
        peak_detect_onset,peak_detect_offset,peak_detect_peak);
    
end



%% Organize: Get indices of visually responsive neurons
% Find neurons that meet a defined criteria
%   Here, we are looking at neurons that have a significantly different
%   firing rate between a baseline and visual period, and also hit the
%   criteria to allow for time extraction (i.e activity > 2SD / 6SD).
neuron_index.visual = find(visual_extract.visual_flag == 1 & ...
    ~isnan(visual_extract.cd_detect_onset));

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
            vline(visual_extract.onset_cd_thresh(neuron_i),'r:')
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

%% Figure: Heatmap comparing visual and non-visual neurons

neuron_index.visual_pos = neuron_index.visual(visual_extract.visual_fac(neuron_index.visual) == 1);
neuron_index.visual_neg = neuron_index.visual(visual_extract.visual_sup(neuron_index.visual) == 1);
[~, latency_order_pos] = sort(visual_extract.cd_detect_onset(neuron_index.visual_pos),'ascend');
[~, latency_order_neg] = sort(visual_extract.cd_detect_onset(neuron_index.visual_neg),'descend');

neuron_order = [neuron_index.visual_pos(latency_order_pos); neuron_index.visual_neg(latency_order_neg)];

fig_out = figure('Renderer', 'painters', 'Position', [100 100 400 600]);
subplot(2,1,1)
imagesc('XData',timewins.sdf,'YData',1:length(neuron_index.visual), 'CData',sdf_visual_maxnorm(neuron_order,:))
xlim([-200 500]); ylim([1 400]);  caxis([-1,1]); colormap(parula); vline(0,'k')
xlabel('Time from Target (ms)'); ylabel('Neuron')

subplot(2,1,2)
imagesc('XData',timewins.sdf,'YData',1:length(neuron_index.nonvisual), 'CData',sdf_visual_maxnorm(neuron_index.nonvisual,:))
xlim([-200 500]); ylim([1 400]); caxis([-1,1]); colormap(parula); vline(0,'k')
xlabel('Time from Target (ms)'); 

% Once we're done with a page, save it and close it.
filename = fullfile(dirs.root,'results','gen_figures','population_figure_heatmap.pdf');
set(fig_out,'PaperSize',[20 10]); %set the paper size to what you want
print(fig_out,filename,'-dpdf') % then print it
close(fig_out)

%% Figure:
clear input_sdf population_sdf
input_sdf = num2cell(sdf_visual_zscore, 2);

for neuron_i = 1:575
    if ismember(neuron_i,neuron_index.visual_pos)
        class_label{neuron_i,1} = '1_visual_pos';
    elseif ismember(neuron_i,neuron_index.visual_neg)
        class_label{neuron_i,1} = '2_visual_neg';
    else
        class_label{neuron_i,1} = '3_other';
    end
end

population_sdf(1,1)=gramm('x',timewins.sdf,'y',input_sdf,'color',class_label);
population_sdf(1,1).stat_summary();
population_sdf(1,1).axe_property('XLim',[-200 500],'YLim',[-15 15]);
population_sdf(1,1).set_names('x','Time from Target (ms)','y','FR (Z-score)');

pop_figure = figure('Renderer', 'painters', 'Position', [100 100 600 300]);
population_sdf.draw();

% Once we're done with a page, save it and close it.
filename = fullfile(dirs.root,'results','gen_figures','population_figure_visual.pdf');
set(pop_figure,'PaperSize',[20 10]); %set the paper size to what you want
print(pop_figure,filename,'-dpdf') % then print it
close(pop_figure)

%% %% Analysis; clustering visual neurons
% 
% timeWin = [-250:500];
% timeWin_inputSDF = timewins.zero + timeWin;
% inputSDF = {sdf_visual_raw(:,timeWin_inputSDF), sdf_visual_raw(:,timeWin_inputSDF)};
% 
% sdfTimes = {timeWin, timeWin}; sdfEpoch = {timeWin, timeWin};
% 
% 
% [sortIDs,idxDist, raw, respSumStruct, rawLink,myK] =...
%     consensusCluster(inputSDF,sdfTimes,'-e',sdfEpoch);
% 
% 
% 
% normResp = scaleResp(inputSDF,sdfTimes,'max');
% 
% nClusters_manual = 3; clusterNeurons = [];
% for i = 1:nClusters_manual
%     clusterNeurons{i} = find(sortIDs(:,nClusters_manual) == i );
% end
% 
% % Plot clustering output %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % Dendrogram
% 
% figure('Renderer', 'painters', 'Position', [100 100 500 400]);
% 
% subplot(1,5,5);
% for ir = 1:size(raw,1)
%     for ic = (ir+1):size(raw,2)
%         raw(ic,ir) = raw(ir,ic);
%     end
% end
% [h,~,outPerm] = dendrogram(rawLink,0,'Orientation','right');
% set(gca,'YDir','Reverse');
% klDendroClustChange(h,rawLink,sortIDs(:,nClusters_manual))
% set(gca,'YTick',[]); xlabel('Similarity')
% subplot(1,5,[1:4]);
% 
% imagesc(raw(outPerm,outPerm));
% colormap(gray);
% xlabel('Unit Number'); set(gca,'YAxisLocation','Left');
% xticks([50:50:500]); yticks([50:50:500])
% 
% close all
% 
% 
% 










