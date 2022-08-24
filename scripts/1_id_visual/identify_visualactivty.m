%% Extract: Get average spiking activity around target onset.

parfor neuron_i = 1:575
    fprintf('Extracting target aligned SDF for neuron %i of %i. \n',neuron_i, 575)
    
    session_i = executiveBeh.neuronMatPosit(neuron_i,1);
    
    spk_in = load(fullfile(dirs.data_spike,['SDF_' int2str(neuron_i)]));
    
    trial_in = []; trial_in = executiveBeh.ttx.GO{session_i};
    epoch_win = []; epoch_win = timewins.baseline_target + timewins.zero;
    
    fr_visual_baseline_mean = nanmean(nanmean(spk_in.SDF.target(trial_in,epoch_win)));
    fr_visual_baseline_std = nanstd(nanmean(spk_in.SDF.target(trial_in,epoch_win)));
    
    if fr_visual_baseline_std == 0
        sdf_visual_raw(neuron_i,:) = nan(1,length(timewins.sdf));
        sdf_visual_zscore(neuron_i,:) = nan(1,length(timewins.sdf));
    else
        sdf_visual_raw(neuron_i,:) = nanmean(spk_in.SDF.target(trial_in,:));
        sdf_visual_zscore(neuron_i,:) = (sdf_visual_raw(neuron_i,:) - fr_visual_baseline_mean)./...
            fr_visual_baseline_std;
    end
    
    
end

%% Analysis: Compare mean firing rate between baseline and visual period
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

%% Analysis: Find peaks
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


%%
table(onset_cd_thresh, offset_cd_thresh, onset_peak_thresh, peak_peak_thresh, offset_peak_thresh)



%% Analysis; clustering visual neurons

timeWin = [-250:500];
timeWin_inputSDF = timewins.zero + timeWin;
inputSDF = {sdf_visual_raw(:,timeWin_inputSDF), sdf_visual_raw(:,timeWin_inputSDF)};

sdfTimes = {timeWin, timeWin}; sdfEpoch = {timeWin, timeWin};


[sortIDs,idxDist, raw, respSumStruct, rawLink,myK] =...
    consensusCluster(inputSDF,sdfTimes,'-e',sdfEpoch);



normResp = scaleResp(inputSDF,sdfTimes,'max');

nClusters_manual = 3; clusterNeurons = [];
for i = 1:nClusters_manual
    clusterNeurons{i} = find(sortIDs(:,nClusters_manual) == i );
end

% Plot clustering output %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Dendrogram

figure('Renderer', 'painters', 'Position', [100 100 500 400]);

subplot(1,5,5);
for ir = 1:size(raw,1)
    for ic = (ir+1):size(raw,2)
        raw(ic,ir) = raw(ir,ic);
    end
end
[h,~,outPerm] = dendrogram(rawLink,0,'Orientation','right');
set(gca,'YDir','Reverse');
klDendroClustChange(h,rawLink,sortIDs(:,nClusters_manual))
set(gca,'YTick',[]); xlabel('Similarity')
subplot(1,5,[1:4]);

imagesc(raw(outPerm,outPerm));
colormap(gray);
xlabel('Unit Number'); set(gca,'YAxisLocation','Left');
xticks([50:50:500]); yticks([50:50:500])

close all













