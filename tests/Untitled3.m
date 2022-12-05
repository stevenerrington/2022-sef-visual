timewins.sdf = [-1000:2000];
timewins.zero = abs(timewins.sdf(1));
timewins.ssrt_baseline = [-100:0];
timewins.ssrt_comp = [-50:250];

%% Extract: Get latency-matched SDF for non-canceled trials
parfor neuron_i = 1:size(acc_map_info,1)
    
    neuralFilename = acc_map_info.session{neuron_i};
    neuronLabel = acc_map_info.unit{neuron_i};
    fprintf('Extracting: %s - %s ... [%i of %i]  \n',neuralFilename,neuronLabel,neuron_i,size(acc_map_info,1))
    
    %... and find the corresponding behavior file index
    behFilename = data_findBehFile(neuralFilename);
    behaviorIdx = find(strcmp(dataFiles_beh,behFilename(1:end-4)));
    
    % Load in pre-processed spike data
    data_in = load(fullfile(dirs.root,'data','SDF',...
        [neuralFilename '_SDF_' neuronLabel '.mat']));
    
    sdf_canceled_targetx = []; sdf_nostop_targetx = [];
    sdf_canceled_ssdx = []; sdf_nostop_ssdx = [];
    sdf_canceled_ssrtx = []; sdf_nostop_ssrtx = [];
    sdf_canceled_tonex = []; sdf_nostop_tonex = [];
    
    % Latency match
    for ssd_i = 1:length(behavior(behaviorIdx).stopSignalBeh.inh_SSD)
        trl_canceled = []; trl_canceled = behavior(behaviorIdx).ttm.C.C{ssd_i};
        trl_nostop = []; trl_nostop = behavior(behaviorIdx).ttm.C.GO{ssd_i};
        
        if length(trl_canceled) < 10 | length(trl_nostop) < 10
            sdf_canceled_ssdx(ssd_i,:) = nan(1,length(-1000:2000));
            sdf_nostop_ssdx(ssd_i,:) = nan(1,length(-1000:2000));
        else
            % Target:
            sdf_canceled_targetx(ssd_i,:) = nanmean(data_in.SDF.target(trl_canceled,:));
            sdf_nostop_targetx(ssd_i,:) = nanmean(data_in.SDF.target(trl_nostop,:));
            % SSD:
            sdf_canceled_ssdx(ssd_i,:) = nanmean(data_in.SDF.stopSignal_artifical(trl_canceled,:));
            sdf_nostop_ssdx(ssd_i,:) = nanmean(data_in.SDF.stopSignal_artifical(trl_nostop,:));
            % SSRT:
            sdf_canceled_ssrtx(ssd_i,:) = nanmean(data_in.SDF.ssrt(trl_canceled,:));
            sdf_nostop_ssrtx(ssd_i,:) = nanmean(data_in.SDF.ssrt(trl_nostop,:));
            % Tone:
            sdf_canceled_tonex(ssd_i,:) = nanmean(data_in.SDF.tone(trl_canceled,:));
            sdf_nostop_tonex(ssd_i,:) = nanmean(data_in.SDF.tone(trl_nostop,:));
        end
    end
    
    % Get mean SDF for:
    % - Target
    sdf_canceled_all_target(neuron_i,:) = nanmean(sdf_canceled_targetx);
    sdf_nostop_all_target(neuron_i,:) = nanmean(sdf_nostop_targetx);
    % - Stop-signal
    sdf_canceled_all_stopsignal(neuron_i,:) = nanmean(sdf_canceled_ssdx);
    sdf_nostop_all_stopsignal(neuron_i,:) = nanmean(sdf_nostop_ssdx);
    % - SSRT
    sdf_canceled_all_ssrt(neuron_i,:) = nanmean(sdf_canceled_ssrtx);
    sdf_nostop_all_ssrt(neuron_i,:) = nanmean(sdf_nostop_ssrtx);
    % - Tone
    sdf_canceled_all_tone(neuron_i,:) = nanmean(sdf_canceled_tonex);
    sdf_nostop_all_tone(neuron_i,:) = nanmean(sdf_nostop_tonex);
    
    
    % Save SSD specific activity, just incase.
    sdf_canceled_ssd{neuron_i} = sdf_canceled_ssdx;
    sdf_nostop_ssd{neuron_i} = sdf_nostop_ssdx;
    
end

%% Analyse: Compare activity for modulation post-stopping
parfor neuron_i = 1:size(acc_map_info,1)
    
    neuralFilename = acc_map_info.session{neuron_i};
    neuronLabel = acc_map_info.unit{neuron_i};
    fprintf('Extracting: %s - %s ... [%i of %i]  \n',neuralFilename,neuronLabel,neuron_i,size(acc_map_info,1))
    
    %... and find the corresponding behavior file index
    behFilename = data_findBehFile(neuralFilename);
    behaviorIdx = find(strcmp(dataFiles_beh,behFilename(1:end-4)));
    
    % Load in pre-processed spike data
    data_in = load(fullfile(dirs.root,'data','SDF',...
        [neuralFilename '_SDF_' neuronLabel '.mat']));
    
    trl_canceled = []; trl_nostop = [];
    trl_canceled = behavior(behaviorIdx).ttx.canceled.all.all;
    trl_nostop = behavior(behaviorIdx).ttx.nostop.all.all;
    
    canc_sdf_mean = []; nostop_sdf_mean = [];
    canc_sdf_mean = nanmean(data_in.SDF.stopSignal_artifical(trl_canceled,timewins.ssrt_comp+timewins.zero),2);
    nostop_sdf_mean = nanmean(data_in.SDF.stopSignal_artifical(trl_nostop,timewins.ssrt_comp+timewins.zero),2);
    
    [h,p,~,stats] = ttest2(canc_sdf_mean,nostop_sdf_mean);
    mean_canceled = nanmean(canc_sdf_mean); mean_nostop = nanmean(nostop_sdf_mean);
    ssrt_dir = mean_canceled > mean_nostop;
    
    ssrt_table(neuron_i,:) = table({neuralFilename},{neuronLabel},h,stats,ssrt_dir,mean_canceled,mean_nostop,...
        'VariableNames',{'neuralFilename','unit','h','stats','c_dir','mean_c','mean_ns'});
    
end



ssrt_neurons = find(ssrt_table.h == 1);


%% Analyse: Clustering approach for errors
sdfWindow = timewins.sdf;
blWindow = [-100:0];
inputNeurons = ssrt_neurons;
inputSDF_target = {sdf_canceled_all_target(inputNeurons,:),sdf_nostop_all_target(inputNeurons,:)};
inputSDF_ssd = {sdf_canceled_all_stopsignal(inputNeurons,:),sdf_nostop_all_stopsignal(inputNeurons,:)};
inputSDF_ssrt = {sdf_canceled_all_ssrt(inputNeurons,:),sdf_nostop_all_ssrt(inputNeurons,:)};
inputSDF_tone = {sdf_canceled_all_tone(inputNeurons,:),sdf_nostop_all_tone(inputNeurons,:)};

sdfTimes = {sdfWindow, sdfWindow};
sdfEpoch = {[-200:600],[-200:600]};

colorMapping = [1,2];

[sortIDs,idxDist, raw, respSumStruct, rawLink,myK] =...
    consensusCluster(inputSDF_ssrt,sdfTimes,'-e',sdfEpoch,'-ei',colorMapping);
normResp_target = scaleResp(inputSDF_target,sdfTimes,'max','-bl',blWindow);
normResp_stopSignal = scaleResp(inputSDF_ssd,sdfTimes,'max','-bl',blWindow);
normResp_ssrt = scaleResp(inputSDF_ssrt,sdfTimes,'max','-bl',blWindow);
normResp_tone = scaleResp(inputSDF_tone,sdfTimes,'max','-bl',blWindow);

nClusters_manual = 4; clusterNeurons = [];
for cluster_i = 1:nClusters_manual
    clusterNeurons{cluster_i} = find(sortIDs(:,nClusters_manual) == cluster_i );
end

% Figure: Dendrogram %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure('Renderer', 'painters', 'Position', [100 100 500 400]);

subplot(1,5,5)
[h,~,outPerm] = dendrogram(rawLink,0,'Orientation','right');
set(gca,'YDir','Reverse');
klDendroClustChange(h,rawLink,sortIDs(:,nClusters_manual))
set(gca,'YTick',[]); xlabel('Similarity')

subplot(1,5,[1:4]);
for ir = 1:size(raw,1)
    for ic = (ir+1):size(raw,2)
        raw(ic,ir) = raw(ir,ic);
    end
end
imagesc(raw(outPerm,outPerm));
colormap(gray);
xlabel('Unit Number'); set(gca,'YAxisLocation','Left');
xticks([0:100:end]); yticks([0:100:end])

% Figure: Cluster population SDF
norm_sdf_list = {normResp_target,normResp_stopSignal,normResp_ssrt,normResp_tone};
xlim_list = {[-500 250],[-250 500],[-250 500],[-500 250]};
ylim_list = {[-0.5 1.5],[-0.5 1.5],[-1 0.5],[-1 0.5]};
% Cluster main
figure('Renderer', 'painters', 'Position', [100 100 1200 900]);hold on
for epoch_i = 1:4
    clear norm_in
    norm_in = norm_sdf_list{epoch_i};
    
    for cluster_i = 1:nClusters_manual
        subplot_pos = epoch_i+(4*(cluster_i-1));
        subplot(nClusters_manual,4,subplot_pos); hold on
        plot(sdfTimes{1},nanmean(norm_in{1}(clusterNeurons{cluster_i},:),1), 'color', [colors.canceled]);
        plot(sdfTimes{2},nanmean(norm_in{2}(clusterNeurons{cluster_i},:),1), 'color', [colors.nostop]);
        vline(0, 'k--'); xlim(xlim_list{epoch_i});ylim(ylim_list{cluster_i})
        
        title(['Cluster ' int2str(cluster_i) ' - n: ' int2str(length(clusterNeurons{cluster_i}))])
        
    end
end

%% Analyse: Spike width
spk_width_cutoff = 250;
figure('Renderer', 'painters', 'Position', [100 100 200 800]);hold on

for cluster_i = 1:nClusters_manual
    subplot(nClusters_manual,1,cluster_i)
    pie([sum(abs(acc_map_info.spk_width(inputNeurons(clusterNeurons{cluster_i}))) >= spk_width_cutoff),...
        sum(abs(acc_map_info.spk_width(inputNeurons(clusterNeurons{cluster_i}))) < spk_width_cutoff)]);
end

%% Analyse: dorsal and ventral MCC

ml_axis_range = [-8:1:8]; ap_axis_range = [20:40]; depth_range = [-3000:50:3000];
ml_axis_zero = abs(min(ml_axis_range)); ap_axis_zero = abs(min(ap_axis_range)); depth_zero = abs(min(depth_range));

figure('Renderer', 'painters', 'Position', [100 100 1400 300]);hold on
for cluster_i = 1:nClusters_manual
    cluster_neurons = []; cluster_neurons = inputNeurons(clusterNeurons{cluster_i});
    
    clear sample_coords
    sample_array = zeros(length(ml_axis_range),length(ap_axis_range), length(depth_range));
    
    for neuron_i = 1:length(cluster_neurons)
        ap_sample_ref = find(ap_axis_range == acc_map_info.ap(cluster_neurons(neuron_i)));
        ml_sample_ref = find(ml_axis_range == acc_map_info.ml(cluster_neurons(neuron_i)));
        depth_sample_ref = find(depth_range == acc_map_info.depth(cluster_neurons(neuron_i)));
        
        sample_array(ap_sample_ref,ml_sample_ref,depth_sample_ref) =...
            sample_array(ap_sample_ref,ml_sample_ref,depth_sample_ref) + 1;
        
        sample_coords(neuron_i,:) = ...
            [acc_map_info.ap(cluster_neurons(neuron_i)),...
            acc_map_info.ml(cluster_neurons(neuron_i)),...
            acc_map_info.depth(cluster_neurons(neuron_i))];
    end
    
    subplot(1,4,cluster_i); hold on
    scatter3(sample_coords(:,1),sample_coords(:,2),sample_coords(:,3),'filled')
    view(71,2.5); ylim([-8 8]); xlim([25 35]); zlim([-3000 3000])
    set(gca,'ZDir','reverse')
    xl = xlim; yl = ylim;
    patch([[1 1]*xl(1) [1 1]*xl(2)], [yl fliplr(yl)], [1 1 1 1]*0, 'k', 'FaceAlpha',0.1);
    
    grid on
end

%% Analyse: activity as a function of SSD
count = 0;
for cluster_i = 1:nClusters_manual
    neuron_list = [];
    neuron_list = inputNeurons(clusterNeurons{cluster_i});
    
    for neuron_list_i = 1:length(neuron_list)
        neuron_i = neuron_list(neuron_list_i);
        count = count + 1;
        behFilename = data_findBehFile(acc_map_info.session{neuron_list(neuron_list_i)});
        behaviorIdx = find(strcmp(dataFiles_beh,behFilename(1:end-4)));
        
        mid_ssd_idx = behavior(behaviorIdx).stopSignalBeh.midSSDidx;
        
        conflict_epoch = [0:100]+timewins.zero;
        
        early_ssd_activity = mean(sdf_canceled_ssd{neuron_i}(mid_ssd_idx-1,conflict_epoch));
        mid_ssd_activity = mean(sdf_canceled_ssd{neuron_i}(mid_ssd_idx,conflict_epoch));
        late_ssd_activity = mean(sdf_canceled_ssd{neuron_i}(mid_ssd_idx+1,conflict_epoch));
        
        all_ssd_activity_norm = [early_ssd_activity;mid_ssd_activity;late_ssd_activity] ./...
            max([early_ssd_activity;mid_ssd_activity;late_ssd_activity]);
        
        
        early_ssd_ssd = behavior(behaviorIdx).stopSignalBeh.inh_SSD(mid_ssd_idx-1);
        mid_ssd_ssd = behavior(behaviorIdx).stopSignalBeh.inh_SSD(mid_ssd_idx);
        late_ssd_ssd = behavior(behaviorIdx).stopSignalBeh.inh_SSD(mid_ssd_idx+1);
        
        early_ssd_pnc = behavior(behaviorIdx).stopSignalBeh.inh_pnc(mid_ssd_idx-1);
        mid_ssd_pnc = behavior(behaviorIdx).stopSignalBeh.inh_pnc(mid_ssd_idx);
        late_ssd_pnc = behavior(behaviorIdx).stopSignalBeh.inh_pnc(mid_ssd_idx+1);
        
        
        functional_stopping_activity(count,:) =...
            table(neuron_i,cluster_i,...
            all_ssd_activity_norm(1), all_ssd_activity_norm(2), all_ssd_activity_norm(3),...
            early_ssd_ssd,mid_ssd_ssd,late_ssd_ssd,early_ssd_pnc,mid_ssd_pnc,late_ssd_pnc,...
            'VariableNames',{'neuron','cluster','activity_early','activity_mid','activity_late',...
            'ssd_early','ssd_mid','ssd_late','pnc_early','pnc_mid','pnc_late'});
        
    end
end

%% Figure: Conflict Activity
cluster_i = 1;

plotData = []; plotDataY = []; labelData = [];

plotData=...
    [functional_stopping_activity.activity_early;...
    functional_stopping_activity.activity_mid;...
    functional_stopping_activity.activity_late];

plotDataY=...
    [functional_stopping_activity.pnc_early;...
    functional_stopping_activity.pnc_mid;...
    functional_stopping_activity.pnc_late];

labelData=...
    [repmat({'1_early'},length(functional_stopping_activity.activity_early),1);...
    repmat({'2_mid'},length(functional_stopping_activity.activity_mid),1);...
    repmat({'3_late'},length(functional_stopping_activity.activity_late),1)];

clusterData = repmat(functional_stopping_activity.cluster,3,1);

testFigure_out = figure('Renderer', 'painters', 'Position', [100 100 1000 200]);
clear testFigure
testFigure(1,1)= gramm('x',labelData,'y',plotData,'color',labelData);
testFigure(1,1).stat_summary('geom',{'bar','line','black_errorbar'});
testFigure.axe_property('YLim',[0.8 1]);
testFigure.facet_grid([],clusterData);
testFigure.draw();

testFigure_out2 = figure('Renderer', 'painters', 'Position', [100 100 1000 200]);
clear testFigure2
testFigure2(1,1)= gramm('x',plotDataY,'y',plotData);
testFigure2(1,1).stat_glm();
testFigure2(1,1).geom_point();
testFigure2.facet_grid([],clusterData);
testFigure2.draw();


%% PCA Analysis: Stop-Signal
for neuron_i = 1:size(acc_map_info,1)
    
    neuralFilename = acc_map_info.session{neuron_i};
    neuronLabel = acc_map_info.unit{neuron_i};    
    %... and find the corresponding behavior file index
    behFilename = data_findBehFile(neuralFilename);
    behaviorIdx = find(strcmp(dataFiles_beh,behFilename(1:end-4)));
    mid_ssd_idx = behavior(behaviorIdx).stopSignalBeh.midSSDidx;
 
    c_ssd1(neuron_i,:) = sdf_canceled_ssd{neuron_i}(mid_ssd_idx-1,:);
    c_ssd2(neuron_i,:) = sdf_canceled_ssd{neuron_i}(mid_ssd_idx,:);
    c_ssd3(neuron_i,:) = sdf_canceled_ssd{neuron_i}(mid_ssd_idx+1,:);
    
    ns_ssd1(neuron_i,:) = sdf_nostop_ssd{neuron_i}(mid_ssd_idx-1,:);
    ns_ssd2(neuron_i,:) = sdf_nostop_ssd{neuron_i}(mid_ssd_idx,:);
    ns_ssd3(neuron_i,:) = sdf_nostop_ssd{neuron_i}(mid_ssd_idx+1,:);   
end

% Set data parameters & windows
ephysWin = [-1000 2000]; winZero = abs(ephysWin(1));
plotWin = [-250:500]; 
analyseTime = [-100:200];
getColors

% Set PCA parameters
samplingRate = 1/1000; % inherent to the data. Do not change
numPCs = 8; % pick a number that will capture most of the variance
timeStep = 10; % smaller values will yield high resolution at the expense of computation time, default will sample at 20ms
withinConditionsOnly = false; % if true, will only analyze tanlging for times within the same condition

clear PCA_mainInput PCA_mainOutput
% Setup data for PCA analysis
PCA_mainInput(1).A = c_ssd1';
PCA_mainInput(1).times = plotWin';
PCA_mainInput(1).analyzeTimes = analyseTime';
PCA_mainInput(1).condition = 'SSD1 - canceled';
PCA_mainInput(1).color = [colors.canceled 0.25];

PCA_mainInput(2).A = c_ssd2';
PCA_mainInput(2).times = plotWin';
PCA_mainInput(2).analyzeTimes = analyseTime';
PCA_mainInput(2).condition = 'SSD2 - canceled';
PCA_mainInput(2).color = [colors.canceled 0.50];

PCA_mainInput(3).A = c_ssd3';
PCA_mainInput(3).times = plotWin';
PCA_mainInput(3).analyzeTimes = analyseTime';
PCA_mainInput(3).condition = 'SSD3 - canceled';
PCA_mainInput(3).color = [colors.canceled 1.00];

PCA_mainInput(4).A = ns_ssd1';
PCA_mainInput(4).times = plotWin';
PCA_mainInput(4).analyzeTimes = analyseTime';
PCA_mainInput(4).condition = 'SSD1 - canceled';
PCA_mainInput(4).color = [colors.nostop 0.25];

PCA_mainInput(5).A = ns_ssd2';
PCA_mainInput(5).times = plotWin';
PCA_mainInput(5).analyzeTimes = analyseTime';
PCA_mainInput(5).condition = 'SSD2 - canceled';
PCA_mainInput(5).color = [colors.nostop 0.50];

PCA_mainInput(6).A = ns_ssd3';
PCA_mainInput(6).times = plotWin';
PCA_mainInput(6).analyzeTimes = analyseTime';
PCA_mainInput(6).condition = 'SSD3 - canceled';
PCA_mainInput(6).color = [colors.nostop 1.00];
 

% Run PCA analysis
[~, PCA_mainOutput] = tangleAnalysis(PCA_mainInput, samplingRate,'numPCs',numPCs,'softenNorm',5 ,'timeStep',timeStep,'withinConditionsOnly',withinConditionsOnly); % soft normalize neural data
pca_figure(PCA_mainInput,PCA_mainOutput)

%% 