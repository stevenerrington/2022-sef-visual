
%% Setup: Load pre-processed data
% Load in pre-processed data
load(fullfile(dirs.root,'src','pran_code','FileNames.mat'));
load(fullfile(dirs.root,'src','pran_code','executiveBeh.mat'));
load(fullfile(dirs.root,'src','pran_code','ChannelDepthMap.mat'));
load(fullfile(dirs.root,'src','pran_code','piadepth_gray.mat'));

%% Setup: Define parameters
Windowtgt=[-200 200];
Windowbl= [-200 0]; %baseline window includes pretarget ramping?

%% Extract: Load in LFP data and extract CSD
for perp_session_i=1:16 %perpendicular recording sessions 1:6 euler; 7:16 xena
    fprintf('Analysing session %s  | (%i)    \n',FileNames{perp_session_i+13}, perp_session_i)
    
    % Load in relevant data files
    data_in = [];
    data_in = load(fullfile(dirs.data_raw,'rawData', FileNames{perp_session_i+13}));
    
    clear LFPdata_raw DepthAlignment LFPdata_depth
    LFPdata_raw = [data_in.AD17 data_in.AD18 data_in.AD19 data_in.AD20 data_in.AD21 data_in.AD22...
        data_in.AD23 data_in.AD24 data_in.AD25 data_in.AD26 data_in.AD27 data_in.AD28...
        data_in.AD29 data_in.AD30 data_in.AD31 data_in.AD32 data_in.AD33 data_in.AD34...
        data_in.AD35 data_in.AD36 data_in.AD37 data_in.AD38 data_in.AD39 data_in.AD40]/1000;
    % Get depth alignment 
    DepthAlignment = [min(ChannelDepthMap.GrayMatterMapping{perp_session_i+13,1}) max(ChannelDepthMap.GrayMatterMapping{perp_session_i+13,1})];
    LFPdata_depth = LFPdata_raw(:,DepthAlignment(1):DepthAlignment(2));
    
    % I considered a filter here?
    %fLFP = SEF_LFP_Filter(inputLFP ,lowerBand ,upperBand ,samplingFreq)
    
    
    for trial_type_i=1:2

        % Get alignment times
        TgtAlign_in = [];
        TgtAlign_in = data_in.Target_ + data_in.TrialStart_;
        
        TgtAlign = [];
        if trial_type_i==1
            TgtAlign = TgtAlign_in(trial_index_matched.lateral{perp_session_i+13}.left,:);
        else
            TgtAlign = TgtAlign_in(trial_index_matched.lateral{perp_session_i+13}.right,:);
        end
        
        % Align LFPs on defined windows and trials
        LFPtgtArray = [];
        for i=1:length(TgtAlign)
            LFPtgtarray = []; LFPblarray = []; 
            LFPtgtarray = LFPdata_depth(TgtAlign(i)+Windowtgt(1):TgtAlign(i)+Windowtgt(2)-1,:);
            LFPblarray = LFPdata_depth(TgtAlign(i)+Windowbl(1):TgtAlign(i)+Windowbl(2)-1,:);
            LFPtgtArray(:,:,i)=LFPtgtarray-mean(LFPblarray);
        end
        
        avgtgtLFP = [];
        avgtgtLFP=transpose(mean(LFPtgtArray,3));
        
        %% Calculate CSD on LFP array
        % Define parameters
        Ne= size(avgtgtLFP,1);
        
        a = 0.1; % [mm] position of the first electrode, shifts sources down, assumed 1 micron
        elec_spacing = 0.15; % [mm] electrode spacing
        ze = a:elec_spacing:((Ne-1)*elec_spacing + a); % electrode positions with respect to the pia surface
        % WARNING: position of the first electrode must be different from zero
        el_pos = ze*1e-3;  % mm to m
        cond = 0.4; %[S/m] gray matter conductance | NOTE: if length(cond)==1, the
        % function iCSDspline_Pettersen() considers con_top = cond (conductance at
        % the top of the cylinder or above the pia matter) %papers with macaques,
        % 0.4
        gauss_sigma = 0.1e-3;   %[m] Gaussian filter std, how smooth you want LFP x depth, increase more smooth
        diam = 3e-3; % [m] cylinder diameter
        Ve = avgtgtLFP;
        
        % Ve == LOCAL FIELD POTENTIALS IN VOLTS
        % Outputs: zs ~ mm & iCSD ~ nA/mm3
        
        % Run CSD
        clear zs iCSD nan_pad_csd
        [zs, iCSD] = iCSDspline_Pettersen(Ve,el_pos,diam,cond,gauss_sigma);
        % iCSD will be a 10 x nCh array
        
        nan_pad_size = (24-(size(iCSD,1)/10))*10;
        nan_pad_csd = nan(nan_pad_size,size(iCSD,2));
        
        iCSD_out = vertcat(iCSD,nan_pad_csd);
        tspan = -200:200;
        
        if trial_type_i==1
            perpCSD.Left.session(perp_session_i) = perp_session_i+13;
            perpCSD.Left.CSD(:,:,perp_session_i) = iCSD_out;
        else
            perpCSD.Right.session(perp_session_i) = perp_session_i+13;
            perpCSD.Right.CSD(:,:,perp_session_i) = iCSD_out;
        end
        
    end
end

%% Figures: Plot CSD for ipsi, contra, and difference
load('csd_data')
load(fullfile(dirs.root,'src','pran_code','FileNames.mat'));
load(fullfile(dirs.root,'src','pran_code','executiveBeh.mat'));
load(fullfile(dirs.root,'src','pran_code','ChannelDepthMap.mat'));
load(fullfile(dirs.root,'src','pran_code','piadepth_gray.mat'));

% Parameters --------------------------------------------------
params.plot.time = -200:200;
params.plot.gauss_smooth = 10;
params.plot.color_csd_max = 75;

% Data --------------------------------------------------
figure_data.csd.left = perpCSD.Left.CSD(1:160,:,:);
figure_data.csd.right = perpCSD.Right.CSD(1:160,:,:);

for perp_session_i = 1:16
    
    figure_data.csd.diff(:,:,perp_session_i)=...
        figure_data.csd.left(:,:,perp_session_i)-...
        figure_data.csd.right(:,:,perp_session_i);
    
end

% Generate --------------------------------------------------
csd_figure = figure('Renderer', 'painters', 'Position', [100 100 950 300]);
ax1 = subplot(1,3,1);
imagesc(params.plot.time, piadepth_gray, smoothdata(nanmean(figure_data.csd.left,3),'gaussian',params.plot.gauss_smooth));
ax1_c = colorbar; colormap(flipud(bluewhitered)); caxis([-params.plot.color_csd_max params.plot.color_csd_max]);
vline(0,'k'); hline(median(piadepth_gray),'k--'), xlim([-50 200])

ax2 = subplot(1,3,2);
imagesc(params.plot.time, piadepth_gray, smoothdata(nanmean(figure_data.csd.right,3),'gaussian',params.plot.gauss_smooth));
ax2_c = colorbar; colormap(flipud(bluewhitered)); caxis([-params.plot.color_csd_max params.plot.color_csd_max]);
vline(0,'k'); hline(median(piadepth_gray),'k--'), xlim([-50 200])

ax3 = subplot(1,3,3);
imagesc(params.plot.time, piadepth_gray, smoothdata(nanmean(figure_data.csd.diff,3),'gaussian',params.plot.gauss_smooth));
ax3_c = colorbar; colormap(flipud(bluewhitered)); caxis([-35 35]);
vline(0,'k'); hline(median(piadepth_gray),'k--'), xlim([-50 200])

% Format --------------------------------------------------
set([ax1 ax2 ax3],'FontSize',9,'LineWidth',1,...
    'xtick',[-200 -100 0 100 200],...
    'xticklabels',{'-200','-100', '0',  '100', '200'},...
    'ytick',[],...
    'xlim',[-50 200],'ylim', [0 max(piadepth_gray)])

ax1_c.FontSize = 9;
ax1_c.Ticks = [-params.plot.color_csd_max 0 params.plot.color_csd_max];
ax1_c.LineWidth = 1;

ax2_c.FontSize = 9;
ax2_c.Ticks = [-params.plot.color_csd_max 0 params.plot.color_csd_max];
ax2_c.LineWidth = 1;

ax3_c.FontSize = 9;
ax3_c.Ticks = [-35 0 35];
ax3_c.LineWidth = 1;

% Once we're done with a page, save it and close it.
filename = fullfile(dirs.root,'results','gen_figures','csd_figure.pdf');
set(csd_figure,'PaperSize',[20 10]); %set the paper size to what you want
print(csd_figure,filename,'-dpdf') % then print it
close(csd_figure)


%% Analysis
csd_mean_analysis_window = 200+[100:200];
csd_mean_baseline_window = 200+[-150:50];

% Extract mean CSD values in a given window
for perp_session_i = 1:16
    
    csd_mean_analysis.baseline.l2(perp_session_i,1) = nanmean(nanmean(figure_data.csd.diff([1:40],csd_mean_baseline_window,perp_session_i)));
    csd_mean_analysis.baseline.l3(perp_session_i,1) = nanmean(nanmean(figure_data.csd.diff([41:80],csd_mean_baseline_window,perp_session_i)));
    csd_mean_analysis.baseline.l5(perp_session_i,1) = nanmean(nanmean(figure_data.csd.diff([81:120],csd_mean_baseline_window,perp_session_i)));
    csd_mean_analysis.baseline.l6(perp_session_i,1) = nanmean(nanmean(figure_data.csd.diff([121:160],csd_mean_baseline_window,perp_session_i)));
    
    csd_mean_analysis.visual.l2(perp_session_i,1) = nanmean(nanmean(figure_data.csd.diff([1:40],csd_mean_analysis_window,perp_session_i)));
    csd_mean_analysis.visual.l3(perp_session_i,1) = nanmean(nanmean(figure_data.csd.diff([41:80],csd_mean_analysis_window,perp_session_i)));
    csd_mean_analysis.visual.l5(perp_session_i,1) = nanmean(nanmean(figure_data.csd.diff([81:120],csd_mean_analysis_window,perp_session_i)));
    csd_mean_analysis.visual.l6(perp_session_i,1) = nanmean(nanmean(figure_data.csd.diff([121:160],csd_mean_analysis_window,perp_session_i)));
    
end

csd_mean_baseline_datain = [csd_mean_analysis.baseline.l2; csd_mean_analysis.baseline.l3; csd_mean_analysis.baseline.l5; csd_mean_analysis.baseline.l6];
csd_mean_baseline_labels = [repmat({'1_L2'},16,1);repmat({'2_L3'},16,1);repmat({'3_L5'},16,1);repmat({'4_L6'},16,1)];
csd_mean_visual_datain = [csd_mean_analysis.visual.l2; csd_mean_analysis.visual.l3; csd_mean_analysis.visual.l5; csd_mean_analysis.visual.l6];
csd_mean_visual_labels = [repmat({'1_L2'},16,1);repmat({'2_L3'},16,1);repmat({'3_L5'},16,1);repmat({'4_L6'},16,1)];

csd_epoch_label = [repmat({'1_baseline'},length(csd_mean_baseline_datain),1);...
    repmat({'2_visual'},length(csd_mean_visual_datain),1)];

% Produce the figure, collapsed across all monkeys
clear test
test(1,1)=gramm('x',[csd_mean_baseline_labels;csd_mean_visual_labels],...
    'y',[csd_mean_baseline_datain;csd_mean_visual_datain],...
    'color',csd_epoch_label);
test(1,1).stat_summary('geom',{'bar','black_errorbar'},'dodge',1,'width',0.5);
test(1,1).geom_jitter('alpha',0.2,'width',0.0,'dodge',1);
test(1,1).geom_hline('yintercept',0,'style','k--');
test(1,1).set_names('x','Cortical Layer','y','Mean CSD (uA)');
test(1,1).no_legend;

csd_bar_figure=figure('Renderer', 'painters', 'Position', [100 100 300 300]);
test.draw

% Once we're done with a page, save it and close it.
filename = fullfile(dirs.root,'results','gen_figures','csd_bar_figure.pdf');
set(csd_bar_figure,'PaperSize',[20 10]); %set the paper size to what you want
print(csd_bar_figure,filename,'-dpdf') % then print it
close(csd_bar_figure)




%%

data_in = []; data_in = nanmean(figure_data.csd.diff(:,:,:),3);

figure; hold on
ax_0 = subplot(1,1,1);
imagesc(params.plot.time, piadepth_gray, smoothdata(nanmean(figure_data.csd.diff,3),'gaussian',params.plot.gauss_smooth));
ax_0_c = colorbar; colormap(flipud(bluewhitered)); caxis([-35 35]);
vline(0,'k'); hline(median(piadepth_gray),'k--'), xlim([-50 200])
set([ax_0],'FontSize',9,'LineWidth',1,...
    'xtick',[-200 -100 0 100 200],...
    'xticklabels',{'-200','-100', '0',  '100', '200'},...
    'ytick',[],...
    'xlim',[-50 200],'ylim', [0 max(piadepth_gray)],...
    'YDir','reverse')


[patchEdges, patch] = patchEdgeDetector_pipeline(data_in, [50 2], 'strict');

figure; hold on
for patch_i = 1:length(patch)
    plot(patchEdges{patch_i}(1,:),patchEdges{patch_i}(2,:),'k-')
end

%% DETERMINE SIGNIFICANT PORTION OF CSD



test_datain = nanmean(figure_data.csd.diff,3);
figure;
imagesc(test_datain)

[allPeakFreq,allPeakTime] = find(imregionalmax(test_datain));
hold on 
scatter(allPeakTime,allPeakFreq,'ro')