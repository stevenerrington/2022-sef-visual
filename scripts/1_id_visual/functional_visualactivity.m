%{
Functional properties of visually responsive neurons.
2022-08-26 | S P Errington

In this script, we will X, Y, and Z

Dependencies: identify_visualactivity (neuron_index and visual_info
variables).
%}
%% Organize: Load in processed data
load(fullfile(dirs.root,'results','mat_files','neuron_index.mat'))
load(fullfile(dirs.root,'results','mat_files','visual_info.mat'))
load(fullfile(dirs.root,'results','mat_files','visual_sdf_1.mat'))
load(fullfile(dirs.root,'data','heatmap_visual_colormap.mat'))

%% Figure: Heatmap comparing visual and non-visual neurons

neuron_index.visual_pos = neuron_index.visual(visual_info.visual_fac(neuron_index.visual) == 1);
neuron_index.visual_neg = neuron_index.visual(visual_info.visual_sup(neuron_index.visual) == 1);
[~, latency_order_pos] = sort(visual_info.cd_detect_onset(neuron_index.visual_pos),'ascend');
[~, latency_order_neg] = sort(visual_info.cd_detect_onset(neuron_index.visual_neg),'descend');

neuron_order = [neuron_index.visual_pos(latency_order_pos); neuron_index.visual_neg(latency_order_neg)];

fig_out = figure('Renderer', 'painters', 'Position', [100 100 400 600]);
subplot(2,1,1)
imagesc('XData',timewins.sdf,'YData',1:length(neuron_index.visual), 'CData',sdf_visual_maxnorm(neuron_order,:))
xlim([-200 500]); ylim([1 400]);  caxis([-1,1]); colormap(heatmap_cmap); vline(0,'k')
xlabel('Time from Target (ms)'); ylabel('Neuron')

subplot(2,1,2)
imagesc('XData',timewins.sdf,'YData',1:length(neuron_index.nonvisual), 'CData',sdf_visual_maxnorm(neuron_index.nonvisual,:))
xlim([-200 500]); ylim([1 400]); caxis([-1,1]); colormap(heatmap_cmap); vline(0,'k')
xlabel('Time from Target (ms)'); 

% Once we're done with a page, save it and close it.
filename = fullfile(dirs.root,'results','gen_figures','population_figure_heatmap.pdf');
set(fig_out,'PaperSize',[20 10]); %set the paper size to what you want
print(fig_out,filename,'-dpdf') % then print it
close(fig_out)

%% Figure: Population spike-density function for visual neurons
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

% Produce the figure, collapsed across all monkeys
population_sdf(1,1)=gramm('x',timewins.sdf,'y',input_sdf,'color',class_label);
population_sdf(1,1).stat_summary();
population_sdf(1,1).axe_property('XLim',[-200 500],'YLim',[-15 15]);
population_sdf(1,1).set_names('x','Time from Target (ms)','y','FR (Z-score)');

% and create a subplot, split by monkey for transparency
population_sdf(2,1)=gramm('x',timewins.sdf,'y',input_sdf,'color',class_label);
population_sdf(2,1).stat_summary();
population_sdf(2,1).axe_property('XLim',[-200 500],'YLim',[-15 15]);
population_sdf(2,1).set_names('x','Time from Target (ms)','y','FR (Z-score)');
population_sdf(2,1).facet_grid([],visual_info.monkey_label); % <- performs monkey split

pop_figure = figure('Renderer', 'painters', 'Position', [100 100 600 600]);
population_sdf.draw();

% Once we're done with a page, save it and close it.
filename = fullfile(dirs.root,'results','gen_figures','population_figure_visual.pdf');
set(pop_figure,'PaperSize',[20 10]); %set the paper size to what you want
print(pop_figure,filename,'-dpdf') % then print it
close(pop_figure)

%% Figure: Proportion of neurons that are faciliated, suppressed, or unmodulated.
spkwidth_cutoff = 249;
neuron_index.vis_fac_narrow = visual_info.neuron_i(ismember(1:575,neuron_index.visual)' & visual_info.visual_fac == 1 & visual_info.spkWidth < spkwidth_cutoff);
neuron_index.vis_fac_broad =  visual_info.neuron_i(ismember(1:575,neuron_index.visual)' & visual_info.visual_fac == 1 & visual_info.spkWidth > spkwidth_cutoff);
neuron_index.vis_sup_narrow = visual_info.neuron_i(ismember(1:575,neuron_index.visual)' & visual_info.visual_sup == 1 & visual_info.spkWidth < spkwidth_cutoff);
neuron_index.vis_sup_broad =  visual_info.neuron_i(ismember(1:575,neuron_index.visual)' & visual_info.visual_sup == 1 & visual_info.spkWidth > spkwidth_cutoff);
neuron_index.nonvis_narrow = visual_info.neuron_i(ismember(1:575,neuron_index.nonvisual)' & visual_info.spkWidth > spkwidth_cutoff);
neuron_index.nonvis_broad = visual_info.neuron_i(ismember(1:575,neuron_index.nonvisual)' & visual_info.spkWidth <= spkwidth_cutoff);

n_visual_fac_narrow = length(neuron_index.vis_fac_narrow);
n_visual_fac_broad = length(neuron_index.vis_fac_broad);
n_visual_sup_narrow = length(neuron_index.vis_sup_narrow);
n_visual_sup_broad = length(neuron_index.vis_sup_broad);
n_nonvisual_narrow = length(neuron_index.nonvis_narrow);
n_nonvisual_broad = length(neuron_index.nonvis_broad);

n_total = sum([n_visual_fac_narrow, n_visual_fac_broad, n_visual_sup_narrow, n_visual_sup_broad, n_nonvisual_narrow, n_nonvisual_broad]);

clear pie_data pie_labels
pie_data = [n_visual_fac_narrow, n_visual_fac_broad, n_visual_sup_narrow, n_visual_sup_broad, n_nonvisual_narrow, n_nonvisual_broad];
pie_labels = {'Visual (Facilitated, Narrow)','Visual (Facilitated, Broad)','Visual (Suppressed, Narrow)','Visual (Suppressed, Broad)','Nonvisual (Narrow)','Nonvisual (Broad)'};

pie_figure = figure('Renderer', 'painters', 'Position', [100 100 600 300]);
pie(pie_data, [1 1 1 1 1 1], pie_labels);

filename = fullfile(dirs.root,'results','gen_figures','population_pie_visual.pdf');
set(pie_figure,'PaperSize',[20 10]); %set the paper size to what you want
print(pie_figure,filename,'-dpdf') % then print it
close(pie_figure)

%% Figure: Cumulative distribution function of onset times.

CDF_onset.vis_fac_narrow = getCDF(visual_info.cd_detect_onset(neuron_index.vis_fac_narrow));
CDF_onset.vis_fac_broad = getCDF(visual_info.cd_detect_onset(neuron_index.vis_fac_broad));
CDF_onset.vis_sup_narrow = getCDF(visual_info.cd_detect_onset(neuron_index.vis_sup_narrow));
CDF_onset.vis_sup_broad = getCDF(visual_info.cd_detect_onset(neuron_index.vis_sup_broad));
CDF_onset.vis_all_narrow = getCDF(visual_info.cd_detect_onset([neuron_index.vis_fac_narrow;neuron_index.vis_sup_narrow]));
CDF_onset.vis_all_broad = getCDF(visual_info.cd_detect_onset([neuron_index.vis_fac_broad;neuron_index.vis_sup_broad]));

clear cdf_fig_data cdf_fig_label cdf_fig_label_spikewidth
cdf_fig_data = [CDF_onset.vis_all_narrow; CDF_onset.vis_all_broad];
cdf_fig_label = [repmat({'narrow'},length(CDF_onset.vis_all_narrow),1);...
    repmat({'broad'},length(CDF_onset.vis_all_broad),1)];

clear cdf_visonset_fig
cdf_visonset_fig(1,1)=gramm('x',cdf_fig_data(:,1),'y',cdf_fig_data(:,2),'color',cdf_fig_label);
cdf_visonset_fig(1,1).geom_line();
cdf_fig = figure('Position',[100 100 400 300]);
cdf_visonset_fig.draw();

filename = fullfile(dirs.root,'results','gen_figures','cdf_visual_onset.pdf');
set(cdf_fig,'PaperSize',[20 10]); %set the paper size to what you want
print(cdf_fig,filename,'-dpdf') % then print it
close(cdf_fig)

%% Analysis: Compare onset time between neuron classes


%% Figure: Example neuron
neuron_index.example_vis_fac = 71;
neuron_index.example_vis_sup = 13;

neuron_plot_list = [neuron_index.example_vis_fac, neuron_index.example_vis_sup];
ylim_list = [60 60];

for neuron_loop_i = 1:length(neuron_plot_list)
    neuron_i = neuron_plot_list(neuron_loop_i);
    fprintf('Extracting target aligned SDF for neuron %i of %i. \n',neuron_i, 575)
    
    session_i = executiveBeh.neuronMatPosit(neuron_i,1);
    
    spk_in = load(fullfile(dirs.data_spike,['SDF_' int2str(neuron_i)]));
    
    trial_in = []; trial_in = executiveBeh.ttx.GO{session_i}; trial_n = length(trial_in);
    baseline_win = []; baseline_win = timewins.baseline_target + timewins.zero;
    target_win = []; target_win = timewins.visual_target + timewins.zero;
    

    clear example_sdf_visual_raw example_sdf_visual_ylim
    for trial_loop_i = 1:trial_n
        trial_i = trial_in(trial_loop_i);
        example_sdf_visual_raw{trial_loop_i,:} = spk_in.SDF.target(trial_i,:);
    end
    
    % Produce the figure, collapsed across all monkeys
    example_sdf(neuron_loop_i,1)=gramm('x',timewins.sdf,'y',example_sdf_visual_raw);
    example_sdf(neuron_loop_i,1).stat_summary();
    example_sdf(neuron_loop_i,1).axe_property('XLim',[-200 500],...
        'YLim',[20 ylim_list(neuron_loop_i)]);
    example_sdf(neuron_loop_i,1).geom_vline('XIntercept',visual_info.cd_detect_onset(neuron_i));
    example_sdf(neuron_loop_i,1).set_names('x','Time from Target (ms)','y','FR (spk/sec)');
    
end
example_sdf_out = figure('Renderer', 'painters', 'Position', [100 100 300 300]);
example_sdf.draw();

% Once we're done with a page, save it and close it.
filename = fullfile(dirs.root,'results','gen_figures','exampleSDF_figure_visual.pdf');
set(example_sdf_out,'PaperSize',[20 10]); %set the paper size to what you want
print(example_sdf_out,filename,'-dpdf') % then print it
close(example_sdf_out)

% Find the relevant session labels
session_label = FileNames{executiveBeh.neuronMatPosit(neuron_plot_list,1)};


%% Organize: Save updated visual neuron categories.
save(fullfile(dirs.root,'results','mat_files','neuron_index.mat'),'neuron_index')





