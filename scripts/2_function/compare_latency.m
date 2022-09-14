%{
Compare latency against other visual neuron activity in cortex and subcortex
2022-05-09 | S P Errington
%}

%% Organize: Load in processed data
load(fullfile(dirs.root,'results','mat_files','neuron_index.mat'))
load(fullfile(dirs.root,'results','mat_files','visual_info.mat'))
load(fullfile(dirs.root,'results','mat_files','visual_sdf_1.mat'))
load(fullfile(dirs.root,'data','heatmap_visual_colormap.mat'))

%% Organize: Read in pre-extracted latencies from csv file (Schmolesky et al., 1998)
latency_table = readtable(fullfile(dirs.data_proc,'visual_latencies','visual_latency_area.csv'));

%% Extract: Get CDF of latencies for each area in table
area_list = unique(latency_table.Area);

clear CDF_onset
for area_i = 1:length(area_list)
    area_obs_idx = [];
    area_obs_idx = find(strcmp(latency_table.Area,area_list{area_i}));
    CDF_onset.(area_list{area_i}(4:end)) = getCDF(latency_table.Latency(area_obs_idx));
end

narrow_onset_latency = visual_info.cd_detect_onset(neuron_index.vis_fac_narrow);
broad_onset_latency = visual_info.cd_detect_onset(neuron_index.vis_fac_broad);

CDF_onset.SEF_vis_narrow = getCDF(narrow_onset_latency);
CDF_onset.SEF_vis_broad = getCDF(broad_onset_latency);

%% Figure: Plot line CDF (p(active | time))
cdf_onset_areas = figure('Renderer', 'painters', 'Position', [100 100 500 300]);hold on
cdf_list = fieldnames(CDF_onset);

for cdf_i = 1:length(cdf_list)
    plot(CDF_onset.(cdf_list{cdf_i})(:,1),CDF_onset.(cdf_list{cdf_i})(:,2))
end

legend(cdf_list,'Location','EastOutside')

filename = fullfile(dirs.root,'results','gen_figures','cdf_onset_areas.pdf');
set(cdf_onset_areas,'PaperSize',[20 10]); %set the paper size to what you want
print(cdf_onset_areas,filename,'-dpdf') % then print it
close(cdf_onset_areas)

%% Figure: Plot line CDF (p(active | time))
Latency = [broad_onset_latency; narrow_onset_latency];
Area = [repmat({'12_SEF_broad'},length(broad_onset_latency),1);...
 repmat({'11_SEF_narrow'},length(narrow_onset_latency),1)];

sef_new_data = table(Latency,Area);
latency_table = [latency_table; sef_new_data];

clear figure_uOnset_area
figure_uOnset_area(1,1)=...
    gramm('x',latency_table.Area,...
    'y',latency_table.Latency);

figure_uOnset_area(1,1).stat_boxplot();
figure_uOnset_area(1,1).geom_jitter('alpha',0.2,'dodge',0.5);

figure_uOnset_area(1,1).no_legend;
figure_uOnset_area(1,1).set_names('x','Model','y','Latency (ms)');
figure_uOnset_area(1,1).axe_property('YLim',[0 200]);

figure_uOnset_area_pdf = figure('Position',[100 100 600 250]);
figure_uOnset_area.draw();

filename = fullfile(dirs.root,'results','gen_figures','boxplot_onset_areas.pdf');
set(figure_uOnset_area_pdf,'PaperSize',[20 10]); %set the paper size to what you want
print(figure_uOnset_area_pdf,filename,'-dpdf') % then print it
close(figure_uOnset_area_pdf)

%% Analysis: Compare latencies between monkeys
monkey_visual_latency_comp = ...
    table(neuron_index.visual_pos,...
    visual_info.cd_detect_onset(neuron_index.visual_pos),...
    visual_info.monkey_label(neuron_index.visual_pos),...
    'VariableNames',{'Neuron','Onset','Monkey'});

eu_onset = []; x_onset = [];
eu_onset = monkey_visual_latency_comp.Onset(strcmp(monkey_visual_latency_comp.Monkey,'Euler'));
x_onset = monkey_visual_latency_comp.Onset(strcmp(monkey_visual_latency_comp.Monkey,'Xena'));

[~,p,~,stat] = ttest2(eu_onset,x_onset);

fprintf('t(%.0f) = %.3f, p = %.3f    \n',stat.df, stat.tstat, p);

fprintf('Mean (+/- SEM) onset (Euler): %.2f +/-  %.2f ms \n',...
    mean(eu_onset), sem(eu_onset))

fprintf('Mean (+/- SEM) onset (Xena): %.2f +/-  %.2f ms \n',...
    mean(x_onset), sem(x_onset))