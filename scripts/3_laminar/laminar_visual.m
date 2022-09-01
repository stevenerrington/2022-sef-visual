neuron_index.visual_fac_perp = neuron_index.visual_pos(neuron_index.visual_pos > 282);

neuronlist = []; neuronlist = neuron_index.visual_fac_perp;
clear timedepth_table
for neuronlist_i = 1:length(neuronlist)
    neuron_i = neuronlist(neuronlist_i);
    session = executiveBeh.neuronMatPosit(neuron_i,1);
    monkey = executiveBeh.nhpSessions.monkeyNameLabel(session);
    site = executiveBeh.bioInfo.sessionSite(session);
    spkwidth = executiveBeh.bioInfo.spkWidth(neuron_i)*25;
    depth_ch = executiveBeh.bioInfo.depthInfo(neuron_i);
    depth_um = depth_ch*150;
    mod_onset_target = visual_info.cd_detect_onset(neuron_i);
    mod_offset_target = visual_info.cd_detect_offset(neuron_i);
    mod_onset_saccade = visual_info.saccade_onset(neuron_i);
    mod_offset_saccade = visual_info.saccade_offset(neuron_i); 
    
    timedepth_table (neuronlist_i,:) =...
        table(neuron_i,session,monkey,site,spkwidth,depth_ch,...
        depth_um,mod_onset_target,mod_offset_target,mod_onset_saccade,mod_offset_saccade);

end

load(fullfile(dirs.root,'data','samplingBias.mat'))
load(fullfile(dirs.root,'data','td_visual_colormap.mat'))
n_depths = 19;

td_plot_win_target = [-500:1000]; td_plot_zero_target = abs(td_plot_win_target(1));
td_plot_array_target = zeros(n_depths,length(td_plot_win_target));

td_plot_win_saccade = [-1000:2000]; td_plot_zero_saccade = abs(td_plot_win_saccade(1));
td_plot_array_saccade = zeros(n_depths,length(td_plot_win_saccade));

for neuronlist_i = 1:length(neuronlist) 
    depth = timedepth_table.depth_ch(neuronlist_i);
    onset_target = timedepth_table.mod_onset_target(neuronlist_i);
    offset_target = timedepth_table.mod_offset_target(neuronlist_i);
    onset_saccade = timedepth_table.mod_onset_saccade(neuronlist_i);
    offset_saccade = timedepth_table.mod_offset_saccade(neuronlist_i);    
    
    if offset_target > max(td_plot_win_target)
        offset_target = max(td_plot_win_target);
    end
    if offset_saccade > max(td_plot_win_saccade)
        offset_saccade = max(td_plot_win_saccade);
    end
    
    td_plot_array_target(depth,onset_target+td_plot_zero_target:offset_target+td_plot_zero_target) = ...
        td_plot_array_target(depth,onset_target+td_plot_zero_target:offset_target+td_plot_zero_target) + 1;
    
    if ~isnan(onset_saccade)
    td_plot_array_saccade(depth,onset_saccade+td_plot_zero_saccade:offset_saccade+td_plot_zero_saccade) = ...
        td_plot_array_saccade(depth,onset_saccade+td_plot_zero_saccade:offset_saccade+td_plot_zero_saccade) + 1;

    end
end

for depth_i = 1:n_depths
    td_plot_array_target(depth_i,:) = td_plot_array_target(depth_i,:)./samplingBias.countPerDepth(depth_i);
    td_plot_array_saccade(depth_i,:) = td_plot_array_saccade(depth_i,:)./samplingBias.countPerDepth(depth_i);
end



timedepth_fig = figure('Renderer', 'painters', 'Position', [100 100 1000 300]);
subplot(1,3,1);hold on
imagesc('XData',td_plot_win_target,'YData',1:n_depths,'CData',td_plot_array_target)
scatter(timedepth_table.mod_onset_target(timedepth_table.spkwidth > 249),...
    timedepth_table.depth_ch(timedepth_table.spkwidth > 249),'^','filled','MarkerEdgeColor','black','MarkerFacecolor','white')
scatter(timedepth_table.mod_onset_target(timedepth_table.spkwidth < 249),...
    timedepth_table.depth_ch(timedepth_table.spkwidth < 249),60,'h','filled','MarkerEdgeColor','black','MarkerFacecolor','white')
set(gca,'XLim',[-200 500],'YLim',[0.5 n_depths+0.5],...
    'YDir','Reverse')
xlabel('Time from Target (ms)'); ylabel('Depth');
colormap(td_visual_colormap)
vline(0,'k'); hline(8.5,'k--')


subplot(1,3,[2 3]);hold on
imagesc('XData',td_plot_win_saccade,'YData',1:n_depths,'CData',td_plot_array_saccade)
% scatter(timedepth_table.mod_onset_saccade(timedepth_table.spkwidth > 249),...
%     timedepth_table.depth_ch(timedepth_table.spkwidth > 249),'^','filled','MarkerEdgeColor','black','MarkerFacecolor','white')
% scatter(timedepth_table.mod_onset_saccade(timedepth_table.spkwidth < 249),...
%     timedepth_table.depth_ch(timedepth_table.spkwidth < 249),60,'h','filled','MarkerEdgeColor','black','MarkerFacecolor','white')
set(gca,'XLim',[-200 1200],'YLim',[0.5 n_depths+0.5],...
    'YDir','Reverse')
xlabel('Time from Saccade (ms)'); ylabel('Depth');
colorbar
colormap(td_visual_colormap)
vline(0,'k'); vline(600,'k:'); hline(8.5,'k--');



filename = fullfile(dirs.root,'results','gen_figures','laminar_td_visual.pdf');
set(timedepth_fig,'PaperSize',[20 10]); %set the paper size to what you want
print(timedepth_fig,filename,'-dpdf') % then print it
close(timedepth_fig)





%% TEST:

visual_long_neurons = timedepth_table.neuron_i(timedepth_table.mod_offset_saccade > 500);

clear input_sdf* population_sdf
input_sdf_target = num2cell(sdf_target_ns_zscore(visual_long_neurons,:), 2);
input_sdf_saccade = num2cell(sdf_saccade_ns_zscore(visual_long_neurons,:), 2);
input_sdf_tone = num2cell(sdf_tone_ns_zscore(visual_long_neurons,:), 2);

% Produce the figure, collapsed across all monkeys
population_sdf(1,1)=gramm('x',timewins.sdf,'y',input_sdf_target);
population_sdf(1,2)=gramm('x',timewins.sdf,'y',input_sdf_saccade);
population_sdf(1,3)=gramm('x',timewins.sdf,'y',input_sdf_tone);

population_sdf(1,1).stat_summary(); population_sdf(1,2).stat_summary(); population_sdf(1,3).stat_summary();
population_sdf(1,1).axe_property('XLim',[-200 800],'YLim',[-2 35]);
population_sdf(1,2).axe_property('XLim',[-200 200],'YLim',[-2 35]);
population_sdf(1,3).axe_property('XLim',[-400 600],'YLim',[-2 35]);
population_sdf(1,1).set_names('x','Time from Target (ms)','y','FR (Z-score)');
population_sdf(1,2).set_names('x','Time from Saccade (ms)','y','FR (Z-score)');
population_sdf(1,3).set_names('x','Time from Tone (ms)','y','FR (Z-score)');

pop_figure = figure('Renderer', 'painters', 'Position', [100 100 1200 300]);
population_sdf.draw();

% Once we're done with a page, save it and close it.
filename = fullfile(dirs.root,'results','gen_figures','population_figure_movement.pdf');
set(pop_figure,'PaperSize',[20 10]); %set the paper size to what you want
print(pop_figure,filename,'-dpdf') % then print it
close(pop_figure)