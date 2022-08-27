neuron_index.visual_fac_perp = neuron_index.visual_pos(neuron_index.visual_pos > 282);

neuronlist = []; neuronlist = neuron_index.visual_fac_perp;

for neuronlist_i = 1:length(neuronlist)
    neuron_i = neuronlist(neuronlist_i);
    session = executiveBeh.neuronMatPosit(neuron_i,1);
    monkey = executiveBeh.nhpSessions.monkeyNameLabel(session);
    site = executiveBeh.bioInfo.sessionSite(session);
    spkwidth = executiveBeh.bioInfo.spkWidth(neuron_i)*25;
    depth_ch = executiveBeh.bioInfo.depthInfo(neuron_i);
    depth_um = depth_ch*150;
    mod_onset = visual_info.cd_detect_onset(neuron_i);
    mod_offset = visual_info.cd_detect_offset(neuron_i);
    
    timedepth_table (neuronlist_i,:) =...
        table(neuron_i,session,monkey,site,spkwidth,depth_ch,...
        depth_um,mod_onset,mod_offset);

end

load(fullfile(dirs.root,'data','samplingBias.mat'))
load(fullfile(dirs.root,'data','td_visual_colormap.mat'))
n_depths = 19;

td_plot_win = [-500:1000];
td_plot_zero = abs(td_plot_win(1));
td_plot_array = zeros(n_depths,length(td_plot_win));

for neuronlist_i = 1:length(neuronlist) 
    depth = timedepth_table.depth_ch(neuronlist_i);
    onset = timedepth_table.mod_onset(neuronlist_i);
    offset = timedepth_table.mod_offset(neuronlist_i);
    
    if offset > max(td_plot_win)
        offset = max(td_plot_win);
    end
    
    td_plot_array(depth,onset+td_plot_zero:offset+td_plot_zero) = ...
        td_plot_array(depth,onset+td_plot_zero:offset+td_plot_zero) + 1;
end


for depth_i = 1:n_depths
    td_plot_array(depth_i,:) = td_plot_array(depth_i,:)./samplingBias.countPerDepth(depth_i);
end

timedepth_fig = figure; hold on
imagesc('XData',td_plot_win,'YData',1:n_depths,'CData',td_plot_array)
scatter(timedepth_table.mod_onset(timedepth_table.spkwidth > 249),...
    timedepth_table.depth_ch(timedepth_table.spkwidth > 249),'^','filled','MarkerEdgeColor','black','MarkerFacecolor','white')
scatter(timedepth_table.mod_onset(timedepth_table.spkwidth < 249),...
    timedepth_table.depth_ch(timedepth_table.spkwidth < 249),60,'h','filled','MarkerEdgeColor','black','MarkerFacecolor','white')
set(gca,'XLim',[-200 500],'YLim',[0.5 n_depths+0.5],...
    'YDir','Reverse')
xlabel('Time from Target (ms)'); ylabel('Depth');
colorbar
colormap(td_visual_colormap)
vline(0,'k'); hline(8.5,'k--')



filename = fullfile(dirs.root,'results','gen_figures','laminar_td_visual.pdf');
set(timedepth_fig,'PaperSize',[20 10]); %set the paper size to what you want
print(timedepth_fig,filename,'-dpdf') % then print it
close(timedepth_fig)

