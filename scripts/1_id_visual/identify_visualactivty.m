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



%%

figure;
plot(timewins.sdf, sdf_visual_zscore(1,:))
hline(0,'r-'); hline([-2 2],'r--'); hline([-6 6],'r--');
vline([min(timewins.baseline_target) max(timewins.baseline_target)],'k--')

%%
