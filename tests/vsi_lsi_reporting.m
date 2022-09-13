sum(vmi_index(neuron_index.visual_pos) > 0.4 )

sum(visual_lateral.lateral_flag_onset == 1 & visual_lateral.LSI_onset < 0)
sum(visual_lateral.lateral_flag_onset == 1 & visual_lateral.LSI_onset > 0)
sum(visual_lateral.lateral_flag_onset == 0)

mean(visual_lateral.LSI_onset)


sum(visual_value.value_flag_onset == 1 & visual_value.VSI_onset < 0)
sum(visual_value.value_flag_onset == 1 & visual_value.VSI_onset > 0)
sum(visual_value.value_flag_onset == 0)

mean(visual_value.VSI_onset)





[~,p,~,stat] = ttest(visual_lateral.LSI_onset);

fprintf('t(%.0f) = %.3f, p = %.3f    \n',stat.df, stat.tstat, p)

[~,p,~,stat] = ttest(visual_value.VSI_onset);

fprintf('t(%.0f) = %.3f, p = %.3f    \n',stat.df, stat.tstat, p)