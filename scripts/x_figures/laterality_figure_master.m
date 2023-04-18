
% visual_lateral.neuron_i(visual_lateral.LSI_onset == min(visual_lateral.LSI_onset))
% visual_lateral.LSI_onset(neuronlist)

example_neuron_i = 124; %124 for ipsi, 293 for contra
example_neuron_session = executiveBeh.neuronMatPosit(example_neuron_i,1);
example_neuron_filename = FileNames{example_neuron_session};

% Example neuron session eye trace
data_in = load(fullfile(dirs.data_raw,'rawData',example_neuron_filename),'EyeX_');
[trial_signal_target] = align_analog_trial(data_in.EyeX_, executiveBeh.TrialEventTimes_Overall{example_neuron_session}, 2, [-1000 2000]);
[trial_signal_saccade] = align_analog_trial(data_in.EyeX_, executiveBeh.TrialEventTimes_Overall{example_neuron_session}, 4, [-1000 2000]);

clear input_eye_trace_left input_eye_trace_right eye_trace_label eye_trace_times

input_eye_trace_left_target = trial_signal_target(executiveBeh.ttx.GO_Left{example_neuron_session},:);
input_eye_trace_right_target = trial_signal_target(executiveBeh.ttx.GO_Right{example_neuron_session},:);
input_eye_trace_left_saccade = trial_signal_saccade(executiveBeh.ttx.GO_Left{example_neuron_session},:);
input_eye_trace_right_saccade = trial_signal_saccade(executiveBeh.ttx.GO_Right{example_neuron_session},:);
eye_trace_label = [repmat({'1_Left'},size(input_eye_trace_left_target,1),1); repmat({'2_Right'},size(input_eye_trace_right_target,1),1)];
eye_trace_times = repmat([-999:2000],size(eye_trace_label,1),1);

% Example neuron raster & sdf
example_neuron_sdf = [];
data_in.spikes = load(fullfile(dirs.data_raw,'spikingData',['Spikes_' int2str(example_neuron_i)]));
data_in.sdf = load(fullfile(dirs.data_raw,'spikingData',['SDF_' int2str(example_neuron_i)]));

clear spk sdf
for trl_i = 1:size(executiveBeh.TrialEventTimes_Overall{example_neuron_session},1)
    spk.target{trl_i,1} = find(data_in.spikes.Spikes.target(trl_i,:) == 1)-1000;
    sdf.target{trl_i,1} = data_in.sdf.SDF.target(trl_i,:);
    spk.saccade{trl_i,1} = find(data_in.spikes.Spikes.saccade(trl_i,:) == 1)-1000;
    sdf.saccade{trl_i,1} = data_in.sdf.SDF.saccade(trl_i,:);
end

% Population SDF.
neuron_idx = [];
neuron_idx = find(visual_lateral.lateral_flag_onset == 1 &...
    visual_lateral.LSI_onset < 0);

clear input_sdf_left* input_sdf_right*
input_sdf_left_target = num2cell(lateral_sdf_left_zscore(neuron_idx,:), 2);
input_sdf_right_target = num2cell(lateral_sdf_right_zscore(neuron_idx,:), 2);
input_sdf_left_saccade = num2cell(lateral_sdf_left_zscore_saccade(neuron_idx,:), 2);
input_sdf_right_saccade = num2cell(lateral_sdf_right_zscore_saccade(neuron_idx,:), 2);

labels_value = [repmat({'1_left'},length(input_sdf_left_target),1);repmat({'2_right'},length(input_sdf_right_target),1)];


%% Create figure
clear laterality_figure

% Input data >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
% Example session eye position
laterality_figure(1,1)=gramm('x',num2cell(eye_trace_times,2),'y',num2cell([input_eye_trace_left_target;input_eye_trace_right_target],2),'color',eye_trace_label);
laterality_figure(1,2)=gramm('x',num2cell(eye_trace_times,2),'y',num2cell([input_eye_trace_left_saccade;input_eye_trace_right_saccade],2),'color',eye_trace_label);

% Example neuron raster
laterality_figure(2,1)=gramm('x',[spk.target(executiveBeh.ttx.GO_Left{example_neuron_session});spk.target(executiveBeh.ttx.GO_Right{example_neuron_session})],...
    'color',eye_trace_label);
laterality_figure(2,2)=gramm('x',[spk.saccade(executiveBeh.ttx.GO_Left{example_neuron_session});spk.saccade(executiveBeh.ttx.GO_Right{example_neuron_session})],...
    'color',eye_trace_label);

% Example neuron SDF
laterality_figure(3,1)=gramm('x',timewins.sdf, ...
    'y',[sdf.target(executiveBeh.ttx.GO_Left{example_neuron_session});sdf.target(executiveBeh.ttx.GO_Right{example_neuron_session})],...
    'color',eye_trace_label);
laterality_figure(3,2)=gramm('x',timewins.sdf, ...
    'y',[sdf.saccade(executiveBeh.ttx.GO_Left{example_neuron_session});sdf.saccade(executiveBeh.ttx.GO_Right{example_neuron_session})],...
    'color',eye_trace_label);

% Population neuron SDF
% Note: extracted in the lateral-visual script
laterality_figure(4,1)=gramm('x',timewins.sdf,...
    'y',[input_sdf_left_target;input_sdf_right_target],...
    'color',labels_value);
laterality_figure(4,2)=gramm('x',timewins.sdf,...
    'y',[input_sdf_left_saccade;input_sdf_right_saccade],...
    'color',labels_value);

% Define figure type >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
laterality_figure(1,1).geom_line('alpha',0.2); laterality_figure(1,2).geom_line('alpha',0.2);
laterality_figure(2,1).geom_raster; laterality_figure(2,2).geom_raster;
laterality_figure(3,1).stat_summary; laterality_figure(3,2).stat_summary;
laterality_figure(4,1).stat_summary; laterality_figure(4,2).stat_summary;

% Parameterize figures >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
xlim_target = [-50 200];
xlim_saccade = [-50 200];

laterality_figure(1,1).axe_property('XLim',xlim_target,'YLim',[-15 15],'XTick',[],'XColor',[1 1 1]);
laterality_figure(1,2).axe_property('XLim',xlim_saccade,'YLim',[-15 15],'XTick',[],'XColor',[1 1 1],'YTick',[],'YColor',[1 1 1]);
laterality_figure(1,1).set_names('X',[]); laterality_figure(1,1).set_names('Y','Eye X (deg)');
laterality_figure(1,2).set_names('X',[]);

laterality_figure(2,1).axe_property('XLim',xlim_target,'XTick',[],'YTick',[],'XColor',[1 1 1],'YColor',[1 1 1]);
laterality_figure(2,2).axe_property('XLim',xlim_saccade,'XTick',[],'YTick',[],'XColor',[1 1 1],'YColor',[1 1 1]);
laterality_figure(2,1).set_names('X',[]); laterality_figure(2,1).set_names('Y','Firing rate (spk/sec)');
laterality_figure(2,2).set_names('X',[]);

laterality_figure(3,1).axe_property('XLim',xlim_target,'YLim',[0 40],'XColor',[1 1 1],'XTick',[]);
laterality_figure(3,2).axe_property('XLim',xlim_saccade,'YLim',[0 40],'XColor',[1 1 1],'XTick',[],'YTick',[],'YColor',[1 1 1]);
laterality_figure(3,1).set_names('X',[]); laterality_figure(3,1).set_names('Y','Firing rate (Z-score)');
laterality_figure(3,2).set_names('X',[]);

laterality_figure(4,1).axe_property('XLim',xlim_target,'YLim',[-5 20]);
laterality_figure(4,2).axe_property('XLim',xlim_saccade,'YLim',[-5 20],'YTick',[],'YColor',[1 1 1]);

laterality_figure(4,1).set_names('X','Time from Target (ms)');
laterality_figure(4,2).set_names('X','Time from Saccade (ms)');


% Setup layout >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
laterality_figure(1,1).set_layout_options('Position',[0.1 0.85 0.4 0.1],... 
    'legend',false,'margin_height',[0.00 0.00],'margin_width',[0.0 0.00],'redraw',false);
laterality_figure(1,2).set_layout_options('Position',[0.52 0.85 0.4 0.1],... 
    'legend',false,'margin_height',[0.00 0.00],'margin_width',[0.0 0.00],'redraw',false);
laterality_figure(2,1).set_layout_options('Position',[0.1 0.72 0.4 0.1],... 
    'legend',false,'margin_height',[0.00 0.00],'margin_width',[0.0 0.00],'redraw',false);
laterality_figure(2,2).set_layout_options('Position',[0.52 0.72 0.4 0.1],... 
    'legend',false,'margin_height',[0.00 0.00],'margin_width',[0.0 0.00],'redraw',false);
laterality_figure(3,1).set_layout_options('Position',[0.1 0.45 0.4 0.25],... 
    'legend',false,'margin_height',[0.00 0.00],'margin_width',[0.0 0.00],'redraw',false);
laterality_figure(3,2).set_layout_options('Position',[0.52 0.45 0.4 0.25],... 
    'legend',false,'margin_height',[0.00 0.00],'margin_width',[0.0 0.00],'redraw',false);
laterality_figure(4,1).set_layout_options('Position',[0.1 0.15 0.4 0.25],... 
    'legend',false,'margin_height',[0.00 0.00],'margin_width',[0.0 0.00],'redraw',false);
laterality_figure(4,2).set_layout_options('Position',[0.52 0.15 0.4 0.25],... 
    'legend',false,'margin_height',[0.00 0.00],'margin_width',[0.0 0.00],'redraw',false);

% Setup colormap & lines
contralateral_color = [243, 175, 37]./255;
ipsilateral_color = [92, 121, 133]./255;
laterality_figure(1,:).set_color_options('map',[ipsilateral_color;contralateral_color]);
laterality_figure(2,:).set_color_options('map',[ipsilateral_color;contralateral_color]);
laterality_figure(3,:).set_color_options('map',[ipsilateral_color;contralateral_color]);
laterality_figure(4,:).set_color_options('map',[ipsilateral_color;contralateral_color]);


figure('Renderer', 'painters', 'Position', [100 100 450 700]);
laterality_figure.draw;





%%

for neuron_i = 1:size(visual_lateral,1)
    
    if visual_lateral.LSI_onset(neuron_i) > 0 & visual_lateral.lateral_p_onset(neuron_i) < 0.05
        LSI_group{neuron_i,1} = '3_contralateral';
    elseif visual_lateral.LSI_onset(neuron_i) < 0 & visual_lateral.lateral_p_onset(neuron_i) < 0.05
        LSI_group{neuron_i,1} = '1_ipsi';
    elseif visual_lateral.lateral_p_onset(neuron_i) > 0.05
        LSI_group{neuron_i,1} = '2_nonsig';
    end
end

clear laterality_histogram
laterality_histogram(1,1)=gramm('x',visual_lateral.LSI_onset,'color',LSI_group);
laterality_histogram(1,1).stat_bin('edges',[-1:0.05:1],'geom','overlaid_bar')
laterality_histogram(1,1).axe_property('YLim',[0 20])
figure('Renderer', 'painters', 'Position', [100 100 400 300]);
laterality_histogram.draw





