
% Load in pre-processed data
load(fullfile(dirs.root,'src','pran_code','FileNames.mat'));
load(fullfile(dirs.root,'src','pran_code','executiveBeh.mat'));
load(fullfile(dirs.root,'src','pran_code','ChannelDepthMap.mat'));
load(fullfile(dirs.root,'src','pran_code','piadepth_gray.mat'));

Windowtgt=[-200 200];
Windowbl= [-200 0]; %baseline window includes pretarget ramping?

for perp_session_i=1:6 %perpendicular recording sessions 1:6 euler; 7:16 xena
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
    
    for trial_type_i=1:2

        % Get alignment times
        TgtAlign_in = [];
        TgtAlign_in = data_in.Target_ + data_in.TrialStart_;
        
        TgtAlign = [];
        if trial_type_i==1
            TgtAlign = TgtAlign_in(executiveBeh.ttx.GO_Left{perp_session_i+13},:);
        else
            TgtAlign = TgtAlign_in(executiveBeh.ttx.GO_Right{perp_session_i+13},:);
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


%%
figure (1)
tspan = -200:200;
font = 9;
CSDdiff=(nanmean(perpCSD.Left.CSD,3)-nanmean(perpCSD.Right.CSD,3));
CSDdiff=smoothdata(CSDdiff,'gaussian',10);

imagesc(tspan, piadepth_gray, CSDdiff);
c = colorbar;
colormap(flip(jet));

bar_max = 100; %Xe
bar_min = -bar_max;
caxis([bar_min bar_max]);
c.Label.FontSize = font+2;
c.FontSize = font+2;
c.FontWeight = 'bold';
c.Ticks = [bar_min 0 bar_max];
% c.TickLabels = {'-350', '', '+350'};
c.LineWidth = 1;
ax = gca;
ax.FontSize = font;
xticks([-200 -150 -100 -50 0 50 100 150 200])
xticklabels({'-200','-150','-100','-50', '0', '50', '100', '150', '200'});
xlim([-200 200])
ylabel('Depth (mm)')
xlabel('Time from target array onset')
set(ax,'fontweight','bold','FontSize',9,'LineWidth',2)
xline(0,'LineWidth',3);
set(gcf, 'Position', get(0, 'Screensize'));
xlim([-50 200])
axis square
