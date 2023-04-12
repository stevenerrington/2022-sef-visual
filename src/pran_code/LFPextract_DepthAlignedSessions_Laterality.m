clear ; close all;
load('FileNames.mat');
load('executiveBeh.mat');
load('ChannelDepthMap.mat');
load('piadepth_gray.mat')
figure;
Windowtgt=[-200 200];
Windowbl= [-200 0]; %baseline window includes pretarget ramping?
cd('C:\Users\pran_\OneDrive - York University\Desktop (Razer)\2022-sef-visual-main\mat_data'); %l
% cd('C:\Users\thirunap\OneDrive - York University\Desktop (Razer)\2022-sef-visual-main\mat_data');
for j=1:2
    for k=1:6 %perpendicular recording sessions 1:6 euler; 7:16 xena
        load(FileNames{k+13})
        LFPdata = [AD17 AD18 AD19 AD20 AD21 AD22 AD23 AD24 AD25 AD26 AD27 AD28 AD29 AD30 AD31 AD32 AD33 AD34 AD35 AD36 AD37 AD38 AD39 AD40]/1000;
        %     LFPdata = horzcat(zeros(length(LFPdata),max(cellfun(@min,ChannelDepthMap.GrayMatterMapping))...
        %         -min(ChannelDepthMap.GrayMatterMapping{k+13,1})-1),...
        %         LFPdata,...
        %         zeros(length(LFPdata),min(ChannelDepthMap.GrayMatterMapping{k+13,1})...
        %         -max((cellfun(@min,ChannelDepthMap.GrayMatterMapping))-min(cellfun(@min,ChannelDepthMap.GrayMatterMapping)))+1));
        DepthAlignment = [min(ChannelDepthMap.GrayMatterMapping{k+13,1}) max(ChannelDepthMap.GrayMatterMapping{k+13,1})];
        LFPdata = LFPdata(:,DepthAlignment(1):DepthAlignment(2));
        %     LFPdata = horzcat(LFPdata, zeros(length(LFPdata),max(cellfun(@length,ChannelDepthMap.GrayMatterMapping))-length(ChannelDepthMap.GrayMatterMapping{k+13,1})));
        size(LFPdata,2)
        TgtAlign = Target_ + TrialStart_;

%         if j==1
%             TgtAlign = TgtAlign(executiveBeh.ttx.GO_Left{k+13},:);
%         else
%             TgtAlign = TgtAlign(executiveBeh.ttx.GO_Right{k+13},:);
%         end
        TgtAlignALL = TgtAlign(executiveBeh.ttx.GO{k+13},:);

        for i=1:length(TgtAlign)
            LFPtgtarray = LFPdata(TgtAlign(i)+Windowtgt(1):TgtAlign(i)+Windowtgt(2)-1,:);
            LFPblarray = LFPdata(TgtAlign(i)+Windowbl(1):TgtAlign(i)+Windowbl(2)-1,:);
            LFPtgtArray(:,:,i)=LFPtgtarray-mean(LFPblarray);
        end

        avgtgtLFP=transpose(mean(LFPtgtArray,3));

        %     clearvars -except avgtgtLFP FileNames executiveBeh Window

        %     plot(transpose(avgtgtLFP));

        %% Parameters %need to scale so they have the same length...
        Ne= size(avgtgtLFP,1);
        %Ne = diff(DepthAlignment)+1; % number of electrodes in the shank
        a = 0.1; % [mm] position of the first electrode, shifts sources down, assumed 1 micron
        elec_spacing = 0.15; % [mm] electrode spacing
        ze = a:elec_spacing:((Ne-1)*elec_spacing + a); % electrode positions wit
        % respect to the pia surface
        % WARNING: position of the first electrode must be different from zero
        el_pos = ze*1e-3;  % mm to m
        cond = 0.4; %[S/m] gray matter conductance | NOTE: if length(cond)==1, the
        % function iCSDspline_Pettersen() considers con_top = cond (conductance at
        % the top of the cylinder or above the pia matter) %papers with macaques,
        % 0.4
        gauss_sigma = 0.1e-3;   %[m] Gaussian filter std, how smooth you want LFP x depth, increase more smooth
        diam = 3e-3; % [m] cylinder diameter
        Ve = avgtgtLFP;
        %%
        % Ve == LOCAL FIELD POTENTIALS IN VOLTS
        % Outputs: zs ~ mm & iCSD ~ nA/mm3
        [zs, iCSD] = iCSDspline_Pettersen(Ve,el_pos,diam,cond,gauss_sigma);
        %%
        iCSD = vertcat(iCSD,NaN(10*max(cellfun(@length,ChannelDepthMap.GrayMatterMapping)-length(ChannelDepthMap.GrayMatterMapping{k+13,1})),size(iCSD,2)));
        tspan = -200:200;
        font = 9;

        figure(1)
        if j==1
         subplot(2,6,k)
        else
         subplot(2,6,k+6)
        end

%         if j==1
%             subplot(2,10,k-6)
%         else
%             subplot(2,10,k+4)
%         end
        % subplot(1,6,k)
        % plot([-200:1:199],transpose(avgtgtLFP));
        % xlim([-200 200])
        % ylim([-5e-5 1e-5])
        % subplot(1,2,2)
        imagesc(tspan, piadepth_gray, iCSD);
        c = colorbar;
        colormap(flip(jet));
        % c.Label.String = {'CSD', 'nA/mm^{3}'};
        bar_max = 200; %Eu
%         bar_max = 130; %Xe
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
        xticks([-50 0 50 100 150 200])
        xticklabels({'-50', '0', '50', '100', '150', '200'});
        xlim([-50 200])
        ylabel('Depth (mm)')
        xlabel('Time from target array onset')
        set(ax,'fontweight','bold','FontSize',9,'LineWidth',2)
        xline(0,'LineWidth',3);
        set(gcf, 'Position', get(0, 'Screensize'));
        colorbar off
        % set(gca,'ytick',[])
        hold on
        subtitle(sprintf('%s',FileNames{k+13}));

        if j==1
            perpCSD.Left.session(k) = k+13;
            perpCSD.Left.CSD(:,:,k) = iCSD;
        else
            perpCSD.Right.session(k) = k+13;
            perpCSD.Right.CSD(:,:,k) = iCSD;
        end
        % saveas(gcf,fullfile('C:\Users\pran_\OneDrive - York University\Desktop (Razer)\2022-sef-visual-main\CSDplots', [FileNames{k+13} '.png']));
        % % close all
        clearvars -except FileNames executiveBeh ChannelDepthMap piadepth_gray Windowtgt Windowbl perpCSD LFPdata Fcs j k
    end

    % close all
    figure(2)
    subplot(1,2,j)
    tspan = -200:200;
    font = 9;
    if j==1
        GrandAveragedCSD=nanmean(perpCSD.Left.CSD,3);
    else
        GrandAveragedCSD=nanmean(perpCSD.Right.CSD,3);
    end
    imagesc(tspan, piadepth_gray, GrandAveragedCSD);
    c = colorbar;
    colormap(flip(jet));
    % c.Label.String = {'CSD', 'nA/mm^{3}'};
    bar_max = 110; %Eu
%     bar_max = 40; %Xe
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
    % axis square
    clearvars -except FileNames executiveBeh ChannelDepthMap piadepth_gray Windowtgt Windowbl j k perpCSD
end
%%
close all
    figure (1)
    tspan = -200:200;
    font = 9;
%   CSDdiff=(nanmean(perpCSD.Left.CSD,3)-nanmean(perpCSD.Right.CSD,3))./(nanmean(perpCSD.Left.CSD,3)+nanmean(perpCSD.Right.CSD,3));
    CSDdiff=(nanmean(perpCSD.Left.CSD,3)-nanmean(perpCSD.Right.CSD,3));
%     CSDdiff(CSDdiff>1)=0;
%     CSDdiff(CSDdiff<-1)=0;
    CSDdiff=smoothdata(CSDdiff,'gaussian',10);
    
    imagesc(tspan, piadepth_gray, CSDdiff);
    c = colorbar;
    colormap(flip(jet));
    % c.Label.String = {'CSD', 'nA/mm^{3}'};
    % bar_max = 110; %Eu
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

    
