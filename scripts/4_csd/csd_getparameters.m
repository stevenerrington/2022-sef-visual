% Load relevant behavioral data
if ispc()
    driveDir = 'C:\Users\Steven\Desktop\Projects\';
    rootDir = 'C:\Users\Steven\Desktop\Projects\2022-sef-visual\';
    dataDir = fullfile(rootDir,'data');    
end



load(fullfile(dataDir, 'behavior','bayesianSSRT'));  % Bayesian SSRT estimates
load(fullfile(dataDir, 'behavior', 'executiveBeh'));  % Behavioural data
load(fullfile(dataDir, 'behavior', 'FileNames'));    % Filenames

%% Map LFP channel to session
sessionLFPmap = table();
channelNames = {'AD17';'AD18';'AD19';'AD20';'AD21';'AD22';'AD23';'AD24';...
    'AD25';'AD26';'AD27';'AD28';'AD29';'AD30';'AD31';'AD32';'AD33';...
    'AD34';'AD35';'AD36';'AD37';'AD38';'AD39';'AD40'};

site_sessionMap = {14:19,20:25,26:29,1:7,8:13};
site_labels = {"p1","p2","p3","np1","np2"};

DepthInfo = [1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1,...
    1, 1, 3, 2, 2, 3, 2, 4, 3, 4, 5, 4, 4, 6, 6];

LFPRange = [5 24; 5 24; 5 24; 5 24; 5 24; 5 24; 5 24; 5 24; 5 24; 5 24; 5 24; 5 24; 5 24;
6 22; 6 22; 8 24; 9 24; 12 24; 6 22; 13 24; 6 22; 8 24; 12 24; 12 24; 10 24;
7 22; 7 22; 8 24; 6 21]; % these are based on Godlove, 2014 figure, 5/6.


for sessionIdx = 1:29
    clear channelN session sessionName cortexFlag laminarFlag monkeyFlag monkeyName tempTable depth site sitelabel
    channelN = [1:24]';
    session = repmat(sessionIdx,24,1);  
    sessionName = repmat(FileNames(sessionIdx),24,1);
    
    if ismember(sessionIdx,site_sessionMap{1}); siteIdx = 1; end
    if ismember(sessionIdx,site_sessionMap{2}); siteIdx = 2; end
    if ismember(sessionIdx,site_sessionMap{3}); siteIdx = 3; end
    if ismember(sessionIdx,site_sessionMap{4}); siteIdx = 4; end
    if ismember(sessionIdx,site_sessionMap{5}); siteIdx = 5; end
        
        
    cortexFlag = zeros(24,1);
    cortexFlag(LFPRange(sessionIdx,1):LFPRange(sessionIdx,2)) = 1;
   
    laminarFlag = repmat(double(sessionIdx > 13),24,1);
    monkeyFlag =  repmat(double(ismember(sessionIdx,...
        executiveBeh.nhpSessions.XSessions)),24,1);
    monkeyName =  repmat(executiveBeh.nhpSessions.monkeyNameLabel(sessionIdx),24,1);
    
    site =  repmat(siteIdx,24,1);
    sitelabel = repmat(site_labels{siteIdx},24,1);
    
    if laminarFlag(1) == 0
        depth = nan(24,1);
    else
        depthCount = 0; %%%% THIS IS VOLATILE - CHECK ALIGNMENT. Originally 0
        for ii = 1:24
            if cortexFlag(ii) == 0
                depth(ii,1) = nan;
            else
               depthCount = depthCount + 1;
               depth(ii,1) = depthCount;
            end
            
        end
        
    end
        
    tempTable = table(channelN, session, sessionName, channelNames,...
        cortexFlag, site, sitelabel, laminarFlag, depth, monkeyFlag, monkeyName);

    sessionLFPmap = [sessionLFPmap; tempTable];
end

clear corticalLFPcontacts
corticalLFPcontacts.all = find(sessionLFPmap.cortexFlag == 1);
corticalLFPcontacts.subset.eu = find(sessionLFPmap.monkeyFlag(corticalLFPcontacts.all) == 0);
corticalLFPcontacts.subset.x = find(sessionLFPmap.monkeyFlag(corticalLFPcontacts.all) == 1);

corticalLFPcontacts.subset.laminar.all = find(sessionLFPmap.laminarFlag(corticalLFPcontacts.all) == 1);
corticalLFPcontacts.subset.laminar.upper = find(sessionLFPmap.laminarFlag(corticalLFPcontacts.all) == 1 & sessionLFPmap.depth(corticalLFPcontacts.all) < 8.5);
corticalLFPcontacts.subset.laminar.lower = find(sessionLFPmap.laminarFlag(corticalLFPcontacts.all) == 1 & sessionLFPmap.depth(corticalLFPcontacts.all) > 8.5);

corticalLFPcontacts.subset.laminar.eu = find(sessionLFPmap.laminarFlag(corticalLFPcontacts.all) == 1 &...
    sessionLFPmap.monkeyFlag(corticalLFPcontacts.all) == 0);
corticalLFPcontacts.subset.laminar.x = find(sessionLFPmap.laminarFlag(corticalLFPcontacts.all) == 1 &...
    sessionLFPmap.monkeyFlag(corticalLFPcontacts.all) == 1);

corticalLFPmap = sessionLFPmap(sessionLFPmap.cortexFlag == 1,:);

%% Analysis parameters:
%  Setup frequencies, channels, site power
csdParameters.samplingFreq = 1000;
csdParameters.frequencies = 15:29;
csdParameters.cycle = 7;

alignmentParameters.eventN = 2; % Align on target
alignmentParameters.alignWin = [-1000 2000];
alignmentParameters.time = alignmentParameters.alignWin(1):...
    alignmentParameters.alignWin(end)-1;

filterBands.all = [1 80];
filterBands.delta = [1 4];
filterBands.theta = [4 9];
filterBands.alpha = [9 15];
filterBands.beta = [15 30];
filterBands.lowGamma = [30 60];
filterBands.highGamma = [60 120];
filterBands.allGamma = [30 120];

filterNames = fieldnames(filterBands);

eventNames = {'fixate','target','stopSignal','saccade','sacc_end','tone','reward','sec_sacc'};

