%% Extract LFP
for eventType = [2,4]
    eventLabel = eventNames{eventType};
    alignmentParameters.eventN = eventType;
    fprintf(['Analysing data aligned on ' eventLabel '. \n']);
    ephysParameters.samplingFreq = 1000;
    
    %% Extract EEG data & calculate power
    for session = 14:29
        % Get session name (to load in relevant file)
        sessionName = FileNames{session};
        
        % Clear workspace
%         clear trials eventTimes inputLFP cleanLFP alignedLFP filteredLFP betaOutput morletLFP pTrl_burst
        
        % Setup key behavior variables
        ssrt = bayesianSSRT.ssrt_mean(session);
        eventTimes = executiveBeh.TrialEventTimes_Overall{session};

        % Load the LFP channels recorded for that session
        fprintf('Analysing session %d of %d... \n', session, 29);
        inputLFP = load(['T:\Users\Steven\dataRepo\2012_Cmand_EuX\rawData\' sessionName],...
            'AD1*','AD2*','AD3*','AD4*');
        
        lfpChannels = fieldnames(inputLFP);
      
        parfor k=1:numel(lfpChannels)
            if (isnumeric(inputLFP.(lfpChannels{k})))
                fprintf('......analysing LFP %d of %d... \n', k, numel(lfpChannels));
                
                LFP_aligned = [];
                
                % Pre-process & filter analog data (EEG/LFP), and align on event
                filter = 'all';
                filterFreq = filterBands.(filter);
                [~, LFP_aligned] = csd_getalignedLFP(inputLFP.(lfpChannels{k}), ephysParameters, filterFreq,...
                    eventTimes, alignmentParameters);

                savename_LFP = ['lfp_session' int2str(session) '_' lfpChannels{k} '_' eventLabel];
                
                parsave_filtered(fullfile(dataDir,'lfp', savename_LFP), LFP_aligned)          
                
            end
        end
    end
end

