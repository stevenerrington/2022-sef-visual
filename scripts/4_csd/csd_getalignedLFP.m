function [filterLFP, alignedLFP] = csd_getalignedLFP(inputData, ephysParameters,...
    filterFreq, eventTimes, alignmentParameters)

% Filter LFP
filterLFP = SEF_LFP_Filter(inputData,filterFreq(1), filterFreq(2), ephysParameters.samplingFreq);

% Align LFP on stop-signal
alignedLFP = SEF_stoppingEEG_trialAlignment(filterLFP,...
    1:size(eventTimes,1),...
    eventTimes, alignmentParameters.eventN, alignmentParameters.alignWin);

end
