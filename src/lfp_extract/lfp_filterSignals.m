function LFP = lfp_filterSignals(tdtLFP,filter)

channelNames = fieldnames(tdtLFP.data);
nChannels = length(channelNames);


parfor channelIdx = 1:nChannels
    channel = channelNames{channelIdx};  
    filterChanData{channelIdx} = util_bandpassFilter(double(tdtLFP.data.(channel)),...
            1017.2526, 1, filter.band, filter.label, 2);
end


for channelIdx = 1:nChannels
    channel = channelNames{channelIdx};
    
    LFP.data.(channel) = ...
        filterChanData{channelIdx}.data;
end