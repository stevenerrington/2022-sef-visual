function  [PerTrialLFP] = SEF_stoppingEEG_trialAlignment(LFP, TrialSubset, TrialEventTimes_all, AlignWith, Win)

% The SDF here can either be spike density function OR LFPs! doesn't
% matter!! 

Win = [Win(1)+1 Win(2)];
PerTrialLFP = zeros(length(TrialSubset),range(Win)+1);

    for i = 1:length(TrialSubset)
        trial_num = TrialSubset(i);
        LowLimit = TrialEventTimes_all(trial_num,AlignWith)+Win(1); % -(length(R)/2);      Note: LowLimit is extended beyond the input dispWin by half of the length of R used for convolution.
        HighLimit = TrialEventTimes_all(trial_num,AlignWith)+Win(2); %+(length(R)/2);      Note: High limit is extended beyond the input dispWin by half of the length of R used for convolution.
        
        if isnan(LowLimit) == 0  && isnan(HighLimit) == 0  && LowLimit > 0  && HighLimit < length(LFP)
             PerTrialLFP(i,:) =  LFP(LowLimit:HighLimit);
        else
             PerTrialLFP(i,:) =  nan(1, range(Win)+1); 
        end
        
    end
    
    PerTrialLFP = PerTrialLFP * 1000;
    
end

    