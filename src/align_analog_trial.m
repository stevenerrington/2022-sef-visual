function  [trial_signal] = align_analog_trial(signal, TrialEventTimes_all, AlignWith, Win)

% The SDF here can either be spike density function OR LFPs! doesn't
% matter!! 

Win = [Win(1)+1 Win(2)];
trial_signal = zeros(length(TrialEventTimes_all),range(Win)+1);

    for i = 1:length(TrialEventTimes_all)
        LowLimit = TrialEventTimes_all(i,AlignWith)+Win(1); % -(length(R)/2);      Note: LowLimit is extended beyond the input dispWin by half of the length of R used for convolution.
        HighLimit = TrialEventTimes_all(i,AlignWith)+Win(2); %+(length(R)/2);      Note: High limit is extended beyond the input dispWin by half of the length of R used for convolution.
        
        if isnan(LowLimit) == 0  && isnan(HighLimit) == 0  && LowLimit > 0  && HighLimit < length(signal)
             trial_signal(i,:) =  signal(LowLimit:HighLimit);
        else
             trial_signal(i,:) =  nan(1, range(Win)+1); 
        end
        
    end
        
end

    