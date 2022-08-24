function [SigModTime, ModTime] = OnsetLatency(SDF, Mean_BaseLine, SD_BaseLine, CP1, CP2, Cdur1)

% take in two inputs: SDF1 and SDF2 (time is not corrected).
% need to know the baseline average corresponding to SDF1 and SDF2.

% Input: Baseline mean, and Baseline SD. This will be used to create your
% cutoff points:

% 1st criteria: Mark all the points which the difference hits the Critical Point 1 mark:


% The output of this function is: 
% SigDiffTime: This is a matrix with 1st row showing timing of START of
% significant difference (significant means all criteria are met). 2nd row
% shows the timing of the END of significant difference

% DiffTime: This is a matrix similar to SigDiffTime, but shows ALL
% difference based on CP1 criteria, regardless of whether or not they meet
% the additional requirements defined by CP2 or Cdur1 (see below).
% 
% SDF1 = A;
% SDF2 = B;
% 
 Mean = Mean_BaseLine;
 SD = SD_BaseLine;
% 
% CP1 = 2;   % 2SD is the cut off.
% CP2 = 6;        % Critical Point 2: In the period where CP1 criteria is met, is the activity hitting  CP2 * SD?  if not, we are not accepting it.
% Cdur1 = 50;     % Critical Duration 1: how long should CP1 be met in order to consider it a change? If it's shorter than Cdur1 we are not accepting it.


DiffReMean = abs( (SDF - Mean));    % absolute value is taken - so only DIFFERENCE matters... doesn't matter what direction the difference is. Also, this value is taken away from the Mean. So, it's centered around the Mean.
TotTime = length(DiffReMean);     % This is the only time window that we care about. 

count_m1 = 0;              % initialize the count for mark 1 (hence count_m1) - Mark 1 is the beginning of difference (defined by CP1 criteria).
count_m2 = 0;              % initialize the count for mark 2 (hence count_m2) - mark 2 is the END of difference (defined by CP1 criteria).
DiffMark1 = [];
DiffMark2 = [];

for T = 2 : TotTime-1      % scanning through time. Checking each time point. I go from 2 to TotTime-1 just because I will use +/- 1 ms window to detect mark1 (start of difference) or mark2 (end of difference)

    if DiffReMean(T) >= CP1*SD  & DiffReMean(T-1) < CP1*SD      % marking when the transition from BELOW to ABOVE CP1 takes place.
        count_m1 = count_m1 + 1;     % count_m1 counts the number of times we are marking beginning of difference. Could remain 0 for the entire interval.
        DiffMark1(count_m1) = T;      %start of the difference > Critical Point 1
    end
    
    if DiffReMean(T) >= CP1*SD  &  DiffReMean(T+1) < CP1*SD        % marking when the transition from ABOVE to BELOW CP1 takes place.
        count_m2 = count_m2 + 1;     % count_m1 counts the number of times we are marking beginning of difference. Could remain 0 for the entire interval.
        DiffMark2(count_m2) = T;     % start of the difference > Critical Point 1
    end
    
end


if isempty(DiffMark1) == 1  & isempty(DiffMark2) == 0
      DiffMark1 = 1;
end
if isempty(DiffMark1) == 0  & isempty(DiffMark2) == 1
     DiffMark2 = length(DiffReMean);    
end
if isempty(DiffMark1) == 1  & isempty(DiffMark2) == 1 & DiffReMean(2) > CP1*SD & DiffReMean(TotTime-1) > CP1*SD &  length(find(DiffReMean(2:TotTime-1)>CP2*SD))>0
     DiffMark1 = 1;
     DiffMark2 = length(DiffReMean);    
end




if isempty(DiffMark1) == 0  &  isempty(DiffMark2) == 0
  
if DiffMark1(1) > DiffMark2(1)
    DiffMark1 = [1 DiffMark1];
end

     
Mark1Num = length(DiffMark1);
Mark2Num = length(DiffMark2);


if Mark1Num == Mark2Num+1
    DiffMark2(Mark1Num) = TotTime;   % The last time-step is taken as the end of diffMark2, of course, it's not really it's end. but we take it as the end since there isn't a real end. Could have been prolonged further.
elseif Mark1Num == Mark2Num - 1
    DiffMark1(2:Mark2Num) = DiffMark1;   % The last time-step is taken as the end of diffMark2, of course, it's not really it's end. but we take it as the end since there isn't a real end. Could have been prolonged further.
    DiffMark1(1) = 1;   % The last time-step is taken as the end of diffMark2, of course, it's not really it's end. but we take it as the end since there isn't a real end. Could have been prolonged further.
end




ModTime = [DiffMark1; DiffMark2];

% % % % % So now we have DiffMark1 and DiffMark2 (but these are based on CP1
% % % % % criteria (likely +/-2SD cutoff).

% % % % % We also have 2 more criteria to consider:
% 
% % % % % if length(DiffMark1) ~= length(DiffMark2)
% % % % % if DiffMark2(1) > DiffMark1(1) 
    

NumOfDiff = length(DiffMark1);
KeepDiff = [];

for DiffOfInterest = 1:NumOfDiff
    Range = [DiffMark1(DiffOfInterest) DiffMark2(DiffOfInterest)];
    DiffDur(DiffOfInterest) = Range(2)-Range(1);
    PeakDiff(DiffOfInterest) = max(DiffReMean(Range(1):Range(2)));
    if PeakDiff(DiffOfInterest) >= SD*CP2  &  DiffDur(DiffOfInterest) >= Cdur1
        KeepDiff(DiffOfInterest) = 1;          % 1 means keep this
    else
        KeepDiff(DiffOfInterest) = 0;          % 0 means difference is not significant. So should be discarded.
    end
end

SigDiffMark1 =   DiffMark1(find(KeepDiff == 1));
SigDiffMark2 =   DiffMark2(find(KeepDiff == 1));



else
    SigDiffMark1 = [];
    SigDiffMark2 = [];
    DiffMark1 = [];
    DiffMark2 = [];
    
    
    
end

SigModTime = [SigDiffMark1; SigDiffMark2];
ModTime = [DiffMark1; DiffMark2];




