function [TT1_RT_LR_match,  TT2_RT_LR_match, TT1, TT2, AllMatchedTrials] = rtmatchTrials_value(TT1, TT2, TrialEventTimes_all, Infos_)


RT = TrialEventTimes_all(:,4) - TrialEventTimes_all(:,2);

%%
TT1_Left = TT1(Infos_.Target_angle(TT1) == 180);
TT1_Right = TT1(Infos_.Target_angle(TT1) == 0);

TT2_Left = TT2(Infos_.Target_angle(TT2) == 180);
TT2_Right = TT2(Infos_.Target_angle(TT2) == 0);


%%
TT1_Left = [TT1_Left  RT(TT1_Left)];
TT1_Right = [TT1_Right  RT(TT1_Right)];
TT2_Left = [TT2_Left  RT(TT2_Left)];
TT2_Right = [TT2_Right  RT(TT2_Right)];


%%
BinSize = 10;

[TT1_Left_Matched, TT2_Left_Matched, AllMatchedTrials.L1, AllMatchedTrials.L2] = rtmatchTrials_laterality(TT1_Left, TT2_Left, BinSize);
[TT1_Right_Matched, TT2_Right_Matched, AllMatchedTrials.R1, AllMatchedTrials.R2] = rtmatchTrials_laterality(TT1_Right, TT2_Right, BinSize);

%%
TT1_Matched = [TT1_Left_Matched; TT1_Right_Matched];
TT2_Matched = [TT2_Left_Matched; TT2_Right_Matched];

%%
TT1_RT_LR_match = TT1_Matched(:,1);
TT2_RT_LR_match = TT2_Matched(:,1);


