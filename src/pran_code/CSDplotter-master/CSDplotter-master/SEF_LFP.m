%Trial Alignment Times at Target Presentation
TargetAlign=TrialStart_+Target_;
TargetAlign=TargetAlign(~isnan(TargetAlign));
Chann_=[AD17 AD18 AD19 AD20 AD21 AD22 AD23 AD24 AD25 AD26 AD27 AD28 AD29 AD30 AD31 AD32 AD33 AD34 AD35 AD36 AD37 AD38 AD39 AD40]; 

%Windowed
win1=-200;
win2=300;

LFPARRAY=[];

for j=1:size(Chann_,2)
    Contact=Chann_(:,j);
    LFPmult=[];
for i=1:length(TargetAlign)
    LFPsing=Contact(TargetAlign(i)+win1:TargetAlign(i)+win2-1);
    LFPmult=[LFPmult LFPsing];
end
    LFPmultmean=mean(LFParray,2)
    LFPARRAY=[LFPARRAY LFPmultmean];
end
LFPARRAY=transpose(LFPARRAY)

figure(1)
for i=1:size(LFPARRAY,1)
    plot(LFPARRAY(i,:))
    hold on
end
