function [OLdist1, OLdist2,  RTmatchedSubtrials_1_noIndexMatch,  RTmatchedSubtrials_2_noIndexMatch] = rtmatchTrials_laterality(Array1, Array2, BinSize)

% This function takes in two arrays. Bins each array and finds the
% intersection of the two arrays, and index of chosen elements from the
% original array:
% 

% The inputs: Array1 and Array2 are the two different distributions.
% Each has to have the following structure: 
%  it should be a N by 2 matrix.  values in the first column are trial numbers and values in the second column are the value that wants to be equated between the two arrays (say RT for instance).

% BinSize indicates the resolution by which we want this "equating" process
% to happen. For RTs I often use 10. So 181ms is treated the same as 189 as
% they are both binned the same if the bin size = 10.

% I highly recommend using BinSize = 10. Test the program for other
% BinSizes and make sure it's operating properly.

% embedded functions: 
% Index_1_2 = MatchTrialIndex(Indexset1, Indexset2)   
      % Imagine the case where in array1 we have 20 trials with an RT of
      % 180-190ms, and in array2 we have only 10 trials. To equate for RT,
      % we need to only select 10 trials from array1 in this RT range so that
      % from both arrays we only have 10 trials. But what 10 trials do we
      % pick? This MatchTrialIndex function keeps those 10 trials from array 1 that have
      % the closest index to trials in array2. This is important for
      % neurophysiology experiments where the isolation of units changes
      % over time.  
      % the reason for this is: imagine a senario where only 1 trial has to
      % be kept out of 2 trials (so one has to be discarded), but they are
      % separated in time by 15 minutes! Which one to keep? We keep the one
      % that's closest in time to the trial in the other array.

           
SubTrials_1 = Array1(:,1);
SubTrials_2 = Array2(:,1);

dist1(:,1) = Array1(:,2);
dist2(:,1) = Array2(:,2);

dist1(:,2) = ceil(dist1(:,1)/BinSize);  % each value on dist(:,2) determines one of the increments used for bar graph or histogram.
dist2(:,2) = ceil(dist2(:,1)/BinSize);

if isempty(dist1) == 0  && isempty(dist2) == 0

for i1 = 1:max(dist1(:,2))           % for all the different X values, from 1 all the way to max value:
    BarArray_1(1,i1) = i1;                           % the X value for the bar is: of course i1
    BarArray_1(2,i1) = sum(dist1(:,2) == i1);        % the y value for the bar is the count of entries with value = i1
end

for i2 = 1:max(dist2(:,2))           % same thing for the second distribition
    BarArray_2(1,i2) = i2;
    BarArray_2(2,i2) = sum(dist2(:,2) == i2);
end


% so far, we have just created the bar graphs for the two distributions. 
% Now we want to pick the portion of the histograms that are overlapping:

for Bin =  1 : min( max(BarArray_2(1,:)),  max(BarArray_1(1,:))  )   % For all the x values that are shared between the two bar graphs:
    BarArray_overlap(1,Bin) = Bin;                                   %  here is the X value
    BarArray_overlap(2,Bin) = min(BarArray_2(2,Bin),BarArray_1(2,Bin));   %  The Y value is the lower of the y-values from the two distributions. That's the part that's overlapping.
end

% Now that we have the overlapping distribution of values, we need to
% extract those trials that match this part of the distribution:

RTmatchedSubtrials_1 = [];
RTmatchedSubtrials_2 = [];
RTmatchedSubtrials_1_noIndexMatch = [];
RTmatchedSubtrials_2_noIndexMatch = [];

RTindex = find(BarArray_overlap(2,:)> 0);        % These are the x values with non-zero y values. This is because our range was from 1 to XXX... 

OLdist1 = [];    % initializing the two lists for the two distributions:
OLdist2 = [];

for p = 1:length(RTindex)      % for all the non-zero y values on the over-lapping distribution do the following:
    S1 = [];               
    S2 = []; 
    IndexSet1=[];
    IndexSet2=[];
    Index_1_2 = [];
     
    i = RTindex(p);
    
    S1 = find(dist1(:,2) == i);      % these are the indeces for the ones with a given value used for pairing
    S2 = find(dist2(:,2) == i);
    
    % Here I can define IndexSet1 and IndexSet2 for Right side, separately
    % from that for Left side!!    
    
    
    IndexSet1= SubTrials_1(S1);      % this extracts the trial indeces
    IndexSet2= SubTrials_2(S2);
    
    Index_1_2 = MatchTrialIndex (IndexSet1, IndexSet2);     % This is the function I wrote which matches indexes to closest so that the overall difference in index # between the two sets is minimized.
                                                  
    OLindex1 = Index_1_2(:,1);  % Trial numbers (from the main array) for overlapped distribution, taken from array 1
    OLindex2 = Index_1_2(:,2);  % Trial numbers (from the main array) for overlapped distribution, taken from array 2
      
    RTmatchedSubtrials_1 = [RTmatchedSubtrials_1; Index_1_2(:,1)];     % Doing this for all the different bin values
    RTmatchedSubtrials_2 = [RTmatchedSubtrials_2; Index_1_2(:,2)];
    
    RTmatchedSubtrials_1_noIndexMatch{p} = IndexSet1;
    RTmatchedSubtrials_2_noIndexMatch{p} = IndexSet2;
    
end
    % So at this point we have the trial numbers. Let's now package it up:

    % we want to have the trial numbers in one column, and RTs in 2nd
    % column
D1 = RTmatchedSubtrials_1;
if isempty(RTmatchedSubtrials_1) == 0  && isempty(Array1) == 0 && isempty(RTmatchedSubtrials_2) == 0  && isempty(Array2) == 0
    Fr = unique(RTmatchedSubtrials_1);
    for er = 1:length(Fr)
        if length( find(RTmatchedSubtrials_1 == Fr(er) )) ~=  length(   find(Array1(:,1) == Fr(er) ))
           
            RTmatchedSubtrials_1(D1 == Fr(er) ) =  [];
            RTmatchedSubtrials_2(D1 == Fr(er) ) =  [];
        end
    end
    
else 
    OLdist1 = [];
end

D2 =  RTmatchedSubtrials_2;
if isempty(RTmatchedSubtrials_1) == 0  && isempty(Array1) == 0 && isempty(RTmatchedSubtrials_2) == 0  && isempty(Array2) == 0
     Fr2 = unique(RTmatchedSubtrials_2);
    for er2 = 1:length(Fr2)
        if length( find(RTmatchedSubtrials_2 == Fr2(er2) )) ~=  length(   find(Array2(:,1) == Fr2(er2) ))
            
            RTmatchedSubtrials_2(D2 == Fr2(er2) ) =  [];
            RTmatchedSubtrials_1(D2 == Fr2(er2) ) =  [];
        end
    end
    
OLdist1(:,1) = RTmatchedSubtrials_1;
OLdist1(:,2) = Array1( ismember(Array1(:,1),RTmatchedSubtrials_1) == 1 , 2 ) ;
OLdist1 = sortrows(OLdist1,1);

OLdist2(:,1) = RTmatchedSubtrials_2;
OLdist2(:,2) = Array2( ismember(Array2(:,1),RTmatchedSubtrials_2) == 1 , 2 ) ;
OLdist2 = sortrows(OLdist2,1);
else
    OLdist2 = [];
end

else
    OLdist1 = [];
    OLdist2 = [];
end

%%
function Index_1_2 = MatchTrialIndex(Indexset1, Indexset2)


% A = randperm(1000)
% Indexset1 = sort(A(1:10)')
% Indexset2 = sort(A(200:220)')

if length(Indexset1)<= length(Indexset2)
    Set1 = Indexset1;
    Set2 = Indexset2;
else
    Set1 = Indexset2;
    Set2 = Indexset1;
end


for i = 1:length(Set1)   % Set1 is the smaller set.
    minDiff = min(abs(Set2 - Set1(i)));
    Index_minDiff = find( abs(Set2 - Set1(i)) == minDiff);
    if length(Index_minDiff) > 1
        Index_minDiff = Index_minDiff(1);
    end
    Index_1_2_unSorted(i,:) = [Set1(i) Set2(Index_minDiff)];
    Set2(Index_minDiff) = [];
end



if length(Indexset1)<= length(Indexset2)
Index_1_2 = [sort(Index_1_2_unSorted(:,1))  sort(Index_1_2_unSorted(:,2))];
else
    Index_1_2 = [sort(Index_1_2_unSorted(:,2))  sort(Index_1_2_unSorted(:,1))];
end


