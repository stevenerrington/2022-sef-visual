function CSD = getCSD(LFP2d, PadType, ChDist)

LFP = LFP2d;
% ChDist = 150;    % spacing between channels in microns.

% the formula is:  CSD = ( LFP(z+h, t) - 2*LFP(z,t) + LFP(z-h, t)  )/ (h^2)
CSD = single(nan( size(LFP2d,1)-2, size(LFP2d,2) ));
if nargin <=2
    if  nargin == 1
        ChDist = 150;    % spacing between channels in microns.
    end
    
    
    for i = 2:(size(LFP2d,1)-1)
        CSD(i-1,:) = (LFP(i+1, :) - 2*LFP(i,:) + LFP(i-1,:)) / (ChDist^2) ;  % Chen, Dhamala, Bollimunta and Ding; book chapter,
    end
    
    
    
elseif nargin == 3
    
    for i = 2:(size(LFP2d,1)-1)
        CSD(i-1,:) = (LFP(i+1, :) - 2*LFP(i,:) + LFP(i-1,:)) / (ChDist^2) ;  % Chen, Dhamala, Bollimunta and Ding; book chapter,
    end
    
    if strcmpi(PadType, 'Zeros') == 1 || strcmpi(PadType, 'Zero') == 1
        CSD = [zeros(1,size(CSD,2)); CSD; zeros(1,size(CSD,2))];
    elseif strcmpi(PadType, 'Duplicate') == 1  || strcmpi(PadType, 'Replicate') || strcmpi(PadType, 'Repeat') == 1
        CSD = [CSD(1,:); CSD; CSD(end,:)];
    end
    
end

CSD = -CSD .* 0.4;   % this converts units to the appropriate unit --> nA/mm^3





