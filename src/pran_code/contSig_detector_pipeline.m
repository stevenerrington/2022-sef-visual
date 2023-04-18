function contSig_patch = contSig_detector_pipeline(statsMatrix, continuityThreshold, sigtestType)
% statsMatrix contains 0's and 1's with 1's being significant.
% figure; imagesc(statsMatrix)
% Threshold_t_d = [50 2];

% sigtestType = 'strict' means that in both directions the patch has to be follow continuity threshold.
% sigtestType = 'relaxed' means that as long as the extreme ends of patch meet the continuity threshold, it's considered significant patch.

if nargin < 4
    fillPixel_flag = 0;
end

input = ~statsMatrix;
t_Thresh = continuityThreshold(1);
d_Thresh = continuityThreshold(2);

%%
if strcmpi(sigtestType,'strict') == 1
    contSig_statsMatrix = zeros( size(statsMatrix) );
    for currRow = 1:(size(input,1)-(d_Thresh-1))
        for currCol = 1:(size(input,2)-(t_Thresh-1))
            if sum(sum( input(currRow + [0:(d_Thresh-1)] , currCol + [0:(t_Thresh-1)]) ) ) == 0
                contSig_statsMatrix(currRow,currCol) = 1;
            else
                contSig_statsMatrix(currRow,currCol) = 0;
            end
        end
    end
    contSig_statsMatrix_full = contSig_statsMatrix;
    for currRow = 1:(size(input,1)-(d_Thresh-1))
        for currCol = 1:(size(input,2)-(t_Thresh-1))
            if contSig_statsMatrix(currRow,currCol) == 1
                contSig_statsMatrix_full(currRow + [0:(d_Thresh-1)] , currCol + [0:(t_Thresh-1)]) = 1;
            end
        end
    end
    contSig_patch = contSig_statsMatrix_full;
    
elseif strcmpi(sigtestType,'relaxed') == 1
    % if the criterion is relaxed, then we don't deal with the pixel
    % refinements here. We will deal with them in the patchEdgeDetector_pipeline.m code.
    contSig_patch = statsMatrix;
end
