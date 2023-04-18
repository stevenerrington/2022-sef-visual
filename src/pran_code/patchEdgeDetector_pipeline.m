function [patchEdges, patch] = patchEdgeDetector_pipeline(patchPixels, threshold_t_d, sigtestType)

% THIS IS CODE TAKEN FROM Codes_Jan2017\CSD-related Codes\CSD_newStats
figure;
[C , ~] = contour(patchPixels, 1);    % C gives us the pixels of significant value.
close;
% Now let's separate out each patch.
% Looking at the C matrix, it appears that whenever a patch ends and a new
% one begins, there is a sharp change in C matrix.
% looking at rate of change of C(1,:) will do it:
RtC = C(1, 2:end) - C(1, 1:end-1);
DtC = C(2, 2:end) - C(2, 1:end-1);
InfTime = find( abs(RtC) > 1  | abs(DtC) > 1 );
if mod(length(InfTime), 2) == 1
    InfTime = [InfTime  size(C,2)];
end
patch(1,:) = 1:length(InfTime)/2;
patch(2,:) = InfTime(patch(1,:)*2 - 1) + 1;
patch(3,:) = InfTime(patch(1,:)*2 ) - 1;

for pc = 1: size(patch,2)
    patch(4,pc) = max( C(1, patch(2,pc):patch(3,pc)+1 )) - min( C(1, patch(2,pc):patch(3,pc)+1 ));
    patch(5,pc) = max( C(2, patch(2,pc):patch(3,pc)+1 )) - min( C(2, patch(2,pc):patch(3,pc)+1 ));
end
patch = int64(patch);
patchtag = 0;
patchEdges = [];

if isempty(patch) == 0
    
    if strcmpi(sigtestType,'strict') == 1
        
        for p = 1:size(patch,2)
            Allpixels = [];
            pixels = [];
            % The C matrix has bunch of XXX.5's in there. In matrix coortinates we
            % can't use 0.5's, so we Ceil.
            Allpixels(:,1) = ceil(C(2, patch(2,p):(patch(3,p)+1) ));
            Allpixels(:,2) = ceil(C(1, patch(2,p):(patch(3,p)+1) ));
            % But Ceil has the problem that it may take edge values over the edge. To
            % avoid this we remove those pixels that now fall off the range.
            pixels = unique(Allpixels,'rows');
            pixels( pixels(:,1) == 0 | pixels(:,1) > size(patchPixels,1) , :) = [];
            pixels( pixels(:,2) == 0 | pixels(:,2) > size(patchPixels,2) , :) = [];
            if size(pixels,1) > 1
                patchtag = patchtag + 1;
                patchEdges{patchtag}(1,:) = C(1, patch(2,p): patch(3,p)+1);
                patchEdges{patchtag}(2,:) = C(2, patch(2,p): patch(3,p)+1);
            end
            
        end
        
        
    elseif strcmpi(sigtestType,'relaxed') == 1
        
        TimeThreshold = threshold_t_d(1);
        DepthThreshold = threshold_t_d(2);
        
        for p = 1:size(patch,2)
            Allpixels = [];
            pixels = [];
                ThresholdTest = (patch(4,p) >= TimeThreshold  &&  patch(5,p) >= DepthThreshold);
            if ThresholdTest == 1  % let's set 50 as out threshold: at least 50 pixels to qualify as a meaningful patch.
                % The C matrix has bunch of XXX.5's in there. In matrix coortinates we
                % can't use 0.5's, so we Ceil.
                Allpixels(:,1) = ceil(C(2, patch(2,p):patch(3,p)+1 ));
                Allpixels(:,2) = ceil(C(1, patch(2,p):patch(3,p)+1 ));
                % But Ceil has the problem that it may take edge values over the edge. To
                % avoid this we remove those pixels that now fall off the range.
                pixels = unique(Allpixels,'rows');
                pixels( pixels(1,:) == 0 | pixels(1,:) > size(patchPixels,1) , :) = [];
                pixels( pixels(2,:) == 0 | pixels(2,:) > size(patchPixels,2) , :) = [];
                
                % Now we extract the data from CSD_all matrix, and plot histogram:
                %         for j = 1:size(pixels,1)
                %             StatsLine(pixels(j,1), pixels(j,2)) = 1;
                %         end
                if size(pixels,1) > 1
                    patchtag = patchtag + 1;
                    patchEdges{patchtag}(1,:) = C(1, patch(2,p): patch(3,p)+1);
                    patchEdges{patchtag}(2,:) = C(2, patch(2,p): patch(3,p)+1);
                end
            end
        end
    end
end
