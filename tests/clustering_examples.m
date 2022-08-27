
%% %% Analysis; clustering visual neurons
% 
% timeWin = [-250:500];
% timeWin_inputSDF = timewins.zero + timeWin;
% inputSDF = {sdf_visual_raw(:,timeWin_inputSDF), sdf_visual_raw(:,timeWin_inputSDF)};
% 
% sdfTimes = {timeWin, timeWin}; sdfEpoch = {timeWin, timeWin};
% 
% 
% [sortIDs,idxDist, raw, respSumStruct, rawLink,myK] =...
%     consensusCluster(inputSDF,sdfTimes,'-e',sdfEpoch);
% 
% 
% 
% normResp = scaleResp(inputSDF,sdfTimes,'max');
% 
% nClusters_manual = 3; clusterNeurons = [];
% for i = 1:nClusters_manual
%     clusterNeurons{i} = find(sortIDs(:,nClusters_manual) == i );
% end
% 
% % Plot clustering output %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % Dendrogram
% 
% figure('Renderer', 'painters', 'Position', [100 100 500 400]);
% 
% subplot(1,5,5);
% for ir = 1:size(raw,1)
%     for ic = (ir+1):size(raw,2)
%         raw(ic,ir) = raw(ir,ic);
%     end
% end
% [h,~,outPerm] = dendrogram(rawLink,0,'Orientation','right');
% set(gca,'YDir','Reverse');
% klDendroClustChange(h,rawLink,sortIDs(:,nClusters_manual))
% set(gca,'YTick',[]); xlabel('Similarity')
% subplot(1,5,[1:4]);
% 
% imagesc(raw(outPerm,outPerm));
% colormap(gray);
% xlabel('Unit Number'); set(gca,'YAxisLocation','Left');
% xticks([50:50:500]); yticks([50:50:500])
% 
% close all
% 
% 
% 










